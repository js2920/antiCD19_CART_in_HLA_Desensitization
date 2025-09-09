#!/usr/bin/env Rscript
# Build eplet_map.csv and save a Registry snapshot for reproducibility.
# Outputs are written to the specified directory (default below).
#
# Usage:
#   Rscript buildepletmap_A.R "/path/to/MFI.csv" "/home/eplet_ref" --verified-only


suppressPackageStartupMessages({
  # no-op, we’ll load packages after ensuring installation
})

# ───────────────────────────── Config / Args ─────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

default_mfi    <- "/home/Luminex_eplet_MFIs.csv" (supplementary data file of Manuscript)
default_outdir <- "/home/eplet_ref"

mfi_path <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_mfi
out_dir  <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_outdir
flags    <- if (length(args) >= 3) args[3:length(args)] else character(0)

verified_only <- any(grepl("^--verified-only$", flags))

cran_mirror <- "https://cloud.r-project.org"

# ───────────────────────────── Helpers ─────────────────────────────
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("[SETUP] Installing ", pkg, " from CRAN…")
    install.packages(pkg, repos = cran_mirror, dependencies = TRUE)
  }
}

ensure_hlapro <- function() {
  if (!requireNamespace("hlapro", quietly = TRUE)) {
    ensure_pkg("remotes")
    message("[SETUP] Installing 'hlapro' from GitHub (lcreteig/hlapro)…")
    remotes::install_github("lcreteig/hlapro")
  }
}

normalize_antigen <- function(x) {
  x <- gsub("\\s*\\([^)]*\\)\\s*$", "", x) # drop trailing "(A2)" etc.
  x <- gsub("\\s+", " ", trimws(x))
  x
}

detect_antigens <- function(df_names) {
  is_ag <- grepl("\\*", df_names) & !grepl("(?i)^(pc_|nc_)", df_names)
  ag_raw <- df_names[is_ag]
  if (!length(ag_raw)) stop("No antigen columns detected (looking for '*').")
  unique(normalize_antigen(ag_raw))
}

to_semicolon <- function(v) if (length(v) == 0) "" else paste(unique(v), collapse = ";")

write_manifest <- function(dir, mapping, params, pkgs) {
  fn <- file.path(dir, "MANIFEST.txt")
  con <- file(fn, "wt")
  on.exit(close(con), add = TRUE)
  cat("eplet_map build manifest\n", file = con)
  cat("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", file = con, sep = "")
  cat("MFI file: ", normalizePath(params$mfi_path, mustWork = FALSE), "\n", file = con, sep = "")
  cat("Verified-only: ", params$verified_only, "\n", file = con, sep = "")
  cat("Output dir: ", normalizePath(dir, mustWork = FALSE), "\n\n", file = con, sep = "")

  cat("Detected antigens (first 20):\n", file = con)
  cat(paste0("  - ", head(mapping$Antigen, 20)), sep = "\n", file = con)
  if (nrow(mapping) > 20) cat("  …\n", file = con)

  cat("\nR session info (key packages):\n", file = con)
  for (p in pkgs) {
    ver <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
    cat(sprintf("  %s: %s\n", p, ver), file = con)
  }
}

save_registry_snapshot <- function(reg, out_dir) {
  # Try to save both a CSV (if possible) and a robust RDS
  csv_path <- file.path(out_dir, sprintf("eplet_registry_snapshot_%s.csv",
                                         format(Sys.Date(), "%Y%m%d")))
  rds_path <- file.path(out_dir, sprintf("eplet_registry_snapshot_%s.rds",
                                         format(Sys.Date(), "%Y%m%d")))
  ok_csv <- FALSE
  if (inherits(reg, "data.frame")) {
    # Avoid gigantic wide tables by writing as-is; users can readRDS for full structure if needed.
    readr::write_csv(reg, csv_path)
    message("[SAVE] ", normalizePath(csv_path, mustWork = FALSE))
    ok_csv <- TRUE
  }
  saveRDS(reg, rds_path)
  message("[SAVE] ", normalizePath(rds_path, mustWork = FALSE))
  invisible(ok_csv)
}

# ───────────────────────────── Ensure packages ─────────────────────────────
options(repos = cran_mirror)
for (p in c("data.table","stringr","dplyr","purrr","readr")) ensure_pkg(p)
ensure_hlapro()

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(dplyr)
  library(purrr)
  library(readr)
  library(hlapro)
})

# ───────────────────────────── I/O ─────────────────────────────
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[INFO] MFI CSV: ", mfi_path)
message("[INFO] Output dir: ", out_dir)
message("[INFO] Verified-only eplets: ", verified_only)

# fread guesses sep and handles decimal commas; keep text as UTF-8
df <- data.table::fread(mfi_path, data.table = FALSE, encoding = "UTF-8", fill = TRUE)

# ───────────────────────────── Detect antigens ─────────────────────────────
antigens <- detect_antigens(names(df))
message("[INFO] Detected ", length(antigens), " antigen columns.")

# ───────────────────────────── Load Registry ─────────────────────────────
message("[INFO] Loading HLA Eplet Registry (hlapro)…")
reg <- hlapro::load_eplet_registry()

# ───────────────────────────── Build mapping ─────────────────────────────
lookup_safe <- function(reg, alleles, verified_only = FALSE) {
  # Try with verified_only argument if supported; otherwise fall back.
  out <- tryCatch(hlapro::lookup_eplets(reg, alleles, verified_only = verified_only),
                  error = function(e) NULL)
  if (is.null(out)) {
    out <- hlapro::lookup_eplets(reg, alleles)
  }
  out
}

eps_list <- lookup_safe(reg, antigens, verified_only = verified_only)

if (is.null(eps_list)) {
  stop("Failed to get eplets from hlapro::lookup_eplets(). Please update 'hlapro' and try again.")
}

mapping <- tibble::tibble(
  Antigen = names(eps_list),
  Eplets  = purrr::map_chr(eps_list, to_semicolon)
) %>%
  arrange(Antigen)

# ───────────────────────────── Save outputs ─────────────────────────────
map_path <- file.path(out_dir, "eplet_map.csv")
readr::write_csv(mapping, map_path)
message("[SAVE] ", normalizePath(map_path, mustWork = FALSE))

# Save a Registry snapshot for reproducibility (CSV if possible + RDS always)
save_registry_snapshot(reg, out_dir)

# Minimal README
readme_path <- file.path(out_dir, "README.txt")
readme_txt <- c(
  "eplet_ref – local reference for HLA eplet mapping",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("Source MFI file: ", normalizePath(mfi_path, mustWork = FALSE)),
  paste0("Verified-only: ", verified_only),
  "",
  "Files:",
  "  • eplet_map.csv                        – Antigen → semicolon-separated eplets",
  "  • eplet_registry_snapshot_YYYYMMDD.rds – Full Registry snapshot (readRDS)",
  "  • eplet_registry_snapshot_YYYYMMDD.csv – Flat snapshot (if available)",
  "  • MANIFEST.txt                         – Build metadata",
  "",
  "Notes:",
  "- If any Antigen has an empty Eplets field, it means no eplets were found in the current Registry.",
  "- Re-run this script whenever the Registry updates, to refresh the local reference."
)
writeLines(readme_txt, con = readme_path)
message("[SAVE] ", normalizePath(readme_path, mustWork = FALSE))

# Manifest with versions
write_manifest(out_dir, mapping,
               params = list(mfi_path = mfi_path, verified_only = verified_only),
               pkgs = c("data.table","stringr","dplyr","purrr","readr","hlapro"))

message("[DONE] eplet_map + Registry reference saved to: ", out_dir)
