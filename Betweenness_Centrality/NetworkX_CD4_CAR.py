#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NetworkX – CAR_CD4 across SCENIC_CAR runs

• Per-run densification (prefer adjacencies CSV; keep TOP_EDGES_PER_TF per TF; hit TARGET_EDGE_RANGE)
• Overlap graphs (relaxed nodes/edges) with presence-weighted edges
• Strict overlap with auto-tuned edge threshold
• Core extraction (k-core)
• Label top nodes by (weighted) betweenness in plots
• Writes node/edge presence counts and multi‑strictness node‑overlap summaries
"""
from __future__ import annotations
import argparse
import glob
import math
import pickle
import random
import re
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

SPRING_SEED = 4_572_321
RNG_SEED = 42
random.seed(RNG_SEED); np.random.seed(RNG_SEED)

def scenic_idx(path_str): 
    m = re.search(r"SCENIC_CAR_(\d+)", Path(path_str).as_posix()); return int(m.group(1)) if m else 10**9

def _extract_edges_from_pickle(path):
    import ctxcore  # noqa: F401
    edges = []
    with open(path, "rb") as fh: regs = pickle.load(fh)
    for reg in regs:
        tf = getattr(reg, "name", None); g2w = getattr(reg, "gene2weight", {}) or {}
        if tf:
            for tgt, w in g2w.items():
                if pd.notna(w): edges.append((str(tf), str(tgt), float(w)))
    return edges

def _extract_edges_from_csv(path):
    df = pd.read_csv(path)
    cmap = {c.lower(): c for c in df.columns}
    TF  = cmap.get("tf", "TF"); TGT = cmap.get("target", "target"); WEI = cmap.get("importance", "importance")
    out = []
    for _, r in df.iterrows():
        if pd.notna(r[WEI]): out.append((str(r[TF]), str(r[TGT]), float(r[WEI])))
    return out

def extract_edges(path_str, cell):
    p = Path(path_str)
    if p.suffix.lower() == ".csv": return _extract_edges_from_csv(p)
    try: return _extract_edges_from_pickle(p)
    except Exception:
        csv_alt = p.with_name(f"adjacencies_{cell}.csv")
        if csv_alt.is_file():
            warnings.warn(f"[INFO] Using CSV fallback: {csv_alt.name}")
            return _extract_edges_from_csv(csv_alt)
        raise

def weight_transform(values, how=None):
    v = np.asarray(values, dtype=float)
    if how is None: return v
    if how == "log1p": return np.log1p(np.maximum(v, 0.0))
    if how == "rank_global":
        order = v.argsort(); ranks = np.empty_like(order, dtype=float); ranks[order] = np.linspace(0.0, 1.0, len(v))
        return ranks
    raise ValueError("Unknown WEIGHT_TRANSFORM")

def build_graph_from_edges(edges, wmin, dmin):
    G = nx.DiGraph()
    for u, v, w in edges:
        if w >= wmin: G.add_edge(u, v, weight=float(w))
    if dmin > 0: G.remove_nodes_from([n for n, d in G.degree() if d < dmin])
    return G

def largest_wcc(G):
    if G.number_of_nodes() == 0: return nx.DiGraph()
    try: comp = max(nx.weakly_connected_components(G), key=len); return G.subgraph(comp).copy()
    except ValueError: return nx.DiGraph()

def bc_and_out(G):
    if G.number_of_nodes() == 0: return {}, {}
    Gu = G.to_undirected()
    for _, _, d in Gu.edges(data=True): d["dist"] = 1.0 / max(d.get("weight", 1e-12), 1e-12)
    bc = nx.betweenness_centrality(Gu, weight="dist"); wout = dict(G.out_degree(weight="weight"))
    return bc, wout

def rank(d): return sorted(d.items(), key=lambda x: x[1], reverse=True)

def draw_graph(G, tag, hubs, bc, out_png, node_size_range=(300, 6000), label_top=50):
    if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
        print(f"[WARN] {tag}: graph empty; skipping plot"); return
    coms = list(nx.community.label_propagation_communities(G.to_undirected()))
    com_idx = {n: i for i, com in enumerate(coms) for n in com}; node_col = [com_idx.get(n, 0) for n in G]
    n_com = (max(node_col) + 1) if node_col else 1
    vals = np.array([bc.get(n, 0.0) for n in G])
    sizes = (np.full(G.number_of_nodes(), node_size_range[0])
             if vals.size == 0 or vals.ptp() == 0
             else np.interp(vals, (vals.min(), vals.max()), node_size_range))
    sz = {n: s for n, s in zip(G.nodes(), sizes)}
    pos = nx.spring_layout(G, k=0.15, seed=SPRING_SEED, weight="weight")
    cmap = plt.colormaps.get_cmap("viridis")
    fig, ax = plt.subplots(figsize=(20, 15))
    nx.draw_networkx_edges(G, pos, edge_color="gainsboro", alpha=0.25, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=node_col, node_size=[sz[n] for n in G], cmap=cmap,
                           linewidths=0, ax=ax)
    texts = nx.draw_networkx_labels(G, pos, labels={n: n for n in hubs[:label_top]},
                                    font_size=14, font_color="white", ax=ax)
    for t in texts.values(): t.set_path_effects([pe.Stroke(linewidth=2.2, foreground="black"), pe.Normal()])
    sm = plt.cm.ScalarMappable(norm=mcolors.Normalize(0, n_com - 1), cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.01); cbar.set_label("Community index", fontsize=12)
    for label, s in zip(["low BC","mid","high"], [node_size_range[0], int(sum(node_size_range)/2), node_size_range[1]]):
        ax.scatter([], [], s=s, c="grey", alpha=0.6, edgecolors="none", label=label)
    ax.legend(scatterpoints=1, frameon=False, fontsize=10, loc="upper left", bbox_to_anchor=(0.02, 0.95),
              title="Weighted\nbetweenness")
    ax.set_title(tag, fontdict={"fontsize": 22, "fontweight": "bold"})
    ax.set_axis_off(); ax.margins(0.1, 0.05); fig.tight_layout()
    out_png = Path(out_png); fig.savefig(out_png, dpi=300); fig.savefig(out_png.with_suffix(".svg")); plt.close()
    print(f"[INFO] Saved {out_png.name} / {out_png.stem}.svg")

def build_per_run_graph(raw_edges, dmin, target_edges=(1500, 20000), top_per_tf=50, transform=None):
    if len(raw_edges) == 0: return nx.DiGraph(), nx.DiGraph(), None, None, 0, 0
    edges = raw_edges
    if top_per_tf is not None:
        by_tf = defaultdict(list)
        for u, v, w in raw_edges: by_tf[u].append((u, v, float(w)))
        edges = []; 
        for tf, lst in by_tf.items():
            lst.sort(key=lambda x: x[2], reverse=True); edges.extend(lst[:top_per_tf])
    if transform is not None:
        ws = [w for _,_,w in edges]; ws_t = weight_transform(ws, transform)
        edges = [(u, v, float(wt)) for (u, v, _), wt in zip(edges, ws_t)]
    weights = np.array([w for _,_,w in edges], dtype=float); uniq = np.unique(weights); uniq.sort()
    chosen_w = None; chosen_G = None
    for wmin in uniq[::-1]:
        G = build_graph_from_edges(edges, wmin, dmin=0); e = G.number_of_edges()
        if target_edges[0] <= e <= target_edges[1]: chosen_w, chosen_G = float(wmin), G; break
    if chosen_G is None: chosen_w = float(uniq[0]); chosen_G = build_graph_from_edges(edges, chosen_w, dmin=0)
    if dmin > 0: chosen_G.remove_nodes_from([n for n, d in chosen_G.degree() if d < dmin])
    G_for_agg = chosen_G; G_for_plot = largest_wcc(chosen_G)
    return G_for_agg, G_for_plot, chosen_w, dmin, len(raw_edges), G_for_agg.number_of_edges()

def aggregate(graphs_dict, keep_fn):
    agg = nx.DiGraph()
    for G in graphs_dict.values():
        for u, v, d in G.edges(data=True):
            if not keep_fn(u, v): continue
            if agg.has_edge(u, v): agg[u][v]["_w"].append(d["weight"])
            else: agg.add_edge(u, v, _w=[d["weight"]])
    for _, _, d in agg.edges(data=True): d["weight"] = float(np.mean(d["_w"])); del d["_w"]
    return agg

def main():
    ap = argparse.ArgumentParser(description="NetworkX consensus/overlap across SCENIC_CAR runs (CAR_CD4).")
    ap.add_argument("--base", required=True, type=Path,
                    help="Base folder containing SCENIC_CAR_* subfolders (as in the Zenodo optional path).")
    ap.add_argument("--cell", default="CAR_CD4")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--top-edges-per-tf", type=int, default=50)
    ap.add_argument("--target-edge-range", nargs=2, type=int, default=[1500, 20000])
    ap.add_argument("--weight-transform", choices=["log1p", "rank_global"], default=None)
    ap.add_argument("--low-deg", type=int, default=0)
    ap.add_argument("--overlap-fraction", type=float, default=0.6)
    ap.add_argument("--edge-overlap-fraction", type=float, default=0.6)
    ap.add_argument("--presence-weight-exp", type=float, default=1.0)
    ap.add_argument("--strict-node-overlap-fraction", type=float, default=0.95)
    ap.add_argument("--strict-auto-min-edges", type=int, default=50)
    ap.add_argument("--label-top", type=int, default=50)
    args = ap.parse_args()

    outdir = args.outdir; outdir.mkdir(parents=True, exist_ok=True)

    runs_all = sorted(glob.glob(str(args.base / "SCENIC_CAR_*")), key=scenic_idx)
    if not runs_all: raise SystemExit(f"No SCENIC_CAR_* folders found under {args.base}")
    PKL_OR_CSV = {}
    for r in runs_all:
        pkl = Path(r) / args.cell / f"regulons_{args.cell}.pkl"
        csv = Path(r) / args.cell / f"adjacencies_{args.cell}.csv"
        if csv.is_file(): PKL_OR_CSV[Path(r).name] = str(csv)
        elif pkl.is_file(): PKL_OR_CSV[Path(r).name] = str(pkl)
    for k in sorted(PKL_OR_CSV, key=scenic_idx): print(f"  {k}: {PKL_OR_CSV[k]}")
    if not PKL_OR_CSV: raise SystemExit(f"No pkl/CSV found for {args.cell} in detected runs.")

    graphs = {}; node_presence = Counter()
    for tag in sorted(PKL_OR_CSV, key=scenic_idx):
        infile = PKL_OR_CSV[tag]; raw_edges = extract_edges(infile, args.cell)
        G_agg, G_plot, wmin_used, dmin_used, n_raw, n_kept = build_per_run_graph(
            raw_edges, dmin=args.low_deg, target_edges=tuple(args.target_edge_range),
            top_per_tf=args.top_edges_per_tf, transform=args.weight_transform)
        if G_agg.number_of_nodes() == 0 or G_agg.number_of_edges() == 0:
            print(f"[WARN] {tag} empty after densification."); continue
        node_presence.update([n for n, deg in G_agg.out_degree() if deg > 0])
        bc_plot, _ = bc_and_out(G_plot)
        hubs_plot = [n for n, _ in rank(bc_plot)[:args.label_top]]
        draw_graph(G_plot, tag, hubs_plot, bc_plot, outdir / f"{tag}_network.png", label_top=args.label_top)
        graphs[tag] = G_agg

    pd.DataFrame(sorted(node_presence.items()), columns=["node","runs_with_outgoing_edge"])\
      .assign(fraction=lambda d: d["runs_with_outgoing_edge"]/max(len(graphs),1))\
      .sort_values(["runs_with_outgoing_edge","node"], ascending=[False, True])\
      .to_csv(outdir / "node_presence_counts.tsv", sep="\t", index=False)

    if not graphs: raise SystemExit("All graphs empty after densification.")

    # consensus
    consensus = aggregate(graphs, lambda u, v: True)
    nx.write_graphml(consensus, outdir / f"consensus_{args.cell}.graphml")
    bc_c, _ = bc_and_out(consensus); hubs_c = [n for n, _ in rank(bc_c)[:args.label_top]]
    draw_graph(consensus, "CONSENSUS", hubs_c, bc_c, outdir / "CONSENSUS_network.png", label_top=args.label_top)

    # relaxed overlaps
    n_runs = len(graphs)
    min_runs_nodes = int(math.ceil(args.overlap_fraction * n_runs))
    occ_nodes = Counter(); [occ_nodes.update(list(G.nodes())) for G in graphs.values()]
    overlap_nodes = {n for n, c in occ_nodes.items() if c >= min_runs_nodes}
    def keep_if_nodes(u, v): return (u in overlap_nodes) and (v in overlap_nodes)
    overlapN = aggregate(graphs, keep_if_nodes)
    for u, v, d in overlapN.edges(data=True):
        pres = sum(G.has_edge(u, v) for G in graphs.values()) / n_runs
        d["presence_fraction"] = pres; d["weight"] *= pres ** args.presence_weight_exp
    nx.write_graphml(overlapN, outdir / f"overlapNodes_{args.cell}.graphml")
    if overlapN.number_of_edges() > 0:
        bc_n, _ = bc_and_out(overlapN); hubs_n = [n for n, _ in rank(bc_n)[:args.label_top]]
        draw_graph(overlapN, f"OVERLAP-NODES (≥{min_runs_nodes}/{n_runs})", hubs_n, bc_n,
                   outdir / "overlapNodes_network.png", label_top=args.label_top)

    min_runs_edges = int(math.ceil(args.edge_overlap_fraction * n_runs))
    def keep_if_edge(u, v): return sum(G.has_edge(u, v) for G in graphs.values()) >= min_runs_edges
    overlapE = aggregate(graphs, keep_if_edge)
    for u, v, d in overlapE.edges(data=True):
        pres = sum(G.has_edge(u, v) for G in graphs.values()) / n_runs
        d["presence_fraction"] = pres; d["weight"] *= pres ** args.presence_weight_exp
    nx.write_graphml(overlapE, outdir / f"overlapEdges_{args.cell}.graphml")
    if overlapE.number_of_edges() > 0:
        bc_e, _ = bc_and_out(overlapE); hubs_e = [n for n, _ in rank(bc_e)[:args.label_top]]
        draw_graph(overlapE, f"OVERLAP-EDGES (≥{min_runs_edges}/{n_runs})", hubs_e, bc_e,
                   outdir / "overlapEdges_network.png", label_top=args.label_top)

    # strict overlap (auto-tuned edge fraction)
    min_runs_nodes_strict = int(math.ceil(args.strict_node_overlap_fraction * n_runs))
    occ_nodes_strict = Counter()
    for G in graphs.values(): occ_nodes_strict.update([n for n, deg in G.out_degree() if deg > 0])
    strict_nodes = {n for n, c in occ_nodes_strict.items() if c >= min_runs_nodes_strict}
    edge_count = Counter(); [edge_count.update(list(G.edges())) for G in graphs.values()]
    pd.DataFrame([(u, v, c, c/n_runs, (u in strict_nodes) and (v in strict_nodes))
                  for (u, v), c in edge_count.items()],
                 columns=["u","v","presence_count","presence_fraction","both_nodes_strict"])\
      .sort_values(["presence_count","u","v"], ascending=[False, True, True])\
      .to_csv(outdir / "edge_presence_counts.tsv", sep="\t", index=False)
    grid = [1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40]
    chosen_frac = None; chosen_edges = 0
    for frac in grid:
        thr_c = int(math.ceil(frac * n_runs))
        n_edges = sum(1 for (u, v), c in edge_count.items()
                      if c >= thr_c and (u in strict_nodes) and (v in strict_nodes))
        if n_edges >= args.strict_auto_min_edges and n_edges > 0:
            chosen_frac, chosen_edges = frac, n_edges; break
    if chosen_frac is None:
        best = []
        for frac in grid:
            thr_c = int(math.ceil(frac * n_runs))
            n_edges = sum(1 for (u, v), c in edge_count.items()
                          if c >= thr_c and (u in strict_nodes) and (v in strict_nodes))
            best.append((n_edges, frac))
        chosen_edges, chosen_frac = max(best); 
        if chosen_edges == 0: chosen_frac = min(grid)
    thr_edges = int(math.ceil(chosen_frac * n_runs))
    def keep_if_strict(u, v):
        c = edge_count[(u, v)]
        return (u in strict_nodes) and (v in strict_nodes) and (c >= thr_edges)
    overlapStrict = aggregate(graphs, keep_if_strict)
    for u, v, d in overlapStrict.edges(data=True):
        pres = edge_count[(u, v)] / n_runs
        d["presence_fraction"] = pres; d["weight"] *= pres ** args.presence_weight_exp
    nx.write_graphml(overlapStrict, outdir / f"overlapStrict_{args.cell}.graphml")
    if overlapStrict.number_of_edges() > 0:
        bc_s, _ = bc_and_out(overlapStrict); hubs_s = [n for n, _ in rank(bc_s)[:args.label_top]]
        draw_graph(overlapStrict, f"OVERLAP-STRICT (nodes ≥{min_runs_nodes_strict}/{n_runs}, edges ≥{thr_edges}/{n_runs})",
                   hubs_s, bc_s, outdir / "overlapStrict_network.png", label_top=args.label_top)

    # core (k-core auto up to k=5)
    def kcore_auto(G, k_target=5):
        if G.number_of_edges() == 0: return nx.DiGraph(), 0
        Gu = G.to_undirected(); core_num = nx.core_number(Gu); kmax = max(core_num.values()) if core_num else 0
        for k in range(min(k_target, kmax), 0, -1):
            nodes_k = nx.k_core(Gu, k=k).nodes
            H_nodes = set(nodes_k)
            if len(H_nodes) >= 2: return G.subgraph(H_nodes).copy(), k
        return nx.DiGraph(), 0
    core_overlap, k_used = kcore_auto(overlapStrict, 5 if overlapStrict.number_of_edges() > 0 else 5)
    if core_overlap.number_of_edges() == 0: core_overlap, k_used = kcore_auto(consensus, 5)
    if core_overlap.number_of_edges() > 0:
        bc_core, _ = bc_and_out(core_overlap); hubs_core = [n for n, _ in rank(bc_core)[:args.label_top]]
        draw_graph(core_overlap, f"OVERLAP-CORE (k={k_used})", hubs_core, bc_core, outdir / "overlapCore_network.png",
                   label_top=args.label_top)
        nx.write_graphml(core_overlap, outdir / f"overlapCore_{args.cell}.graphml")
    print(f"\n[ALL DONE] – Outputs are in {outdir}")

if __name__ == "__main__":
    main()
