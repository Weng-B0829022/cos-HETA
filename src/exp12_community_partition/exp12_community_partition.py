# -*- coding: utf-8 -*-
"""
實驗 9 — Fig. 12：4 個代表網路（rdgam, women, camp92, karate）的階層社群分割
對應論文 §3.4, Fig. 12a–d
節點以社群顏色區分，連結以 HETA 類型顏色區分。
輸出：fig12.png
"""
import os
import networkx as nx
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from heta import HETA as heta

# 本腳本位於 HETA/src/exp12_community_partition/
# 輸入資料：HETA 根目錄下的 data/
# 輸出：HETA/experiment/experiment_<timestamp>/
#   - 由 auto.py 以 HETA_OUTPUT_DIR 環境變數指定
#   - 獨立執行時自動建立 experiment_<timestamp>/ 備援目錄
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_HETA_ROOT  = os.path.abspath(os.path.join(_SCRIPT_DIR, "..", ".."))
DATA_DIR = os.path.join(_HETA_ROOT, "data", "net")
POS_DIR  = os.path.join(_HETA_ROOT, "data", "pos")

_out = os.environ.get("HETA_OUTPUT_DIR")
if not _out:
    from datetime import datetime as _dt
    _out = os.path.join(_HETA_ROOT, "experiment",
                        f"experiment_{_dt.now().strftime('%Y%m%d_%H%M%S')}")
os.makedirs(_out, exist_ok=True)
os.chdir(_out)

TIMES = int(os.environ.get("HETA_TIMES", 1000))

NETWORKS = [
    ("(a) rdgam",  "14p.net",    "14p.csv"),
    ("(b) women",  "women.net",  "women.csv"),
    ("(c) camp92", "USA.net",    "USA.csv"),
    ("(d) karate", "karate.net", "karate.csv"),
]

EDGE_COLOR = {
    "Bond1": "blue", "Bond2": "blue", "Bond3": "blue",
    "Local_Bridge1": "red", "Local_Bridge2": "red", "Local_Bridge3": "red",
    "Global_Bridge": "green",
    "Silk": "gold",
}
COMMUNITY_COLORS = ["blue", "green", "red", "purple", "orange",
                    "cyan", "magenta", "brown", "olive", "teal"]


def load_positions(csv_path, G):
    if not os.path.exists(csv_path):
        return nx.spring_layout(G, seed=42)
    df = pd.read_csv(csv_path)
    pos = {str(row["node"]): (row["xx"], row["yy"]) for _, row in df.iterrows()}
    missing = set(G.nodes()) - set(pos.keys())
    if missing:
        sub = G.subgraph(missing)
        pos.update(nx.spring_layout(sub, seed=42))
    return pos


def flatten_communities(communities):
    """heta.community_sorting() 回傳結構不一定統一，將其展平為 [set(), set(), ...]

    community_sorting 的正式回傳是 4-tuple：
        (nodes_sorted, h_tree_dic, flat_cdic, kmax)
    其中 flat_cdic = {leader: [member, ...]}。舊版通用展平邏輯會把
    flat_cdic 的 list 值拆成 singleton，導致社群偵測失效（同 Table 2 bug）。
    這裡優先走 4-tuple 專用分支。
    """
    if (isinstance(communities, tuple) and len(communities) == 4
            and isinstance(communities[2], dict)):
        flat_cdic = communities[2]
        return [set(str(n) for n in nodes) for nodes in flat_cdic.values()]

    out = []
    if communities is None:
        return out
    if isinstance(communities, (set, frozenset)):
        return [set(str(n) for n in communities)]
    if isinstance(communities, dict):
        for v in communities.values():
            out.extend(flatten_communities(v))
        return out
    if isinstance(communities, (list, tuple)):
        for item in communities:
            if isinstance(item, (set, frozenset)):
                out.append(set(str(n) for n in item))
            elif isinstance(item, (list, tuple)):
                out.extend(flatten_communities(item))
            elif isinstance(item, dict):
                out.extend(flatten_communities(item))
            else:
                # 單一節點
                out.append({str(item)})
    return out


def main():
    print(f"HETA_TIMES = {TIMES}")
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))

    for idx, (title, net_file, pos_file) in enumerate(NETWORKS):
        ax = axes[idx // 2][idx % 2]
        print(f"  {title} ...")

        G = nx.Graph(nx.read_pajek(os.path.join(DATA_DIR, net_file)))
        G.remove_edges_from(nx.selfloop_edges(G))
        G, _, _ = heta.bridge_or_bond(G, times=TIMES)

        try:
            communities = heta.community_sorting(G)
        except Exception as e:
            print(f"    community_sorting 失敗：{e}")
            communities = []
        comm_list = [c for c in flatten_communities(communities) if len(c) > 1]

        node_community = {str(n): -1 for n in G.nodes()}
        for cid, comm in enumerate(comm_list):
            for n in comm:
                if str(n) in node_community:
                    node_community[str(n)] = cid

        pos = load_positions(os.path.join(POS_DIR, pos_file), G)
        node_colors = []
        for n in G.nodes():
            cid = node_community.get(str(n), -1)
            if cid == -1:
                node_colors.append("lightgray")
            else:
                node_colors.append(COMMUNITY_COLORS[cid % len(COMMUNITY_COLORS)])

        edge_colors = [EDGE_COLOR.get(G[u][v].get("type"), "gray")
                       for u, v in G.edges()]
        nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors,
                               width=0.6, alpha=0.7)
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=60,
                               node_color=node_colors,
                               edgecolors="black", linewidths=0.5)
        ax.set_title(f"{title}  (groups={len(comm_list)})", fontsize=12)
        ax.set_aspect("equal")
        ax.axis("off")

    plt.suptitle("Fig. 12: Hierarchical community partition (rdgam, women, camp92, karate)",
                 fontsize=13)
    plt.tight_layout()
    out = "fig12.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"產生：{os.path.abspath(out)}")


if __name__ == "__main__":
    main()
