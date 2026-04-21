# -*- coding: utf-8 -*-
"""
實驗 4 — Fig. 7：16 個實證網路的連結識別
對應論文 §3.2, Fig. 7
4×4 子圖網格；每格一個網路。
盡量用 data/pos/*.csv 預存座標做像素級還原；缺漏節點以 spring_layout 補齊。
輸出：fig7.png
"""
import os
import networkx as nx
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from heta import HETA as heta

os.chdir(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = "../data/net"
POS_DIR  = "../data/pos"

TIMES = int(os.environ.get("HETA_TIMES", 1000))

# 順序 = Fig. 7（依 bond 比例由高至低）
NETWORKS = [
    ("rdgam",       "14p.net",         "14p.csv"),
    ("women",       "women.net",       "women.csv"),
    ("football",    "football.net",    "football.csv"),
    ("jazz",        "jazz.net",        "jazz.csv"),
    ("lesmis",      "lesmis.net",      "lesmis.csv"),
    ("camp92",      "USA.net",         "USA.csv"),
    ("prisonInter", "prisonInter.net", "prisonInter.csv"),
    ("karate",      "karate.net",      "karate.csv"),
    ("k-core",      "k-core.net",      "k-core.csv"),
    ("ragusa16",    "Ragusa16.net",    "Ragusa16.csv"),
    ("dolphins",    "dolphins.net",    "dolphins.csv"),
    ("celegans",    "celegans.net",    "celegans.csv"),
    ("s208",        "s208.net",        "s208.csv"),
    ("ba_sfn",      "ba_sfn.net",      "ba_sfn.csv"),
    ("leader",      "leader.net",      "leader.csv"),
    ("florentine",  "families.net",    "families.csv"),
]

COLOR_MAP = {
    "Bond1": "blue", "Bond2": "blue", "Bond3": "blue",
    "Local_Bridge1": "red", "Local_Bridge2": "red", "Local_Bridge3": "red",
    "Global_Bridge": "green",
    "Silk": "gold",
}


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


def main():
    print(f"HETA_TIMES = {TIMES}")
    fig, axes = plt.subplots(4, 4, figsize=(18, 18))

    for i, (label, net_file, pos_file) in enumerate(NETWORKS):
        row, col = i // 4, i % 4
        ax = axes[row][col]
        print(f"  [{i+1:>2}/16] {label} ...")
        G = nx.Graph(nx.read_pajek(os.path.join(DATA_DIR, net_file)))
        G.remove_edges_from(nx.selfloop_edges(G))
        G, _, _ = heta.bridge_or_bond(G, times=TIMES)

        pos = load_positions(os.path.join(POS_DIR, pos_file), G)
        edge_colors = [COLOR_MAP.get(G[u][v].get("type"), "gray")
                       for u, v in G.edges()]
        nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors,
                               width=0.5, alpha=0.7)
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=15,
                               node_color="lightgray")
        ax.set_title(f"({chr(97+i)}) {label}", fontsize=10)
        ax.set_aspect("equal")
        ax.axis("off")

    plt.suptitle("Fig. 7: Link type identification of 16 networks", fontsize=14)
    plt.tight_layout()
    out = "fig7.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"產生：{os.path.abspath(out)}")


if __name__ == "__main__":
    main()
