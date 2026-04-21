# -*- coding: utf-8 -*-
"""
實驗 2 — Fig. 5：16 個實證網路的 R¹ 箱型圖
對應論文 §3.1, Fig. 5
輸出：fig5.png
驗收：箱子由中位數高至低排序，可由兩條紅虛線分為 7 / 4 / 5 三群
"""
import os
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 本腳本位於 HETA/src/exp5_real16_r1_boxplot/
# 輸入資料：HETA 根目錄下的 data/
# 輸出：HETA/experiment/experiment_<timestamp>/
#   - 由 auto.py 以 HETA_OUTPUT_DIR 環境變數指定
#   - 獨立執行時自動建立 experiment_<timestamp>/ 備援目錄
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_HETA_ROOT  = os.path.abspath(os.path.join(_SCRIPT_DIR, "..", ".."))
DATA_DIR = os.path.join(_HETA_ROOT, "data", "net")

_out = os.environ.get("HETA_OUTPUT_DIR")
if not _out:
    from datetime import datetime as _dt
    _out = os.path.join(_HETA_ROOT, "experiment",
                        f"experiment_{_dt.now().strftime('%Y%m%d_%H%M%S')}")
os.makedirs(_out, exist_ok=True)
os.chdir(_out)

# 順序 = Fig. 5（依中位數由高至低）
NETWORKS = [
    ("rdgam",       "14p.net"),
    ("women",       "women.net"),
    ("lesmis",      "lesmis.net"),
    ("jazz",        "jazz.net"),
    ("karate",      "karate.net"),
    ("camp92",      "USA.net"),
    ("ragusa16",    "Ragusa16.net"),
    ("football",    "football.net"),
    ("dolphins",    "dolphins.net"),
    ("celegans",    "celegans.net"),
    ("leader",      "leader.net"),
    ("prisonInter", "prisonInter.net"),
    ("k-core",      "k-core.net"),
    ("florentine",  "families.net"),
    ("ba_sfn",      "ba_sfn.net"),
    ("s208",        "s208.net"),
]


def first_layer_cn_ratio(G, u, v):
    nbr_u = set(G.neighbors(u)) - {v}
    nbr_v = set(G.neighbors(v)) - {u}
    if not nbr_u or not nbr_v:
        return 0.0
    return len(nbr_u & nbr_v) / min(len(nbr_u), len(nbr_v))


def main():
    data, labels, medians = [], [], []
    for name, filename in NETWORKS:
        G = nx.Graph(nx.read_pajek(os.path.join(DATA_DIR, filename)))
        G.remove_edges_from(nx.selfloop_edges(G))
        r1 = [first_layer_cn_ratio(G, u, v) for u, v in G.edges()]
        data.append(r1)
        labels.append(name)
        medians.append(np.median(r1) if r1 else 0.0)

    fig, ax = plt.subplots(figsize=(13, 5.5))
    try:
        ax.boxplot(data, tick_labels=labels)
    except TypeError:
        ax.boxplot(data, labels=labels)
    ax.set_ylabel(r"$R^1_{i,j}$")
    ax.set_xlabel("Network")
    ax.set_ylim(-0.05, 1.05)
    ax.tick_params(axis="x", rotation=45)

    # 兩條紅虛線：7 / 4 / 5 分群（與論文一致）
    ax.axvline(x=7.5,  color="red", linestyle="--", alpha=0.6, linewidth=1.2)
    ax.axvline(x=11.5, color="red", linestyle="--", alpha=0.6, linewidth=1.2)

    plt.tight_layout()
    out = "fig5.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"產生：{os.path.abspath(out)}")

    print("\n[驗收] 16 個網路的 R¹ 中位數：")
    for name, med in zip(labels, medians):
        print(f"  {name:<13} median = {med:.3f}")


if __name__ == "__main__":
    main()
