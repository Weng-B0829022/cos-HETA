# -*- coding: utf-8 -*-
"""
實驗 1 — Fig. 4：WS/NW 32 個小世界網路 R¹ 箱型圖
對應論文 §3.1, Fig. 4a (WS) 與 Fig. 4b (NW)
輸入：data/net/ws_swn_*.net、data/net/nws_swn_*.net（各 16 個 p 值）
輸出：fig4.png
"""
import os
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 本腳本位於 HETA/src/exp4_smallworld_r1_boxplot/
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

P_VALUES = [
    (0.0,   "0.0"),   (0.001, "0.001"), (0.002, "0.002"), (0.004, "0.004"),
    (0.008, "0.008"), (0.016, "0.016"), (0.032, "0.032"), (0.064, "0.064"),
    (0.128, "0.128"), (0.256, "0.256"), (0.384, "0.384"), (0.512, "0.512"),
    (0.640, "0.640"), (0.768, "0.768"), (0.896, "0.896"), (1.0,   "1.0"),
]


def first_layer_cn_ratio(G, u, v):
    """論文 Eq. 1：第一層共同鄰居比率"""
    nbr_u = set(G.neighbors(u)) - {v}
    nbr_v = set(G.neighbors(v)) - {u}
    if not nbr_u or not nbr_v:
        return 0.0
    return len(nbr_u & nbr_v) / min(len(nbr_u), len(nbr_v))


def compute_r1_distribution(G):
    return [first_layer_cn_ratio(G, u, v) for u, v in G.edges()]


def plot_boxplot(ax, mechanism_label, filename_prefix):
    data, labels = [], []
    for _, p_str in P_VALUES:
        path = os.path.join(DATA_DIR, f"{filename_prefix}_swn_{p_str}.net")
        G = nx.Graph(nx.read_pajek(path))
        G.remove_edges_from(nx.selfloop_edges(G))
        data.append(compute_r1_distribution(G))
        labels.append(p_str)
    try:
        ax.boxplot(data, tick_labels=labels)
    except TypeError:
        ax.boxplot(data, labels=labels)  # 老版 matplotlib 兼容
    title_full = "WS small-world network" if mechanism_label == "a" else "NW small-world network"
    xlabel = "Random rewiring probability" if mechanism_label == "a" \
             else "Shortcut addition probability"
    ax.set_title(f"({mechanism_label}) {title_full}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$R^1_{i,j}$")
    ax.set_ylim(-0.05, 1.05)
    ax.tick_params(axis="x", rotation=45)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    plot_boxplot(axes[0], "a", "ws")
    plot_boxplot(axes[1], "b", "nws")
    plt.tight_layout()
    out = "fig4.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"產生：{os.path.abspath(out)}")

    # 驗收用：印出幾個關鍵 p 的中位數
    print("\n[驗收] WS / NW 各 p 的中位數對照論文 Fig. 4：")
    print(f"  {'p':>6}  {'WS_median':>10}  {'NW_median':>10}")
    for _, p_str in P_VALUES:
        if p_str not in {"0.0", "0.064", "0.128", "0.256", "0.640", "1.0"}:
            continue
        Gw = nx.Graph(nx.read_pajek(os.path.join(DATA_DIR, f"ws_swn_{p_str}.net")))
        Gn = nx.Graph(nx.read_pajek(os.path.join(DATA_DIR, f"nws_swn_{p_str}.net")))
        Gw.remove_edges_from(nx.selfloop_edges(Gw))
        Gn.remove_edges_from(nx.selfloop_edges(Gn))
        mw = np.median(compute_r1_distribution(Gw))
        mn = np.median(compute_r1_distribution(Gn))
        print(f"  {p_str:>6}  {mw:>10.3f}  {mn:>10.3f}")


if __name__ == "__main__":
    main()
