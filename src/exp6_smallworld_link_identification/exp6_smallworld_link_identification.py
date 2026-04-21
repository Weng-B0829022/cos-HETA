# -*- coding: utf-8 -*-
"""
實驗 3 — Fig. 6：WS / NW 在 p = 0.004 / 0.064 / 0.256 的連結識別
對應論文 §3.2, Fig. 6a–f（圓形佈局）
六張子圖 = WS（上排 a-c） + NW（下排 d-f）
連結顏色：藍 = bond、紅 = local bridge、綠 = global bridge、黃 = silk
輸出：fig6.png
"""
import os
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from heta import HETA as heta

# 本腳本位於 HETA/src/exp6_smallworld_link_identification/
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

P_VALUES = ["0.004", "0.064", "0.256"]
TIMES = int(os.environ.get("HETA_TIMES", 1000))

COLOR_MAP = {
    "Bond1":         "blue", "Bond2":         "blue", "Bond3":         "blue",
    "Local_Bridge1": "red",  "Local_Bridge2": "red",  "Local_Bridge3": "red",
    "Global_Bridge": "green",
    "Silk":          "gold",
}


def draw_circular_with_types(G, ax, title):
    pos = nx.circular_layout(G)
    edge_colors = [COLOR_MAP.get(G[u][v].get("type"), "gray")
                   for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors,
                           width=0.6, alpha=0.7)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=20,
                           node_color="lightgray")
    ax.set_title(title, fontsize=11)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    print(f"HETA_TIMES = {TIMES}（每網路產生 {TIMES} 個隨機化網路）")
    fig, axes = plt.subplots(2, 3, figsize=(16, 11))

    for col, p_str in enumerate(P_VALUES):
        # 上排：WS
        ws_path = os.path.join(DATA_DIR, f"ws_swn_{p_str}.net")
        G = nx.Graph(nx.read_pajek(ws_path))
        G.remove_edges_from(nx.selfloop_edges(G))
        print(f"  [{col+1}/3] WS p={p_str} ...")
        G, _, _ = heta.bridge_or_bond(G, times=TIMES)
        draw_circular_with_types(G, axes[0][col],
                                 f"({chr(97+col)}) WS  p={p_str}")

        # 下排：NW
        nw_path = os.path.join(DATA_DIR, f"nws_swn_{p_str}.net")
        G = nx.Graph(nx.read_pajek(nw_path))
        G.remove_edges_from(nx.selfloop_edges(G))
        print(f"  [{col+1}/3] NW p={p_str} ...")
        G, _, _ = heta.bridge_or_bond(G, times=TIMES)
        draw_circular_with_types(G, axes[1][col],
                                 f"({chr(100+col)}) NW  p={p_str}")

    plt.suptitle("Fig. 6: WS (top) vs NW (bottom) small-world link identification",
                 fontsize=13)
    plt.tight_layout()
    out = "fig6.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"產生：{os.path.abspath(out)}")


if __name__ == "__main__":
    main()
