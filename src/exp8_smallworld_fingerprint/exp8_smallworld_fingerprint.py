# -*- coding: utf-8 -*-
"""
實驗 5 — Fig. 8：WS / NW 32 個小世界網路的指紋（fingerprint）堆疊長條圖
對應論文 §3.3, Fig. 8a (WS) 與 Fig. 8b (NW)
指紋 = [bond%, local%, global%, silk%]
快取：../../_cache/fp_ws.pkl, fp_nw.pkl（給 exp10 重複利用）
輸出：fig8.png
"""
import os
import pickle
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from heta import HETA as heta

# 本腳本位於 HETA/src/exp8_smallworld_fingerprint/
# 輸入資料 + 共用快取：HETA 根目錄下的 data/ 與 _cache/
# 輸出：HETA/experiment/experiment_<timestamp>/
#   - 由 auto.py 以 HETA_OUTPUT_DIR 環境變數指定
#   - 獨立執行時自動建立 experiment_<timestamp>/ 備援目錄
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_HETA_ROOT  = os.path.abspath(os.path.join(_SCRIPT_DIR, "..", ".."))
DATA_DIR  = os.path.join(_HETA_ROOT, "data", "net")
CACHE_DIR = os.path.join(_HETA_ROOT, "_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

_out = os.environ.get("HETA_OUTPUT_DIR")
if not _out:
    from datetime import datetime as _dt
    _out = os.path.join(_HETA_ROOT, "experiment",
                        f"experiment_{_dt.now().strftime('%Y%m%d_%H%M%S')}")
os.makedirs(_out, exist_ok=True)
os.chdir(_out)

TIMES = int(os.environ.get("HETA_TIMES", 1000))

P_VALUES = [
    (0.0,   "0.0"),   (0.001, "0.001"), (0.002, "0.002"), (0.004, "0.004"),
    (0.008, "0.008"), (0.016, "0.016"), (0.032, "0.032"), (0.064, "0.064"),
    (0.128, "0.128"), (0.256, "0.256"), (0.384, "0.384"), (0.512, "0.512"),
    (0.640, "0.640"), (0.768, "0.768"), (0.896, "0.896"), (1.0,   "1.0"),
]


def compute_or_load_fingerprints(prefix, label):
    """回傳 16 個 (bond, local, global, silk) 比例向量"""
    cache_path = os.path.join(CACHE_DIR, f"fp_{label}.pkl")
    if os.path.exists(cache_path):
        print(f"  從快取載入：{cache_path}")
        with open(cache_path, "rb") as f:
            return pickle.load(f)

    fps = []
    for _, p_str in P_VALUES:
        print(f"  跑 {label} p={p_str} ...")
        G = nx.Graph(nx.read_pajek(
            os.path.join(DATA_DIR, f"{prefix}_swn_{p_str}.net")))
        G.remove_edges_from(nx.selfloop_edges(G))
        G, _, _ = heta.bridge_or_bond(G, times=TIMES)
        fps.append(heta.fingerprint(G))
    with open(cache_path, "wb") as f:
        pickle.dump(fps, f)
    print(f"  快取已寫入：{cache_path}")
    return fps


def plot_stacked_bar(ax, fingerprints, labels, title, xlabel):
    fp = np.array(fingerprints)
    bond, local, gbr, silk = fp[:, 0], fp[:, 1], fp[:, 2], fp[:, 3]
    x = np.arange(len(labels))
    ax.bar(x, bond,  color="blue",  label="Bond link")
    ax.bar(x, local, bottom=bond,            color="red",   label="Local bridge")
    ax.bar(x, gbr,   bottom=bond+local,      color="green", label="Global bridge")
    ax.bar(x, silk,  bottom=bond+local+gbr,  color="gold",  label="Silk link")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Percentage")
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.legend(loc="lower left", fontsize=8)


def main():
    print(f"HETA_TIMES = {TIMES}")
    print("處理 WS ...")
    ws_fps = compute_or_load_fingerprints("ws",  "ws")
    print("處理 NW ...")
    nw_fps = compute_or_load_fingerprints("nws", "nw")

    labels = [p[1] for p in P_VALUES]
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    plot_stacked_bar(axes[0], ws_fps, labels,
                     "(a) WS small-world network",
                     "Random rewiring probability")
    plot_stacked_bar(axes[1], nw_fps, labels,
                     "(b) NW small-world network",
                     "Shortcut addition probability")
    plt.tight_layout()
    out = "fig8.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"產生：{os.path.abspath(out)}")


if __name__ == "__main__":
    main()
