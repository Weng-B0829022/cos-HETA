# -*- coding: utf-8 -*-
"""
實驗 7 — Fig. 10：WS / NW 16 個指紋兩兩 Pearson 相關係數矩陣
對應論文 §3.3, Fig. 10a (WS) 與 Fig. 10b (NW)
依賴：../../_cache/fp_ws.pkl 與 fp_nw.pkl（由 exp8 建立）
若快取不存在則自動呼叫 exp8 的相同邏輯重新計算。
輸出：fig10.png
"""
import os
import pickle
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 本腳本位於 HETA/src/exp10_smallworld_correlation/
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

P_STR = ["0.0", "0.001", "0.002", "0.004", "0.008", "0.016", "0.032",
         "0.064", "0.128", "0.256", "0.384", "0.512", "0.640", "0.768",
         "0.896", "1.0"]


def _ensure_cache(prefix, label):
    """若快取不存在，重跑指紋計算（與 exp8 一致）"""
    cache_path = os.path.join(CACHE_DIR, f"fp_{label}.pkl")
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            return np.array(pickle.load(f))

    print(f"  快取不存在，重新計算 {label} 指紋（耗時較長）...")
    import networkx as nx
    from heta import HETA as heta
    fps = []
    for p_str in P_STR:
        print(f"    {label} p={p_str}")
        G = nx.Graph(nx.read_pajek(
            os.path.join(DATA_DIR, f"{prefix}_swn_{p_str}.net")))
        G.remove_edges_from(nx.selfloop_edges(G))
        G, _, _ = heta.bridge_or_bond(G, times=TIMES)
        fps.append(heta.fingerprint(G))
    with open(cache_path, "wb") as f:
        pickle.dump(fps, f)
    return np.array(fps)


def main():
    ws = _ensure_cache("ws",  "ws")   # shape (16, 4)
    nw = _ensure_cache("nws", "nw")
    ws_corr = np.corrcoef(ws)         # (16, 16)
    nw_corr = np.corrcoef(nw)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    sns.heatmap(ws_corr, ax=axes[0], cmap="RdBu_r", center=0,
                vmin=-0.8, vmax=0.8, square=True,
                xticklabels=P_STR, yticklabels=P_STR,
                cbar_kws={"shrink": 0.8})
    axes[0].set_title("(a) WS small-world network")
    axes[0].set_xlabel("Random rewiring probability")
    axes[0].set_ylabel("Random rewiring probability")

    sns.heatmap(nw_corr, ax=axes[1], cmap="Reds", vmin=0.75, vmax=1.00,
                square=True, xticklabels=P_STR, yticklabels=P_STR,
                cbar_kws={"shrink": 0.8})
    axes[1].set_title("(b) NW small-world network")
    axes[1].set_xlabel("Shortcut addition probability")
    axes[1].set_ylabel("Shortcut addition probability")

    plt.tight_layout()
    out = "fig10.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"產生：{os.path.abspath(out)}")

    print("\n[驗收] WS corr 子矩陣對角項與最低項：")
    print(f"  WS  min corr ≈ {ws_corr.min():+.3f}")
    print(f"  NW  min corr ≈ {nw_corr.min():+.3f}（應 ≥ 0.7）")


if __name__ == "__main__":
    main()
