# -*- coding: utf-8 -*-
"""
實驗 8 — Fig. 11：16 個實證網路指紋相關係數矩陣 + 階層分群 (clustermap)
對應論文 §3.3, Fig. 11
依賴：../../_cache/fp_real16.pkl（由 exp9 建立）
若快取不存在則自動呼叫 exp9 的相同邏輯重新計算。
輸出：fig11.png
"""
import os
import pickle
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.chdir(os.path.dirname(os.path.abspath(__file__)))
# 本資料夾位於 HETA/experiment/exp11_.../，因此 data/、_cache/ 要往上兩層
DATA_DIR = "../../data/net"
CACHE_DIR = "../../_cache"
os.makedirs(CACHE_DIR, exist_ok=True)

TIMES = int(os.environ.get("HETA_TIMES", 1000))

LABELS = [
    "rdgam", "women", "football", "jazz", "lesmis", "camp92",
    "prisonInter", "karate", "k-core", "ragusa16", "dolphins",
    "celegans", "s208", "ba_sfn", "leader", "florentine",
]
FILES = {
    "rdgam":       "14p.net",          "women":       "women.net",
    "football":    "football.net",     "jazz":        "jazz.net",
    "lesmis":      "lesmis.net",       "camp92":      "USA.net",
    "prisonInter": "prisonInter.net",  "karate":      "karate.net",
    "k-core":      "k-core.net",       "ragusa16":    "Ragusa16.net",
    "dolphins":    "dolphins.net",     "celegans":    "celegans.net",
    "s208":        "s208.net",         "ba_sfn":      "ba_sfn.net",
    "leader":      "leader.net",       "florentine":  "families.net",
}


def _ensure_cache():
    cache_path = os.path.join(CACHE_DIR, "fp_real16.pkl")
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            return np.array(pickle.load(f))

    print("  快取不存在，重新計算 16 實證網路指紋（耗時較長）...")
    import networkx as nx
    from heta import HETA as heta
    fps = []
    for i, name in enumerate(LABELS, 1):
        print(f"    [{i:>2}/16] {name}")
        G = nx.Graph(nx.read_pajek(os.path.join(DATA_DIR, FILES[name])))
        G.remove_edges_from(nx.selfloop_edges(G))
        G, _, _ = heta.bridge_or_bond(G, times=TIMES)
        fps.append(heta.fingerprint(G))
    with open(cache_path, "wb") as f:
        pickle.dump(fps, f)
    return np.array(fps)


def main():
    fps = _ensure_cache()
    corr = np.corrcoef(fps)  # 16×16

    cg = sns.clustermap(
        corr, method="average", metric="euclidean",
        cmap="RdBu_r", center=0, vmin=-0.8, vmax=0.8,
        xticklabels=LABELS, yticklabels=LABELS,
        figsize=(11, 11),
    )
    out = "fig11.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"產生：{os.path.abspath(out)}")

    row_order = cg.dendrogram_row.reordered_ind
    print("\n[驗收] 本次階層聚類順序（與論文 Fig. 11 對照）：")
    for i, idx in enumerate(row_order, 1):
        print(f"  {i:>2}. {LABELS[idx]}")
    print("論文 Fig. 11 順序：k-core, ba_sfn, florentine, s208, celegans, "
          "leader, karate, ragusa16, dolphins, prisonInter, camp92, jazz, "
          "lesmis, football, rdgam, women")


if __name__ == "__main__":
    main()
