# -*- coding: utf-8 -*-
"""
實驗 6 — Fig. 9：16 個實證網路的指紋（fingerprint）堆疊長條圖
對應論文 §3.3, Fig. 9
快取：../_cache/fp_real16.pkl（給 exp11 重複利用）
輸出：fig9.png
"""
import os
import pickle
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from heta import HETA as heta

os.chdir(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = "../data/net"
CACHE_DIR = "../_cache"
os.makedirs(CACHE_DIR, exist_ok=True)

TIMES = int(os.environ.get("HETA_TIMES", 1000))

# 順序 = Fig. 9（依 bond 比例由高至低）
NETWORKS = [
    ("rdgam",       "14p.net"),
    ("women",       "women.net"),
    ("football",    "football.net"),
    ("jazz",        "jazz.net"),
    ("lesmis",      "lesmis.net"),
    ("camp92",      "USA.net"),
    ("prisonInter", "prisonInter.net"),
    ("karate",      "karate.net"),
    ("k-core",      "k-core.net"),
    ("ragusa16",    "Ragusa16.net"),
    ("dolphins",    "dolphins.net"),
    ("celegans",    "celegans.net"),
    ("s208",        "s208.net"),
    ("ba_sfn",      "ba_sfn.net"),
    ("leader",      "leader.net"),
    ("florentine",  "families.net"),
]


def compute_or_load():
    cache_path = os.path.join(CACHE_DIR, "fp_real16.pkl")
    if os.path.exists(cache_path):
        print(f"  從快取載入：{cache_path}")
        with open(cache_path, "rb") as f:
            return pickle.load(f)

    fps = []
    for i, (name, filename) in enumerate(NETWORKS, 1):
        print(f"  [{i:>2}/16] 處理 {name} ...")
        G = nx.Graph(nx.read_pajek(os.path.join(DATA_DIR, filename)))
        G.remove_edges_from(nx.selfloop_edges(G))
        G, _, _ = heta.bridge_or_bond(G, times=TIMES)
        fps.append(heta.fingerprint(G))
    with open(cache_path, "wb") as f:
        pickle.dump(fps, f)
    print(f"  快取已寫入：{cache_path}")
    return fps


def main():
    print(f"HETA_TIMES = {TIMES}")
    fps = np.array(compute_or_load())
    labels = [n[0] for n in NETWORKS]
    bond, local, gbr, silk = fps[:, 0], fps[:, 1], fps[:, 2], fps[:, 3]
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(14, 5.5))
    ax.bar(x, bond,  color="blue",  label="Bond link")
    ax.bar(x, local, bottom=bond,           color="red",   label="Local bridge")
    ax.bar(x, gbr,   bottom=bond+local,     color="green", label="Global bridge")
    ax.bar(x, silk,  bottom=bond+local+gbr, color="gold",  label="Silk link")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Percentage")
    ax.set_xlabel("Network")
    ax.set_title("Fig. 9: Network fingerprints of 16 real-world networks")
    ax.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    out = "fig9.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"產生：{os.path.abspath(out)}")

    print("\n[驗收] 16 個網路指紋（bond / local / global / silk）：")
    for name, fp in zip(labels, fps):
        print(f"  {name:<13} {fp[0]:.2f}  {fp[1]:.2f}  {fp[2]:.2f}  {fp[3]:.2f}")


if __name__ == "__main__":
    main()
