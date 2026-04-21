# -*- coding: utf-8 -*-
"""
Table 1：16 個實證網路基本指標（非實驗，純資料規格）
對應論文：Huang et al. 2019 Physica A 536, 121027 — Table 1
列出每個網路的 n（節點數）、m（連結數）、density（密度）、
ASP（平均最短路徑）、CC（平均叢集係數），並與論文值比對。
輸出：
    table1.csv      — 機器可讀格式
    table1_compare.txt — 與論文 Table 1 的逐項比對
"""
import os
import sys

import networkx as nx
import pandas as pd

# 切換到本腳本所在資料夾，使輸出落在這裡
os.chdir(os.path.dirname(os.path.abspath(__file__)))

DATA_DIR = "../data/net"

# ──────────────────────────────────────────────────────────────
# 論文 Table 1 原始數據
# ──────────────────────────────────────────────────────────────
TABLE1 = {
    "rdgam":       {"n": 12,  "m": 28,   "density": 0.4242, "asp": 2.0909, "cc": 0.8306},
    "women":       {"n": 16,  "m": 50,   "density": 0.4167, "asp": 1.8167, "cc": 0.8182},
    "football":    {"n": 115, "m": 613,  "density": 0.0935, "asp": 2.5082, "cc": 0.4032},
    "jazz":        {"n": 198, "m": 2742, "density": 0.1406, "asp": 2.2350, "cc": 0.6175},
    "lesmis":      {"n": 77,  "m": 254,  "density": 0.0868, "asp": 2.6411, "cc": 0.5731},
    "camp92":      {"n": 18,  "m": 35,   "density": 0.2288, "asp": 2.6601, "cc": 0.5685},
    "prisonInter": {"n": 67,  "m": 142,  "density": 0.0642, "asp": 3.3546, "cc": 0.3099},
    "karate":      {"n": 34,  "m": 78,   "density": 0.1390, "asp": 2.4082, "cc": 0.5706},
    "leader":      {"n": 32,  "m": 80,   "density": 0.1613, "asp": 2.2964, "cc": 0.3266},
    "ragusa16":    {"n": 24,  "m": 68,   "density": 0.2464, "asp": 2.0942, "cc": 0.3427},
    "dolphins":    {"n": 62,  "m": 159,  "density": 0.0841, "asp": 3.3570, "cc": 0.2590},
    "florentine":  {"n": 15,  "m": 20,   "density": 0.1905, "asp": 2.4857, "cc": 0.1600},
    "celegans":    {"n": 297, "m": 2148, "density": 0.0489, "asp": 2.4553, "cc": 0.2924},
    "s208":        {"n": 122, "m": 189,  "density": 0.0256, "asp": 4.9278, "cc": 0.0591},
    "ba_sfn":      {"n": 100, "m": 196,  "density": 0.0396, "asp": 3.0313, "cc": 0.1343},
    "k-core":      {"n": 26,  "m": 31,   "density": 0.0954, "asp": 3.4708, "cc": 0.1607},
}

# repo 檔名對照
FILE_MAP = {
    "rdgam":       "14p.net",
    "women":       "women.net",
    "football":    "football.net",
    "jazz":        "jazz.net",
    "lesmis":      "lesmis.net",
    "camp92":      "USA.net",
    "prisonInter": "prisonInter.net",
    "karate":      "karate.net",
    "leader":      "leader.net",
    "ragusa16":    "Ragusa16.net",
    "dolphins":    "dolphins.net",
    "florentine":  "families.net",
    "celegans":    "celegans.net",
    "s208":        "s208.net",
    "ba_sfn":      "ba_sfn.net",
    "k-core":      "k-core.net",
}

DISCIPLINE = {
    **{k: "Social Sciences" for k in
       ["rdgam","women","football","jazz","lesmis","camp92","prisonInter",
        "karate","leader","ragusa16","dolphins","florentine"]},
    "celegans": "Biology",
    "s208":     "Electronics",
    "ba_sfn":   "Network Theory",
    "k-core":   "Network Theory",
}


def calc_metrics(G):
    """計算 Table 1 的五項指標"""
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    n = G.number_of_nodes()
    m = G.number_of_edges()
    density = nx.density(G)
    if nx.is_connected(G):
        asp = nx.average_shortest_path_length(G)
    else:
        largest_cc = max(nx.connected_components(G), key=len)
        sub = G.subgraph(largest_cc)
        asp = nx.average_shortest_path_length(sub)
    cc = nx.average_clustering(G)
    return n, m, density, asp, cc


def main():
    rows = []
    print("\n" + "=" * 130)
    print(f"{'Name':<12} {'File':<20} {'Discipline':<18} "
          f"{'N':>4}/{'N*':>4}  {'M':>5}/{'M*':>5}  "
          f"{'Dens':>7}/{'Dens*':>7}  {'ASP':>6}/{'ASP*':>6}  "
          f"{'CC':>6}/{'CC*':>6}  Match")
    print("-" * 130)

    for name, filename in FILE_MAP.items():
        path = os.path.join(DATA_DIR, filename)
        if not os.path.exists(path):
            print(f"{name:<12} {filename:<20}  MISSING")
            rows.append({"name": name, "file": filename, "status": "MISSING"})
            continue
        try:
            G = nx.read_pajek(path)
            n, m, density, asp, cc = calc_metrics(G)
            exp = TABLE1[name]
            match_n  = (n == exp["n"])
            match_m  = (m == exp["m"])
            match_d  = (abs(density - exp["density"]) < 0.002)
            match_a  = (abs(asp - exp["asp"]) < 0.05)
            match_c  = (abs(cc - exp["cc"]) < 0.005)
            all_match = match_n and match_m and match_d and match_a and match_c
            mk = lambda x: "✓" if x else "✗"
            print(f"{name:<12} {filename:<20} {DISCIPLINE[name]:<18} "
                  f"{n:>4}/{exp['n']:>4}{mk(match_n)}  "
                  f"{m:>5}/{exp['m']:>5}{mk(match_m)}  "
                  f"{density:>7.4f}/{exp['density']:>7.4f}{mk(match_d)}  "
                  f"{asp:>6.4f}/{exp['asp']:>6.4f}{mk(match_a)}  "
                  f"{cc:>6.4f}/{exp['cc']:>6.4f}{mk(match_c)}  "
                  f"{'ALL' if all_match else 'PARTIAL'}")
            rows.append({
                "name": name, "file": filename, "discipline": DISCIPLINE[name],
                "n": n, "n_paper": exp["n"],
                "m": m, "m_paper": exp["m"],
                "density": round(density, 4),  "density_paper": exp["density"],
                "asp": round(asp, 4),          "asp_paper": exp["asp"],
                "cc": round(cc, 4),            "cc_paper": exp["cc"],
                "all_match": all_match,
            })
        except Exception as e:
            print(f"{name:<12} {filename:<20}  ERROR: {e}")
            rows.append({"name": name, "file": filename, "status": f"ERROR: {e}"})

    print("-" * 130)
    n_match = sum(1 for r in rows if r.get("all_match"))
    print(f"完全匹配：{n_match}/{len(rows)}（women / Ragusa16 為已知差異版本，"
          "詳見 EXPERIMENT_MANUAL.txt 附錄 E）")
    print("=" * 130 + "\n")

    df = pd.DataFrame(rows)
    df.to_csv("table1.csv", index=False, encoding="utf-8-sig")
    print(f"已存至：{os.path.abspath('table1.csv')}")


if __name__ == "__main__":
    main()
