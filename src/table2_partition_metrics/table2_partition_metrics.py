# -*- coding: utf-8 -*-
"""
實驗 10 — Table 2：16 個實證網路社群分割的 6 項量化指標
對應論文 §3.4, Table 2

6 項指標
--------
1. group_num         群數
2. independent_num   獨立節點數（不屬於任何群）
3. total_in_groups   各群節點數總和
4. density (avg±std) 群密度的平均與標準差
5. within  (avg±std) 群內連結 R¹ 的加權平均（以群節點數加權）
6. between (avg±std) 群間連結 R¹ 的平均

輸出
----
table2.csv          機器可讀
table2_summary.txt  人類可讀（對齊論文表格樣式）
"""
import os
import numpy as np
import networkx as nx
import pandas as pd

from heta import HETA as heta

# 本腳本位於 HETA/src/table2_partition_metrics/
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

TIMES = int(os.environ.get("HETA_TIMES", 1000))

# 順序 = Table 2 排列（與 Fig. 9 / Fig. 7 一致）
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


def first_layer_cn_ratio(G, u, v):
    nbr_u = set(G.neighbors(u)) - {v}
    nbr_v = set(G.neighbors(v)) - {u}
    if not nbr_u or not nbr_v:
        return 0.0
    return len(nbr_u & nbr_v) / min(len(nbr_u), len(nbr_v))


def flatten_communities(communities):
    """將 heta.community_sorting() 的回傳結構展平為 [set, set, ...]

    heta.community_sorting(G) 官方回傳是 4-tuple：
        (nodes_sorted, h_tree_dic, flat_cdic, kmax)
    其中 flat_cdic = {leader_node: [member_node, ...]}，
    正是 Table 2 所要的「社群列表」。早期版本將整個 tuple 當作
    通用容器遞迴展平，會把 flat_cdic 的 value（節點字串 list）
    錯誤地拆成 len==1 的 singleton set，造成後續 `len(c)>1` 過濾
    後全部為空、導致 table2.csv 每一欄皆為 0。

    正確作法：若輸入明顯是 community_sorting 的 4-tuple，直接取
    第 3 個元素 flat_cdic；否則才退回原本的保守展平邏輯。
    """
    # 專用分支：community_sorting 的 4-tuple → 取 flat_cdic
    if (isinstance(communities, tuple) and len(communities) == 4
            and isinstance(communities[2], dict)):
        flat_cdic = communities[2]
        return [set(str(n) for n in nodes) for nodes in flat_cdic.values()]

    # 通用分支（保留舊邏輯作為 fallback）
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
                out.append({str(item)})
    return out


def compute_metrics(G, communities):
    """計算 Table 2 的 6 項指標"""
    # 過濾單節點群
    communities = [c for c in communities if len(c) > 1]

    group_num = len(communities)

    in_group = set()
    for comm in communities:
        in_group.update(str(n) for n in comm)
    independent = [str(n) for n in G.nodes() if str(n) not in in_group]
    independent_num = len(independent)

    total_in_groups = sum(len(c) for c in communities)

    # 群密度
    densities = []
    for comm in communities:
        nodes = [n for n in (str(x) for x in comm) if n in G.nodes()]
        sub = G.subgraph(nodes)
        if sub.number_of_nodes() > 1:
            densities.append(nx.density(sub))
    density_avg = float(np.mean(densities)) if densities else 0.0
    density_std = float(np.std(densities))  if densities else 0.0

    # 群內連結 R¹（以群節點數為權重）
    within_r1, within_w = [], []
    for comm in communities:
        comm_set = set(str(n) for n in comm)
        sub_edges = [(u, v) for u, v in G.edges()
                     if u in comm_set and v in comm_set]
        if sub_edges:
            r1s = [first_layer_cn_ratio(G, u, v) for u, v in sub_edges]
            within_r1.append(np.mean(r1s))
            within_w.append(len(comm_set))
    if within_r1:
        within_avg = float(np.average(within_r1, weights=within_w))
        within_std = float(np.sqrt(np.average(
            [(x - within_avg) ** 2 for x in within_r1], weights=within_w)))
    else:
        within_avg, within_std = 0.0, 0.0

    # 群間連結 R¹
    comm_lookup = {}
    for i, comm in enumerate(communities):
        for n in comm:
            comm_lookup[str(n)] = i
    between_edges = [(u, v) for u, v in G.edges()
                     if u in comm_lookup and v in comm_lookup
                     and comm_lookup[u] != comm_lookup[v]]
    if between_edges:
        r1s = [first_layer_cn_ratio(G, u, v) for u, v in between_edges]
        between_avg = float(np.mean(r1s))
        between_std = float(np.std(r1s))
    else:
        between_avg, between_std = 0.0, 0.0

    return {
        "group_num":       group_num,
        "independent_num": independent_num,
        "total_in_groups": total_in_groups,
        "density_avg":     round(density_avg, 4),
        "density_std":     round(density_std, 4),
        "within_avg":      round(within_avg, 4),
        "within_std":      round(within_std, 4),
        "between_avg":     round(between_avg, 4),
        "between_std":     round(between_std, 4),
    }


def main():
    print(f"HETA_TIMES = {TIMES}")
    results = []
    for i, (name, filename) in enumerate(NETWORKS, 1):
        print(f"  [{i:>2}/16] 處理 {name} ...")
        G = nx.Graph(nx.read_pajek(os.path.join(DATA_DIR, filename)))
        G.remove_edges_from(nx.selfloop_edges(G))
        G, _, _ = heta.bridge_or_bond(G, times=TIMES)
        try:
            communities = heta.community_sorting(G)
        except Exception as e:
            print(f"    community_sorting 失敗：{e}")
            communities = []
        comm_list = flatten_communities(communities)

        m = compute_metrics(G, comm_list)
        m["name"] = name
        results.append(m)

    df = pd.DataFrame(results)
    df = df[["name", "group_num", "independent_num", "total_in_groups",
             "density_avg", "density_std",
             "within_avg",  "within_std",
             "between_avg", "between_std"]]
    df.to_csv("table2.csv", index=False, encoding="utf-8-sig")

    # 同步輸出對齊版的 txt 摘要
    lines = []
    lines.append("Table 2 — 16 networks community partition metrics "
                 f"(HETA_TIMES={TIMES})")
    header = (f"{'Network':<13} {'Grp#':>4} {'Indep':>5} {'Total':>5}  "
              f"{'Density':>14}  {'Within R1':>14}  {'Between R1':>14}")
    lines.append(header)
    lines.append("-" * len(header))
    for _, r in df.iterrows():
        lines.append(
            f"{r['name']:<13} {int(r['group_num']):>4d} "
            f"{int(r['independent_num']):>5d} {int(r['total_in_groups']):>5d}  "
            f"{r['density_avg']:.2f} ({r['density_std']:.2f})    "
            f"{r['within_avg']:.2f} ({r['within_std']:.2f})    "
            f"{r['between_avg']:.2f} ({r['between_std']:.2f})"
        )
    text = "\n".join(lines) + "\n"
    print("\n" + text)
    with open("table2_summary.txt", "w", encoding="utf-8") as f:
        f.write(text)

    print(f"已輸出：{os.path.abspath('table2.csv')}")
    print(f"已輸出：{os.path.abspath('table2_summary.txt')}")


if __name__ == "__main__":
    main()
