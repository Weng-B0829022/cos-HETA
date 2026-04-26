"""
HETA - Hierarchical Edge Type Analysis
=======================================
根據論文：Huang et al. (2019), Physica A 536, 121027
"Beyond bond links in complex networks: Local bridges, global bridges and silk links"
"""

import networkx as nx
import numpy as np
import random
import sys
import os
import time
import math

# ════════════════════════════════════════════════════════
# 1. 共同鄰居計算
# ════════════════════════════════════════════════════════

def get_kth_layer_neighbors(G, node, exclude, k):
    """回傳節點 node 的第 k 層鄰居集合（排除 exclude 節點及前幾層已訪問的節點）"""
    if k == 1:
        return set(G.neighbors(node)) - {exclude}

    visited = {node, exclude}
    frontier = {node}

    for _ in range(1, k):
        next_frontier = set()
        for n in frontier:
            for nb in G.neighbors(n):
                if nb not in visited:
                    next_frontier.add(nb)
                    visited.add(nb)
        frontier = next_frontier

    kth = set()
    for n in frontier:
        for nb in G.neighbors(n):
            if nb not in visited:
                kth.add(nb)
                visited.add(nb)
    return kth

def compute_common_neighbor_ratio(G, u, v, k):
    """計算邊 (u, v) 的第 k 層共同鄰居比例 R^k_{u,v}

    註：原論文分母採用 min(|A|,|B|)，本實作改用 sqrt(|A|·|B|)
    （幾何平均，等價於 Salton/cosine 相似度的歸一化方式）。
    """
    if k == 1:
        V1_u = set(G.neighbors(u)) - {v}
        V1_v = set(G.neighbors(v)) - {u}
        if not V1_u or not V1_v:
            return 0.0
        return len(V1_u & V1_v) / math.sqrt(len(V1_u) * len(V1_v))

    Vk_u  = get_kth_layer_neighbors(G, u, v, k)
    Vk_v  = get_kth_layer_neighbors(G, v, u, k)
    Vk1_u = get_kth_layer_neighbors(G, u, v, k - 1)
    Vk1_v = get_kth_layer_neighbors(G, v, u, k - 1)

    union = (Vk_u & Vk_v) | (Vk1_u & Vk_v) | (Vk_u & Vk1_v)
    if not union:
        return 0.0

    num = len(Vk_u & Vk_v) + len(Vk1_u & Vk_v) + len(Vk_u & Vk1_v)
    den = (math.sqrt(len(Vk_u)  * len(Vk_v)) +
           math.sqrt(len(Vk1_u) * len(Vk_v)) +
           math.sqrt(len(Vk_u)  * len(Vk1_v)))
    return num / den if den > 0 else 0.0

# ════════════════════════════════════════════════════════
# 2. kmax 計算
# ════════════════════════════════════════════════════════

def compute_kmax(G):
    """kmax = floor( 平均最短路徑長度 / 2 )，最小值為 1"""
    if not nx.is_connected(G):
        G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    try:
        avg_spl = nx.average_shortest_path_length(G)
    except Exception:
        return 1
    return max(1, int(avg_spl / 2))

# ════════════════════════════════════════════════════════
# 3. 隨機化網路生成（switching algorithm）
# ════════════════════════════════════════════════════════

def switching_randomize(G, Q=100):
    """使用 switching randomization 生成一個隨機化網路"""
    Gr = G.copy()
    edges = list(Gr.edges())
    m = len(edges)

    for _ in range(Q * m):
        if len(edges) < 2:
            break
        i, j = random.sample(range(len(edges)), 2)
        a, b = edges[i]
        c, d = edges[j]
        if a == d or c == b:
            continue
        if Gr.has_edge(a, d) or Gr.has_edge(c, b):
            continue
        Gr.remove_edge(a, b)
        Gr.remove_edge(c, d)
        Gr.add_edge(a, d)
        Gr.add_edge(c, b)
        edges[i] = (a, d)
        edges[j] = (c, b)
    return Gr

# ════════════════════════════════════════════════════════
# 4. 閾值計算
# ════════════════════════════════════════════════════════

def compute_external_threshold(G, k, n_random=1000):
    """計算第 k 層的外部閾值 T_E^k"""
    rand_ratios = []
    for _ in range(n_random):
        Gr = switching_randomize(G)
        for (u, v) in Gr.edges():
            rand_ratios.append(compute_common_neighbor_ratio(Gr, u, v, k))
    if not rand_ratios:
        return 0.0
    return float(np.mean(rand_ratios) + 2 * np.std(rand_ratios))

def compute_internal_threshold(candidate_ratios):
    """計算第 k 層的內部閾值 T_I^k"""
    if not candidate_ratios:
        return 0.0
    arr = np.array(candidate_ratios)
    return float(np.mean(arr) - np.std(arr))

# ════════════════════════════════════════════════════════
# 5. 主演算法：HETA
# ════════════════════════════════════════════════════════

def heta(G, n_random=1000):
    """
    對網路 G 執行 HETA，回傳每條邊的類型與層級。
    回傳: { (u, v): (type_string, layer_int) }
    """
    if not nx.is_connected(G):
        G = G.subgraph(max(nx.connected_components(G), key=len)).copy()

    kmax  = compute_kmax(G)
    edges = list(G.edges())

    link_types = {e: None  for e in edges}
    pass_flag  = {e: True  for e in edges}

    # Step 3.1：識別 silk links
    for (u, v) in edges:
        if G.degree(u) == 1 or G.degree(v) == 1:
            link_types[(u, v)] = ('silk', 0)
            pass_flag[(u, v)]  = False

    # Step 2：計算各層外部閾值
    ext_thresholds = {}
    for k in range(1, kmax + 1):
        ext_thresholds[k] = compute_external_threshold(G, k, n_random)

    # Step 3.2：逐層識別 bond links 與 local bridges
    for k in range(1, kmax + 1):
        T_E = ext_thresholds[k]
        candidate_edges  = []
        candidate_ratios = []

        for (u, v) in edges:
            if not pass_flag[(u, v)]:
                continue
            r = compute_common_neighbor_ratio(G, u, v, k)
            if r >= T_E:
                link_types[(u, v)] = ('bond', k)
                pass_flag[(u, v)]  = False
            else:
                candidate_edges.append((u, v))
                candidate_ratios.append(r)

        if candidate_edges:
            T_I = compute_internal_threshold(candidate_ratios)
            for (u, v), r in zip(candidate_edges, candidate_ratios):
                if r > T_I:
                    link_types[(u, v)] = ('local_bridge', k)
                    pass_flag[(u, v)]  = False

    # Step 3.3：剩餘未分類 = global bridges
    for (u, v) in edges:
        if pass_flag[(u, v)]:
            link_types[(u, v)] = ('global_bridge', 0)

    return link_types, kmax

# ════════════════════════════════════════════════════════
# 6. 輔助：統計與輸出
# ════════════════════════════════════════════════════════

def summarize(link_types):
    counts = {'silk': 0, 'bond': 0, 'local': 0, 'global': 0}
    for (t_str, layer) in link_types.values():
        if t_str == 'bond':           counts['bond'] += 1
        elif t_str == 'local_bridge': counts['local'] += 1
        elif t_str == 'global_bridge':counts['global'] += 1
        elif t_str == 'silk':         counts['silk'] += 1
        
    total = sum(counts.values())
    return {
        k: {'count': v, 'pct': round(v / total * 100, 1) if total else 0.0}
        for k, v in counts.items()
    }

def read_edge_list(filepath):
    G = nx.Graph()
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                G.add_edge(parts[0], parts[1])
    return G

# ════════════════════════════════════════════════════════
# 7. 主程式
# ════════════════════════════════════════════════════════

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("用法: python3 HETA.py <edge_file>")
        sys.exit(1)

    filepath = sys.argv[1]
    n_random = 1000 # 依論文預設為 1000
    
    if not os.path.exists(filepath):
        print(f"找不到檔案: {filepath}")
        sys.exit(1)

    # 讀取網路圖
    G = read_edge_list(filepath)
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()

    # 執行演算法並計時
    start_time = time.time()
    link_types, kmax = heta(G, n_random=n_random)
    exec_time = time.time() - start_time

    # 統計結果
    stats = summarize(link_types)

    # 準備輸出內容
    out_lines = []
    out_lines.append("HETA 分析結果")
    out_lines.append("=" * 80)
    out_lines.append(f"節點數     : {num_nodes}")
    out_lines.append(f"邊數       : {num_edges}")
    out_lines.append(f"kmax       : {kmax}")
    out_lines.append(f"n_random   : {n_random}")
    out_lines.append(f"執行時間   : {exec_time:.3f} 秒")
    out_lines.append("=" * 80)
    out_lines.append("")
    out_lines.append("分類統計")
    out_lines.append("-" * 80)
    out_lines.append(f"  silk link      : {stats['silk']['count']:>4} 條  ({stats['silk']['pct']:>5}%)")
    out_lines.append(f"  bond link      : {stats['bond']['count']:>4} 條  ({stats['bond']['pct']:>5}%)")
    out_lines.append(f"  local bridge   : {stats['local']['count']:>4} 條  ({stats['local']['pct']:>5}%)")
    out_lines.append(f"  global bridge  : {stats['global']['count']:>4} 條  ({stats['global']['pct']:>5}%)")
    out_lines.append("-" * 80)
    out_lines.append("")
    out_lines.append("逐邊分類結果")
    out_lines.append("-" * 80)
    out_lines.append(f"  {'u':<6}  {'v':<6}  {'type':<22}  {'layer':<5}")
    out_lines.append(f"  {'──────':<6}  {'──────':<6}  {'──────────────────────':<22}  {'─────':<5}")

    for (u, v), (t_str, layer) in link_types.items():
        # 依照範例格式化 type 名稱
        display_type = f"local_bridge_k{layer}" if t_str == 'local_bridge' else t_str
        out_lines.append(f"  {u:<6}  {v:<6}  {display_type:<22}  {layer:<5}")

    output_text = "\n".join(out_lines) + "\n"

    # 建立目錄並輸出檔案
    out_dir = "python_result"
    os.makedirs(out_dir, exist_ok=True)
    
    filename = os.path.basename(filepath)
    file_prefix = os.path.splitext(filename)[0]
    out_path = os.path.join(out_dir, f"{file_prefix}_python_result.txt")

    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(output_text)

    print(f"分析完成！結果已儲存至: {out_path}")