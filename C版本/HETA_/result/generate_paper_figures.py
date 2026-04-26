#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_paper_figures.py
==========================
從 HETA 結果檔（*_result.txt 與 summary.csv）一次產生論文 Fig. 4 ~ Fig. 12 的圖。

使用方式：
    把這個檔案放進跟所有 *_result.txt 同一個資料夾，然後執行
        python3 generate_paper_figures.py

    或指定資料夾與輸出資料夾：
        python3 generate_paper_figures.py --result-dir ./result --out-dir ./figures

需求套件：
    matplotlib, numpy, networkx, scipy

對應論文：
    Huang, Chin, Fu, Tsai (2019). "Beyond bond links in complex networks:
    Local bridges, global bridges and silk links." Physica A 536, 121027.
"""

import argparse
import re
from pathlib import Path

# 強制使用 non-GUI backend，避免在沒有圖形介面（或 Qt 異常）的環境啟動失敗
# 必須在 `import matplotlib.pyplot` 之前呼叫
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform


# =============================================================================
# 全域設定
# =============================================================================

# 論文掃描的 16 個機率值
P_VALUES = ['0.0', '0.001', '0.002', '0.004', '0.008', '0.016', '0.032',
            '0.064', '0.128', '0.256', '0.384', '0.512', '0.640', '0.768',
            '0.896', '1.0']

# 論文 Table 1 代號 → 你的檔名
NETWORK_MAP = [
    ('rdgam',       'rdgam'),
    ('women',       'women'),
    ('football',    'football'),
    ('jazz',        'jazz'),
    ('lesmis',      'lesmis'),
    ('camp92',      'camp92'),
    ('prisonInter', 'prisonInter'),
    ('karate',      'karate'),
    ('leader',      'leader'),
    ('ragusa16',    'Ragusa16'),
    ('dolphins',    'dolphins'),
    ('florentine',  'florentine'),
    ('celegans',    'celegans'),
    ('s208',        's208'),
    ('ba_sfn',      'ba_sfn'),
    ('k-core',      'k-core'),
]

# 連結類型對應顏色（與論文圖例一致：純藍/純紅/純綠/純黃）
EDGE_COLORS = {
    'bond':          '#0000FF',   # 純藍
    'local_bridge':  '#FF0000',   # 純紅
    'global_bridge': '#008000',   # 純綠（深綠，跟論文一致）
    'silk':          '#E6C700',   # 金黃（純黃在白底會看不到，用稍深的金黃近似論文圖例黃色）
}

# 社群上色用調色盤
COMMUNITY_PALETTE = ['#1f77b4', '#2ca02c', '#ff7f0e', '#d62728', '#9467bd',
                     '#8c564b', '#e377c2', '#bcbd22', '#17becf', '#7f7f7f']

KEYS = ['bond', 'local_bridge', 'global_bridge', 'silk']
LEGEND_LABELS = ['Bond link', 'Local bridge', 'Global bridge', 'Silk link']


# =============================================================================
# 解析 *_result.txt
# =============================================================================

def parse_result_file(filepath):
    """
    解析單一 HETA 結果檔，回傳 dict：
        meta       : nodes, edges, kmax, n_random
        stats      : 各類型計數與百分比
        edges      : list of {u, v, type, layer, r1}
    """
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # meta
    meta = {}
    for key, pat in [
        ('nodes',    r'節點數\s*:\s*(\d+)'),
        ('edges',    r'邊數\s*:\s*(\d+)'),
        ('kmax',     r'kmax\s*:\s*(\d+)'),
        ('n_random', r'n_random\s*:\s*(\d+)'),
    ]:
        m = re.search(pat, content)
        if m:
            meta[key] = int(m.group(1))

    # 分類統計
    stats = {}
    for key, label in [
        ('silk',          'silk link'),
        ('bond',          'bond link'),
        ('local_bridge',  'local bridge'),
        ('global_bridge', 'global bridge'),
    ]:
        m = re.search(rf'{label}\s*:\s*(\d+)\s*條\s*\(\s*([\d.]+)%\)', content)
        if m:
            stats[key] = {'count': int(m.group(1)), 'pct': float(m.group(2))}

    # 逐邊資料（含 R^1）
    edges = []
    block_match = re.search(r'逐邊分類結果.*?\n(.*?)(?:\n-{5,}\s*$|\Z)',
                            content, re.DOTALL)
    if block_match:
        for line in block_match.group(1).split('\n'):
            line = line.strip()
            if not line or line.startswith('-') or line.startswith('u') \
                    or 'type' in line.lower():
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    edges.append({
                        'u':     int(parts[0]),
                        'v':     int(parts[1]),
                        'type':  parts[2],
                        'layer': int(parts[3]),
                        'r1':    float(parts[4]) if len(parts) >= 5 else None,
                    })
                except ValueError:
                    pass

    return {'meta': meta, 'stats': stats, 'edges': edges}


def normalize_type(t):
    """把 raw type 字串歸到 4 類之一"""
    if t == 'silk':          return 'silk'
    if t == 'bond':          return 'bond'
    if 'local_bridge' in t:  return 'local_bridge'
    if t == 'global_bridge': return 'global_bridge'
    return None


def get_r1_list(filepath):
    """回傳所有邊的 R^1 值"""
    r = parse_result_file(filepath)
    return [e['r1'] for e in r['edges'] if e['r1'] is not None]


def get_pct_vec(filepath):
    """回傳 4 維百分比向量 [bond, local_bridge, global_bridge, silk]"""
    r = parse_result_file(filepath)
    edges = r['edges']
    if not edges:
        return None
    counts = dict.fromkeys(KEYS, 0)
    for e in edges:
        nt = normalize_type(e['type'])
        if nt:
            counts[nt] += 1
    n = len(edges)
    return np.array([counts[k] / n for k in KEYS])


def get_pct_dict(filepath):
    vec = get_pct_vec(filepath)
    if vec is None:
        return None
    return {k: vec[i] for i, k in enumerate(KEYS)}


def build_graph(filepath):
    """從結果檔建 networkx Graph 並回傳 (G, type_map)
    type_map[(u_sorted, v_sorted)] = (normalized_type, layer)
    """
    r = parse_result_file(filepath)
    G = nx.Graph()
    type_map = {}
    for e in r['edges']:
        u, v = e['u'], e['v']
        G.add_edge(u, v)
        nt = normalize_type(e['type']) or 'bond'
        layer = e['layer']
        # local_bridge 的 layer 從 type 字串取（local_bridge_k1 → 1）
        if nt == 'local_bridge':
            try:
                layer = int(e['type'].split('_k')[-1])
            except ValueError:
                pass
        type_map[tuple(sorted([u, v]))] = (nt, layer)
    return G, type_map


# =============================================================================
# 通用畫圖工具
# =============================================================================

def make_legend_handles():
    return [Line2D([0], [0], color=EDGE_COLORS[k], lw=2, label=lbl)
            for k, lbl in zip(KEYS, LEGEND_LABELS)]


def safe_boxplot(ax, data, labels, **kwargs):
    """相容新舊 matplotlib：3.9+ 用 tick_labels、3.8- 用 labels"""
    try:
        return ax.boxplot(data, tick_labels=labels, **kwargs)
    except TypeError:
        return ax.boxplot(data, labels=labels, **kwargs)


def plot_stacked(ax, data_list, labels, title, xlabel, show_legend=True):
    """堆疊長條圖。論文配色：bond 藍、local_bridge 紅、global_bridge 綠、silk 黃。"""
    n = len(data_list)
    x = np.arange(n)
    bottoms = np.zeros(n)
    for key, label in zip(KEYS, LEGEND_LABELS):
        vals = np.array([d.get(key, 0) for d in data_list])
        ax.bar(x, vals, bottom=bottoms, color=EDGE_COLORS[key],
               label=label, width=0.85, edgecolor='white', linewidth=0.3)
        bottoms += vals
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylim(0, 1.0)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylabel('Percentage', fontsize=11)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_title(title, fontsize=12)
    if show_legend:
        # 論文風格：legend 放在繪圖區頂部一橫排，不遮住長條也不蓋住 title
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.10),
                  ncol=4, fontsize=10, framealpha=0.95,
                  handlelength=1.5, columnspacing=1.2)


def correlation_matrix(vectors):
    """計算 N×N Pearson 相關矩陣（向量為 4 維）"""
    n = len(vectors)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            v1, v2 = vectors[i], vectors[j]
            if np.std(v1) == 0 or np.std(v2) == 0:
                M[i, j] = 1.0 if np.allclose(v1, v2) else 0.0
            else:
                M[i, j] = np.corrcoef(v1, v2)[0, 1]
    return M


def draw_network(ax, G, type_map, title, layout='spring', seed=42):
    """畫單一網路：邊按類型上色、節點黑點"""
    if len(G.nodes()) == 0:
        ax.set_title(f"{title} (empty)")
        ax.axis('off')
        return

    if layout == 'spring':
        pos = nx.spring_layout(G, seed=seed,
                               k=1 / np.sqrt(len(G.nodes())))
    elif layout == 'circular_by_id':
        # 按節點編號等距排在圓周上（適合 ring lattice）
        sorted_nodes = sorted(G.nodes())
        n = len(sorted_nodes)
        pos = {v: (np.cos(2 * np.pi * i / n), np.sin(2 * np.pi * i / n))
               for i, v in enumerate(sorted_nodes)}
    else:
        pos = nx.spring_layout(G, seed=seed)

    # bond 在底層，bridges 在上層
    for et in ['bond', 'local_bridge', 'global_bridge', 'silk']:
        edge_list = [(u, v) for (u, v) in G.edges()
                     if type_map.get(tuple(sorted([u, v])), ('bond', 0))[0] == et]
        if edge_list:
            lw = 0.6 if et == 'bond' else 1.2
            alpha = 0.55 if et == 'bond' else 0.9
            nx.draw_networkx_edges(G, pos, edgelist=edge_list,
                                   edge_color=EDGE_COLORS[et],
                                   width=lw, alpha=alpha, ax=ax)

    node_size = max(3, min(40, 200 / np.sqrt(len(G.nodes()))))
    nx.draw_networkx_nodes(G, pos, node_size=node_size,
                           node_color='black', ax=ax)
    ax.set_title(title, fontsize=10)
    ax.axis('off')


# =============================================================================
# Algorithm 2: 階層式社群分割
# =============================================================================

def hierarchical_partition(G, type_map):
    """回傳 dict: node -> community_id（outliers = -1）"""
    G_work = G.copy()
    outliers = set()

    # Step a：移除 silk
    for u, v in [k for k, val in type_map.items() if val[0] == 'silk']:
        if G_work.has_edge(u, v):
            G_work.remove_edge(u, v)
    for n in list(G_work.nodes()):
        if G_work.degree(n) == 0:
            outliers.add(n); G_work.remove_node(n)

    # Step b：移除 global bridges
    for u, v in [k for k, val in type_map.items() if val[0] == 'global_bridge']:
        if G_work.has_edge(u, v):
            G_work.remove_edge(u, v)
    for n in list(G_work.nodes()):
        if G_work.degree(n) == 0:
            outliers.add(n); G_work.remove_node(n)

    # Step c, d：對每個連通分量遞迴
    max_layer = max([v[1] for v in type_map.values() if v[0] == 'local_bridge'],
                    default=0)

    def recurse(nodes, layer):
        sub = G_work.subgraph(nodes).copy()
        if layer < 1 or sub.number_of_edges() == 0:
            return [set(nodes)]
        lb_edges = [(u, v) for (u, v) in sub.edges()
                    if type_map.get(tuple(sorted([u, v])), ('', 0))
                       == ('local_bridge', layer)]
        for u, v in lb_edges:
            sub.remove_edge(u, v)
        result = []
        for comp in nx.connected_components(sub):
            if len(comp) == 1:
                outliers.update(comp)
            else:
                result.extend(recurse(comp, layer - 1))
        return result

    final = []
    for comp in nx.connected_components(G_work):
        if len(comp) == 1:
            outliers.update(comp)
        else:
            final.extend(recurse(comp, max_layer))

    node2comm = {}
    for cid, comm in enumerate(final):
        for n in comm:
            node2comm[n] = cid
    for n in outliers:
        node2comm[n] = -1
    for n in G.nodes():
        node2comm.setdefault(n, -1)
    return node2comm


def draw_partition(ax, G, type_map, title):
    if len(G.nodes()) == 0:
        ax.set_title(title + " (empty)"); ax.axis('off'); return

    node2comm = hierarchical_partition(G, type_map)
    pos = nx.spring_layout(G, seed=7, k=1.5 / np.sqrt(len(G.nodes())))

    node_colors = ['lightgray' if node2comm[n] == -1
                   else COMMUNITY_PALETTE[node2comm[n] % len(COMMUNITY_PALETTE)]
                   for n in G.nodes()]

    for et in ['bond', 'local_bridge', 'global_bridge', 'silk']:
        edge_list = [(u, v) for (u, v) in G.edges()
                     if type_map.get(tuple(sorted([u, v])), ('bond', 0))[0] == et]
        if edge_list:
            lw = 1.5 if et != 'bond' else 1.0
            nx.draw_networkx_edges(G, pos, edgelist=edge_list,
                                   edge_color=EDGE_COLORS[et],
                                   width=lw, alpha=0.85, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                           node_size=180, edgecolors='black',
                           linewidths=0.6, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)
    ax.set_title(title, fontsize=11)
    ax.axis('off')


# =============================================================================
# 9 張圖
# =============================================================================

def fig4(result_dir, out_dir):
    """Fig. 4: WS / NW small-world R^1 箱型圖"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, prefix, title, xlabel, color in [
        (axes[0], 'ws_swn',  '(a) WS small-world networks',
         'Random rewiring probability',  'lightblue'),
        (axes[1], 'nws_swn', '(b) NW small-world networks',
         'Shortcut addition probability', 'lightcoral'),
    ]:
        data, labs = [], []
        for p in P_VALUES:
            fp = result_dir / f'{prefix}_{p}_result.txt'
            if fp.exists():
                data.append(get_r1_list(fp)); labs.append(p)
        safe_boxplot(ax, data, labs, patch_artist=True,
                     boxprops=dict(facecolor=color),
                     medianprops=dict(color='red'))
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel('First-layer common neighbor ratio $R^1$', fontsize=11)
        ax.set_ylim(-0.05, 1.05)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / 'Fig4_R1_boxplot_smallworld.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print('  Fig. 4 saved')


def fig5(result_dir, out_dir):
    """Fig. 5: 16 個真實網路 R^1 箱型圖（橫軸順序按使用者指定，固定）"""
    # 固定的橫軸順序
    fixed_order = ['rdgam', 'women', 'lesmis', 'jazz', 'karate', 'camp92',
                   'ragusa16', 'football', 'dolphins', 'celegans', 'leader',
                   'prisonInter', 'k-core', 'florentine', 'ba_sfn', 's208']
    name_to_file = dict(NETWORK_MAP)

    data = []
    for paper in fixed_order:
        fname = name_to_file.get(paper)
        if fname is None:
            continue
        fp = result_dir / f'{fname}_result.txt'
        if fp.exists():
            r1 = get_r1_list(fp)
            data.append((paper, r1))

    fig, ax = plt.subplots(figsize=(13, 5.5))
    safe_boxplot(ax, [d[1] for d in data], [d[0] for d in data],
                 patch_artist=True,
                 boxprops=dict(facecolor='lightyellow'),
                 medianprops=dict(color='red'))
    ax.set_title('Fig. 5 — First-layer common neighbor ratio distributions for the 16 networks',
                 fontsize=12)
    ax.set_xlabel('Network', fontsize=11)
    ax.set_ylabel('First-layer common neighbor ratio $R^1$', fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.tick_params(axis='x', rotation=45)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / 'Fig5_R1_boxplot_16networks.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print('  Fig. 5 saved')


def fig6(result_dir, out_dir):
    """Fig. 6: WS / NW 三個 p 的視覺化（環狀排列；bond 沿圓周畫弧形強化環骨架）"""
    from matplotlib.patches import Arc

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    settings = [
        ('ws_swn_0.004',  '(a) WS, p = 0.004'),
        ('ws_swn_0.064',  '(b) WS, p = 0.064'),
        ('ws_swn_0.256',  '(c) WS, p = 0.256'),
        ('nws_swn_0.004', '(d) NW, p = 0.004'),
        ('nws_swn_0.064', '(e) NW, p = 0.064'),
        ('nws_swn_0.256', '(f) NW, p = 0.256'),
    ]
    for idx, (fname, title) in enumerate(settings):
        ax = axes[idx // 3, idx % 3]
        fp = result_dir / f'{fname}_result.txt'
        if fp.exists():
            G, tm = build_graph(fp)
            sorted_nodes = sorted(G.nodes())
            n = len(sorted_nodes)
            # 節點角度（單位：度，順時針從 12 點方向開始）
            node_angle_deg = {}
            pos = {}
            for i, v in enumerate(sorted_nodes):
                # 從 12 點方向（90°）順時針排
                deg = 90.0 - 360.0 * i / n
                rad = np.deg2rad(deg)
                pos[v] = (np.cos(rad), np.sin(rad))
                node_angle_deg[v] = deg

            # 1) 先畫 bond：用沿圓周的弧線（Arc patch）
            bond_edges = [(u, v) for (u, v) in G.edges()
                          if tm.get(tuple(sorted([u, v])), ('bond', 0))[0] == 'bond']
            for u, v in bond_edges:
                a1 = node_angle_deg[u]
                a2 = node_angle_deg[v]
                # 取較短弧；以 (a1+a2)/2 為中點，但要處理跨越 0 度的情況
                diff = (a2 - a1) % 360
                if diff > 180:
                    diff -= 360
                # 弧的兩端角度（matplotlib Arc 角度為「逆時針從 +x 軸」）
                if diff >= 0:
                    theta1, theta2 = a1, a1 + diff
                else:
                    theta1, theta2 = a1 + diff, a1
                arc = Arc((0, 0), 2.0, 2.0, angle=0,
                          theta1=theta1, theta2=theta2,
                          color=EDGE_COLORS['bond'], linewidth=1.2, alpha=0.85)
                ax.add_patch(arc)

            # 2) bridges、silk 用直線（穿越圓內部）
            for et in ['local_bridge', 'global_bridge', 'silk']:
                edge_list = [(u, v) for (u, v) in G.edges()
                             if tm.get(tuple(sorted([u, v])), ('bond', 0))[0] == et]
                if edge_list:
                    nx.draw_networkx_edges(G, pos, edgelist=edge_list,
                                           edge_color=EDGE_COLORS[et],
                                           width=1.0, alpha=0.9, ax=ax)
            # 3) 節點
            nx.draw_networkx_nodes(G, pos, node_size=4,
                                   node_color='black', ax=ax)
            ax.set_title(title, fontsize=11)
            ax.set_xlim(-1.15, 1.15); ax.set_ylim(-1.15, 1.15)
            ax.set_aspect('equal')
            ax.axis('off')
        else:
            ax.set_title(f"{title} (missing)"); ax.axis('off')
    fig.legend(handles=make_legend_handles(), loc='lower center',
               ncol=4, fontsize=11, bbox_to_anchor=(0.5, -0.02))
    plt.tight_layout()
    plt.savefig(out_dir / 'Fig6_visualization_smallworld.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print('  Fig. 6 saved')


def fig7(result_dir, out_dir):
    """Fig. 7: 16 個真實網路視覺化"""
    def bond_pct(fp):
        d = get_pct_dict(fp); return d['bond'] if d else 0

    ordered = sorted(NETWORK_MAP,
                     key=lambda nm: -bond_pct(result_dir / f'{nm[1]}_result.txt'))

    fig, axes = plt.subplots(4, 4, figsize=(16, 16))
    for idx, (paper, fname) in enumerate(ordered):
        ax = axes[idx // 4, idx % 4]
        fp = result_dir / f'{fname}_result.txt'
        if fp.exists():
            G, tm = build_graph(fp)
            draw_network(ax, G, tm, paper, layout='spring')
        else:
            ax.set_title(f"{paper} (missing)"); ax.axis('off')

    fig.legend(handles=make_legend_handles(), loc='lower center',
               ncol=4, fontsize=12, bbox_to_anchor=(0.5, -0.005))
    fig.suptitle('Fig. 7 — Link type identification for the 16 networks (sorted by bond ratio)',
                 fontsize=14, y=0.995)
    plt.tight_layout()
    plt.savefig(out_dir / 'Fig7_visualization_16networks.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print('  Fig. 7 saved')


def fig8(result_dir, out_dir):
    """Fig. 8: small-world 指紋掃描"""
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    plot_settings = [
        (axes[0], 'ws_swn',  '(a) WS small-world networks',
         'Random rewiring probability',  True),
        (axes[1], 'nws_swn', '(b) NW small-world networks',
         'Shortcut addition probability', False),
    ]
    for ax, prefix, title, xlabel, show_leg in plot_settings:
        pcts, labs = [], []
        for p in P_VALUES:
            fp = result_dir / f'{prefix}_{p}_result.txt'
            if fp.exists():
                d = get_pct_dict(fp)
                if d: pcts.append(d); labs.append(p)
        plot_stacked(ax, pcts, labs, title, xlabel, show_legend=show_leg)
    plt.tight_layout()
    plt.savefig(out_dir / 'Fig8_fingerprint_smallworld.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print('  Fig. 8 saved')


def fig9(result_dir, out_dir):
    """Fig. 9: 16 個真實網路指紋"""
    data = []
    for paper, fname in NETWORK_MAP:
        fp = result_dir / f'{fname}_result.txt'
        if fp.exists():
            d = get_pct_dict(fp)
            if d: data.append((paper, d))
    data.sort(key=lambda x: -x[1]['bond'])

    fig, ax = plt.subplots(figsize=(13, 6.3))
    plot_stacked(ax, [d[1] for d in data], [d[0] for d in data],
                 'Fig. 9 — Network fingerprints for the 16 networks (sorted by bond ratio)',
                 'Network')
    plt.tight_layout()
    plt.savefig(out_dir / 'Fig9_fingerprint_16networks.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print('  Fig. 9 saved')


def fig10(result_dir, out_dir):
    """Fig. 10: small-world 指紋相關矩陣
    (a) WS：色階 -0.8~0.8（紅藍對比，因為有正負相關）
    (b) NW：色階 0.75~1.0（單一深紅，因為全部高度相關）
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    plot_settings = [
        (axes[0], 'ws_swn',  '(a) WS small-world network',
         'Random rewiring probability',  'RdBu_r', -0.8, 0.8),
        (axes[1], 'nws_swn', '(b) NW small-world network',
         'Shortcut addition probability', 'Reds',   0.75, 1.0),
    ]
    for ax, prefix, title, xlabel, cmap, vmin, vmax in plot_settings:
        vecs, labs = [], []
        for p in P_VALUES:
            fp = result_dir / f'{prefix}_{p}_result.txt'
            if fp.exists():
                v = get_pct_vec(fp)
                if v is not None: vecs.append(v); labs.append(p)
        M = correlation_matrix(vecs)
        im = ax.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
        ax.set_xticks(range(len(labs))); ax.set_yticks(range(len(labs)))
        ax.set_xticklabels(labs, rotation=45, ha='right', fontsize=9)
        ax.set_yticklabels(labs, fontsize=9)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(xlabel, fontsize=11); ax.set_ylabel(xlabel, fontsize=11)
        plt.colorbar(im, ax=ax, shrink=0.8)
    plt.tight_layout()
    plt.savefig(out_dir / 'Fig10_corrmatrix_smallworld.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print('  Fig. 10 saved')


def fig11(result_dir, out_dir):
    """Fig. 11: 16 個真實網路指紋相關矩陣 + 階層聚類"""
    vecs, labs = [], []
    for paper, fname in NETWORK_MAP:
        fp = result_dir / f'{fname}_result.txt'
        if fp.exists():
            v = get_pct_vec(fp)
            if v is not None: vecs.append(v); labs.append(paper)

    M = correlation_matrix(vecs)
    dist = np.clip((1 - M + (1 - M).T) / 2, 0, 2)
    np.fill_diagonal(dist, 0)
    Z = linkage(squareform(dist, checks=False), method='average')

    fig = plt.figure(figsize=(11, 9))
    gs = GridSpec(2, 2, width_ratios=[1, 4], height_ratios=[1, 4],
                  hspace=0.02, wspace=0.02)
    ax_top = fig.add_subplot(gs[0, 1])
    dn_top = dendrogram(Z, no_labels=True, ax=ax_top, color_threshold=0,
                        above_threshold_color='black')
    ax_top.axis('off'); order = dn_top['leaves']

    ax_left = fig.add_subplot(gs[1, 0])
    dendrogram(Z, no_labels=True, ax=ax_left, orientation='left',
               color_threshold=0, above_threshold_color='black')
    ax_left.axis('off'); ax_left.invert_yaxis()

    ax_main = fig.add_subplot(gs[1, 1])
    M_reord = M[np.ix_(order, order)]
    labs_reord = [labs[i] for i in order]
    im = ax_main.imshow(M_reord, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    ax_main.set_xticks(range(len(labs_reord)))
    ax_main.set_yticks(range(len(labs_reord)))
    ax_main.set_xticklabels(labs_reord, rotation=90, fontsize=10)
    ax_main.set_yticklabels(labs_reord, fontsize=10)
    ax_main.yaxis.tick_right()
    cb = fig.colorbar(im, ax=ax_main, shrink=0.6, pad=0.15)
    cb.set_label('Correlation', fontsize=10)
    fig.suptitle('Fig. 11 — Correlation matrix for the 16 networks',
                 fontsize=13, y=0.98)
    plt.savefig(out_dir / 'Fig11_corrmatrix_16networks.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print('  Fig. 11 saved')


def fig12(result_dir, out_dir):
    """Fig. 12: 4 個示範網路階層式社群分割"""
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    networks = [
        ('rdgam',  '(a) rdgam'),
        ('women',  '(b) women'),
        ('camp92', '(c) camp92'),
        ('karate', '(d) karate'),
    ]
    for idx, (fname, title) in enumerate(networks):
        ax = axes[idx // 2, idx % 2]
        fp = result_dir / f'{fname}_result.txt'
        if fp.exists():
            G, tm = build_graph(fp)
            draw_partition(ax, G, tm, title)
        else:
            ax.set_title(f"{title} (missing)"); ax.axis('off')
    fig.legend(handles=make_legend_handles(), loc='lower center',
               ncol=4, fontsize=11, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle('Fig. 12 — Hierarchical community partition for 4 networks',
                 fontsize=13, y=0.995)
    plt.tight_layout()
    plt.savefig(out_dir / 'Fig12_partition.png', dpi=150, bbox_inches='tight')
    plt.close()
    print('  Fig. 12 saved')


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='生成論文 Fig. 4~12')
    parser.add_argument('--result-dir', type=str, default='.',
                        help='HETA 結果檔目錄（預設：當前目錄）')
    parser.add_argument('--out-dir', type=str, default='./figures',
                        help='輸出圖片目錄（預設：./figures）')
    args = parser.parse_args()

    result_dir = Path(args.result_dir).resolve()
    out_dir    = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f'Result dir : {result_dir}')
    print(f'Output dir : {out_dir}')
    print()
    print('Generating figures...')

    for fn in [fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12]:
        try:
            fn(result_dir, out_dir)
        except Exception as e:
            print(f'  [!] {fn.__name__} failed: {e}')

    print()
    print(f'Done. Figures saved to: {out_dir}')


if __name__ == '__main__':
    main()