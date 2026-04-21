import os, sys, networkx as nx

HETA_ROOT = os.path.abspath(os.path.dirname(__file__))
DATA = os.path.join(HETA_ROOT, 'data', 'net')
sys.path.insert(0, HETA_ROOT)
from heta import HETA as heta

# inline the table2 helpers so we don't trigger its os.chdir side-effect on import
import numpy as np

def first_layer_cn_ratio(G, u, v):
    nbr_u = set(G.neighbors(u)) - {v}
    nbr_v = set(G.neighbors(v)) - {u}
    if not nbr_u or not nbr_v:
        return 0.0
    return len(nbr_u & nbr_v) / min(len(nbr_u), len(nbr_v))

def flatten_communities(communities):
    out = []
    if communities is None: return out
    if isinstance(communities, (set, frozenset)):
        return [set(str(n) for n in communities)]
    if isinstance(communities, dict):
        for v in communities.values(): out.extend(flatten_communities(v))
        return out
    if isinstance(communities, (list, tuple)):
        for item in communities:
            if isinstance(item, (set, frozenset)): out.append(set(str(n) for n in item))
            elif isinstance(item, (list, tuple)): out.extend(flatten_communities(item))
            elif isinstance(item, dict): out.extend(flatten_communities(item))
            else: out.append({str(item)})
    return out

def compute_metrics(G, communities):
    communities = [c for c in communities if len(c) > 1]
    group_num = len(communities)
    in_group = set()
    for comm in communities: in_group.update(str(n) for n in comm)
    independent = [str(n) for n in G.nodes() if str(n) not in in_group]
    independent_num = len(independent)
    total_in_groups = sum(len(c) for c in communities)
    densities = []
    for comm in communities:
        nodes = [n for n in (str(x) for x in comm) if n in G.nodes()]
        sub = G.subgraph(nodes)
        if sub.number_of_nodes() > 1: densities.append(nx.density(sub))
    dA = float(np.mean(densities)) if densities else 0.0
    dS = float(np.std(densities))  if densities else 0.0
    return dict(group_num=group_num, independent_num=independent_num,
                total_in_groups=total_in_groups,
                density_avg=round(dA,4), density_std=round(dS,4))

for name, fn in [('karate','karate.net'), ('women','women.net'), ('camp92','USA.net')]:
    G = nx.Graph(nx.read_pajek(os.path.join(DATA, fn)))
    G.remove_edges_from(nx.selfloop_edges(G))
    G, _, _ = heta.bridge_or_bond(G, times=20)
    ret = heta.community_sorting(G)

    cur = flatten_communities(ret)
    print(f'=== {name}: N={G.number_of_nodes()} ===')
    print(f'  CURRENT  ({type(ret).__name__}→flatten): kept={[len(c) for c in cur if len(c)>1]}')
    flat_cdic = ret[2]
    correct = [set(str(n) for n in nodes) for nodes in flat_cdic.values()]
    print(f'  FIXED    (use flat_cdic): {len(correct)} groups, sizes={sorted((len(c) for c in correct), reverse=True)}')
    print(f'  metrics  :', compute_metrics(G, correct))
    print()
