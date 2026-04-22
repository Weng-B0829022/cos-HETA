#from heta.HETA import *

# -*- coding: utf-8 -*-
"""cos_heta 套件入口。

cos-HETA 版本：把論文 Eq. 1 的 R¹ 分母由 min 改為幾何平均，
與 cos_heta/HETA.py、cos_heta/HETA2.py 內的 HETA 距離公式
分母定義一致（sqrt(|N_s[l]| · |N_t[l]|)）。

此函式會在 cos_auto.py 執行時，經由 _cos_shim/heta/__init__.py
以 `heta.first_layer_cn_ratio` 的名義轉接出去，讓 src/ 腳本
統一透過 `from heta import first_layer_cn_ratio` 取用。
"""

import math


def first_layer_cn_ratio(G, u, v):
    """cos-HETA 版 Eq. 1：第一層共同鄰居比率 R¹（分母＝幾何平均）。

    Parameters
    ----------
    G : networkx.Graph
    u, v : hashable
        邊的兩端點。

    Returns
    -------
    float
        len(N(u)\\{v} ∩ N(v)\\{u}) / sqrt(|N(u)\\{v}| · |N(v)\\{u}|)
        若任一端沒有其他鄰居則回傳 0.0。
    """
    nbr_u = set(G.neighbors(u)) - {v}
    nbr_v = set(G.neighbors(v)) - {u}
    if not nbr_u or not nbr_v:
        return 0.0
    return len(nbr_u & nbr_v) / math.sqrt(len(nbr_u) * len(nbr_v))
