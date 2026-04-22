#from heta.HETA import *

# -*- coding: utf-8 -*-
"""heta 套件入口。

此處除了保留原有的 HETA 主模組匯出外，另外提供一個論文 Eq. 1 的
第一層共同鄰居比例 R¹ 的共用函式，供 src/ 下各實驗腳本統一匯入。

原版（heta/）分母 = min(|N_u \\ {v}|, |N_v \\ {u}|)
幾何平均版（cos_heta/）分母 = sqrt(|N_u \\ {v}| · |N_v \\ {u}|)

兩套套件必須暴露同名函式，才能讓 cos_auto.py 的 shim 以相同介面
將呼叫從 heta 轉接到 cos_heta，讓 src/ 腳本透過
    from heta import first_layer_cn_ratio
就能自動跟著主控腳本（auto.py 或 cos_auto.py）切換分母定義。
"""


def first_layer_cn_ratio(G, u, v):
    """論文 Eq. 1：第一層共同鄰居比率 R¹（分母＝min；原版 HETA）。

    Parameters
    ----------
    G : networkx.Graph
    u, v : hashable
        邊的兩端點。

    Returns
    -------
    float
        len(N(u)\\{v} ∩ N(v)\\{u}) / min(|N(u)\\{v}|, |N(v)\\{u}|)
        若任一端沒有其他鄰居則回傳 0.0。
    """
    nbr_u = set(G.neighbors(u)) - {v}
    nbr_v = set(G.neighbors(v)) - {u}
    if not nbr_u or not nbr_v:
        return 0.0
    return len(nbr_u & nbr_v) / min(len(nbr_u), len(nbr_v))
