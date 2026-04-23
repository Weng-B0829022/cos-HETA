# -*- coding: utf-8 -*-
"""
Skewed_exp.py — Experiment 5.1：R¹ 全網均值分佈常態性檢定（CLT 版）
=====================================================================

設計變更（2026-04-23 修訂）
--------------------------
舊版以「把每個 null 的所有邊 R¹ pool 成一個大樣本」方式做常態性檢定。
該做法在稀疏 k 值（如 k=6）下嚴重受零膨脹影響，edge-level 分佈永遠
偏離常態，無法有效驗證 Mean+2SD 閾值公式的抽樣論據。

本版改為「每個 null 產生 1 個 graph-level 均值標量」：

    μ_b = (1 / m) · Σ_{(u,v) ∈ E(G_b)} R^1(u, v; G_b)

對每個目標網路 G 生成 B 個 null G_1, ..., G_B，收集 B 個 μ 值，
檢驗這 B 個均值的抽樣分佈（sampling distribution of the sample mean）
是否服從常態——這才是中央極限定理（CLT）可直接保障的對象，也是
論文 Mean+2SD 閾值公式背後真正要用到的隨機變數。

兩種分母仍並行量測：

    (A) HPI / min 版本   （論文原版）
        R^1_{i,j} = |N(i)\\{j} ∩ N(j)\\{i}| / min(|N(i)\\{j}|, |N(j)\\{i}|)

    (B) Salton / sqrt 版（cos-HETA）
        R^1_{i,j} = |N(i)\\{j} ∩ N(j)\\{i}| / sqrt( |N(i)\\{j}| · |N(j)\\{i}| )

實驗步驟（依論文 5.1，CLT 版修訂）
----------------------------------
Step 1. 對每個目標網路 G 以 switching randomization 生成 B 個
        null network G_1, ..., G_B（保留度序列、打破結構相關）。
Step 2. 對每個 null G_b 計算全網邊 R¹ 均值 μ_b（min 與 sqrt 各一）。
        同一批 null 供兩版共用，以確保比較公平。
Step 3. 對 B 個均值的樣本計算：skewness、excess kurtosis、
        Shapiro-Wilk p-value、Anderson-Darling A²。
Step 4. 以 Mean + 2·SD 計算閾值（此即論文 T^k_E），並計算樣本中
        μ_b < 閾值的百分比，與理論常態近似 95.0% 相比。

視覺化
------
每個子圖：
  * histogram（density scale）
  * gaussian_kde 曲線疊加（紅色實線）
  * mean（點虛線）與 mean+2SD（紅虛線）垂直線
  * 右上角統計量文字框：skew / kurt / SW p / coverage%
  * 邊框色：|skew|<0.5 綠、0.5~1 橘、≥1 紅

輸出（每次執行產生 6 張圖 + 6 張 CSV + 1 個 run_info.txt）
---------------------------------------------------------
HETA/experiment/Skewed_exp/
  [min 版 — HPI，論文原版]
  - fig_ws_skew_min.png       WS 小世界 16 個 p 值（4×4）
  - fig_nw_skew_min.png       NW 小世界 16 個 p 值（4×4）
  - fig_real16_skew_min.png   16 個實證網路（4×4）
  - stats_ws_min.csv / stats_nw_min.csv / stats_real16_min.csv

  [sqrt 版 — Salton / cos-HETA]
  - fig_ws_skew_sqrt.png      WS 小世界 16 個 p 值（4×4）
  - fig_nw_skew_sqrt.png      NW 小世界 16 個 p 值（4×4）
  - fig_real16_skew_sqrt.png  16 個實證網路（4×4）
  - stats_ws_sqrt.csv / stats_nw_sqrt.csv / stats_real16_sqrt.csv

  [共同]
  - run_info.txt              執行參數紀錄 + 兩版偏態等級次數彙整

執行
----
    python3 Skewed_exp.py [--nulls B] [--Q Q] [--seed SEED]

預設 B=1000（論文 formal setting；pilot 可用 --nulls 200）。
Q 為每個 null network 進行 double edge swap 的倍數（Q × m 次）。
為使 min 與 sqrt 兩版使用同一批 null networks，兩版在同一次
Step 2 的 for-loop 內共用 Gr。
"""
from __future__ import annotations

import os
import sys
import csv
import math
import time
import argparse
from datetime import datetime

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from tqdm import tqdm


# ──────────────────────────────────────────────────────────────────────
# 路徑設定
# ──────────────────────────────────────────────────────────────────────
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_THIS_DIR, "data", "net")
OUT_DIR = os.path.join(_THIS_DIR, "experiment", "Skewed_exp")
os.makedirs(OUT_DIR, exist_ok=True)


# ──────────────────────────────────────────────────────────────────────
# 網路清單
# ──────────────────────────────────────────────────────────────────────
# WS / NW 的 16 個 p 值（對應 data/net/ws_swn_{p}.net 與 nws_swn_{p}.net）
P_VALUES = [
    "0.0",   "0.001", "0.002", "0.004",
    "0.008", "0.016", "0.032", "0.064",
    "0.128", "0.256", "0.384", "0.512",
    "0.640", "0.768", "0.896", "1.0",
]

# 16 個實證網路（對應論文 Fig. 5 的中位數由高至低排序）
REAL16_NETWORKS = [
    ("rdgam",       "14p.net"),
    ("women",       "women.net"),
    ("lesmis",      "lesmis.net"),
    ("jazz",        "jazz.net"),
    ("karate",      "karate.net"),
    ("camp92",      "USA.net"),
    ("ragusa16",    "Ragusa16.net"),
    ("football",    "football.net"),
    ("dolphins",    "dolphins.net"),
    ("celegans",    "celegans.net"),
    ("leader",      "leader.net"),
    ("prisonInter", "prisonInter.net"),
    ("k-core",      "k-core.net"),
    ("florentine",  "families.net"),
    ("ba_sfn",      "ba_sfn.net"),
    ("s208",        "s208.net"),
]


# ──────────────────────────────────────────────────────────────────────
# 核心函式
# ──────────────────────────────────────────────────────────────────────
def first_layer_cn_ratio_min(G: nx.Graph, u, v) -> float:
    """論文 Eq. (1)：HPI 版 R¹（分母 = min）。

    R^1_{i,j} = |N(i)\\{j} ∩ N(j)\\{i}| / min(|N(i)\\{j}|, |N(j)\\{i}|)

    若任一端扣除對端後無鄰居，回傳 0.0。
    """
    nbr_u = set(G.neighbors(u)) - {v}
    nbr_v = set(G.neighbors(v)) - {u}
    if not nbr_u or not nbr_v:
        return 0.0
    return len(nbr_u & nbr_v) / min(len(nbr_u), len(nbr_v))


def first_layer_cn_ratio_sqrt(G: nx.Graph, u, v) -> float:
    """Salton / cos-HETA 版 R¹（分母 = sqrt(|N(u)\\{v}| · |N(v)\\{u}|)）。

    R^1_{i,j} = |N(i)\\{j} ∩ N(j)\\{i}| / sqrt( |N(i)\\{j}| · |N(j)\\{i}| )

    幾何平均分母對兩端節點對稱，消除 HPI 在度異質網路上的結構偏誤。
    若任一端扣除對端後無鄰居，回傳 0.0。
    """
    nbr_u = set(G.neighbors(u)) - {v}
    nbr_v = set(G.neighbors(v)) - {u}
    if not nbr_u or not nbr_v:
        return 0.0
    return len(nbr_u & nbr_v) / math.sqrt(len(nbr_u) * len(nbr_v))


# metric key → (函式, 人類可讀的版本名稱, 公式字串)
METRIC_FUNCS = {
    "min":  (first_layer_cn_ratio_min,
             "HPI / min",
             r"$R^1 = |N(i)\backslash\{j\} \cap N(j)\backslash\{i\}| \; / \; \min(|N(i)\backslash\{j\}|, |N(j)\backslash\{i\}|)$"),
    "sqrt": (first_layer_cn_ratio_sqrt,
             "Salton / sqrt",
             r"$R^1 = |N(i)\backslash\{j\} \cap N(j)\backslash\{i\}| \; / \; \sqrt{|N(i)\backslash\{j\}| \cdot |N(j)\backslash\{i\}|}$"),
}


def load_graph(path: str) -> nx.Graph:
    """讀取 .net（Pajek）並移除自環，取最大連通成分副本。"""
    G = nx.Graph(nx.read_pajek(path))
    G.remove_edges_from(nx.selfloop_edges(G))
    if G.number_of_edges() == 0:
        return G
    largest = max(nx.connected_components(G), key=len)
    return G.subgraph(largest).copy()


def generate_r1_means(
    G: nx.Graph,
    num_nulls: int,
    Q: int,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """Step 1–2：switching randomization 產生 B 個 null，**每個 null 回傳 1 個均值**。

    這是 CLT 版的核心修訂：過去版本把所有 null 的所有邊 R¹ pool 成
    一個大樣本，本版改為每個 null 計算一個全網均值，B 個 null 共
    產出 B 個樣本點，以檢驗「樣本均值」的抽樣分佈（CLT 保證近常態）。

    為確保兩版的比較公平，兩種 R¹ 皆在同一批 null networks 上計算。

    Parameters
    ----------
    G : 目標網路（已取最大連通成分）
    num_nulls : null network 數量 B
    Q : 每個 null 執行 Q × m 次 double edge swap（swap 強度）
    rng : np.random.Generator，供可重現性

    Returns
    -------
    dict : {"min": ndarray_shape_(B,), "sqrt": ndarray_shape_(B,)}
           每個元素為對應 null 的 graph-level 均值
    """
    means: dict[str, list[float]] = {"min": [], "sqrt": []}
    m = G.number_of_edges()
    if m < 2:
        return {k: np.array(v, dtype=float) for k, v in means.items()}

    for _ in range(num_nulls):
        Gr = G.copy()
        swap_seed = int(rng.integers(1, 2**31 - 1))
        try:
            # connected_double_edge_swap 保證每次 swap 後圖仍連通
            nx.connected_double_edge_swap(Gr, nswap=Q * m, seed=swap_seed)
        except (nx.NetworkXAlgorithmError, nx.NetworkXError):
            # 對極小網路 swap 可能失敗，退而以原網路做樣本
            Gr = G.copy()
        # 在同一 null 上計算兩種版本的 R¹，分別取平均
        r1_min_vals = []
        r1_sqrt_vals = []
        for u, v in Gr.edges():
            r1_min_vals.append(first_layer_cn_ratio_min(Gr, u, v))
            r1_sqrt_vals.append(first_layer_cn_ratio_sqrt(Gr, u, v))
        if r1_min_vals:
            means["min"].append(float(np.mean(r1_min_vals)))
            means["sqrt"].append(float(np.mean(r1_sqrt_vals)))
        else:
            means["min"].append(float("nan"))
            means["sqrt"].append(float("nan"))
    return {k: np.asarray(v, dtype=float) for k, v in means.items()}


def compute_stats(vals: np.ndarray) -> dict:
    """Step 3–4：對 B 個均值樣本計算四組常態性統計量 + Mean+2SD 涵蓋率。"""
    vals = np.asarray(vals, dtype=float)
    vals = vals[~np.isnan(vals)]
    n = int(vals.size)
    out = {
        "n": n,
        "mean": float("nan"),
        "std": float("nan"),
        "skew": float("nan"),
        "excess_kurt": float("nan"),
        "sw_p": float("nan"),
        "ad_A2": float("nan"),
        "threshold": float("nan"),
        "cov_pct": float("nan"),
        "gap_from_95": float("nan"),
    }
    if n < 4:
        return out

    mean = float(np.mean(vals))
    std = float(np.std(vals, ddof=1))
    out["mean"] = mean
    out["std"] = std
    # 常數序列 std=0 時，skew/kurt/SW/AD 皆不可用；回傳 0 或保留 NaN
    if std > 0:
        out["skew"] = float(stats.skew(vals, bias=False))
        out["excess_kurt"] = float(stats.kurtosis(vals, fisher=True, bias=False))
    else:
        out["skew"] = 0.0
        out["excess_kurt"] = 0.0

    # Shapiro-Wilk（n=1000 遠小於 5000 上限）
    try:
        if std > 0:
            if n > 5000:
                sub = np.random.default_rng(42).choice(vals, size=5000, replace=False)
                sw = stats.shapiro(sub)
            else:
                sw = stats.shapiro(vals)
            out["sw_p"] = float(sw.pvalue)
    except Exception:
        pass

    # Anderson-Darling 對常態的 A² 統計量
    try:
        if std > 0:
            ad = stats.anderson(vals, dist="norm")
            out["ad_A2"] = float(ad.statistic)
    except Exception:
        pass

    # Mean + 2·SD 閾值與其涵蓋率
    threshold = mean + 2.0 * std
    out["threshold"] = threshold
    cov = float(np.mean(vals < threshold) * 100.0)
    out["cov_pct"] = cov
    out["gap_from_95"] = cov - 95.0
    return out


def skew_level(skew: float) -> tuple[str, str]:
    """回傳 (等級文字, 邊框顏色) — 用於圖上視覺標示偏態嚴重度。"""
    if np.isnan(skew):
        return "NA", "#888888"
    a = abs(skew)
    if a < 0.5:
        return "near-sym", "#2ca02c"     # 接近對稱，綠
    if a < 1.0:
        return "moderate", "#ff7f0e"     # 中度偏態，橘
    return "high", "#d62728"             # 高度偏態，紅


# ──────────────────────────────────────────────────────────────────────
# 繪圖
# ──────────────────────────────────────────────────────────────────────
def _annotate(ax, label: str, st: dict) -> None:
    """在單一子圖上標註名稱與統計量，並依偏度調整邊框顏色。"""
    level_txt, frame_color = skew_level(st["skew"])
    # 標題
    ax.set_title(label, fontsize=10, pad=3)
    # 邊框顏色顯示偏態等級
    for spine in ax.spines.values():
        spine.set_edgecolor(frame_color)
        spine.set_linewidth(1.8)

    # 右上角統計量文字框
    sw_p = st["sw_p"]
    sw_str = f"{sw_p:.2e}" if not np.isnan(sw_p) else "NA"
    txt = (
        f"skew={st['skew']:.2f}\n"
        f"kurt={st['excess_kurt']:.2f}\n"
        f"SW p={sw_str}\n"
        f"cov={st['cov_pct']:.1f}%"
    )
    ax.text(
        0.97, 0.96, txt,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=7.5, family="monospace",
        bbox=dict(boxstyle="round,pad=0.25",
                  facecolor="white", edgecolor=frame_color, alpha=0.85),
    )
    # 偏態等級 badge
    ax.text(
        0.03, 0.96, level_txt,
        transform=ax.transAxes, ha="left", va="top",
        fontsize=7.5, color=frame_color, weight="bold",
    )


def _plot_single_hist(ax, vals: np.ndarray, st: dict) -> None:
    """於子圖上繪 B 個 graph-level 均值的 histogram + KDE + 垂直參考線。"""
    vals = np.asarray(vals, dtype=float)
    vals = vals[~np.isnan(vals)]
    if vals.size == 0:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        return

    mu = st["mean"]
    sd = st["std"]
    # 自適應 x 範圍：以 mean ± 4·SD 為視窗（常態下涵蓋 >99.99%）
    if np.isfinite(mu) and np.isfinite(sd) and sd > 0:
        lo = max(0.0, mu - 4 * sd)
        hi = mu + 4 * sd
    else:
        vmin = float(np.min(vals))
        vmax = float(np.max(vals))
        if vmax - vmin < 1e-12:
            lo = max(0.0, vmin - 1e-3)
            hi = vmax + 1e-3
        else:
            pad = (vmax - vmin) * 0.1
            lo = max(0.0, vmin - pad)
            hi = vmax + pad

    # histogram
    ax.hist(vals, bins=30, range=(lo, hi), density=True,
            color="#4C78A8", alpha=0.70, edgecolor="white", linewidth=0.3)

    # KDE 疊加
    try:
        if vals.size >= 2 and sd > 0:
            kde = stats.gaussian_kde(vals)
            xs = np.linspace(lo, hi, 200)
            ax.plot(xs, kde(xs), color="#e45756", linewidth=1.3)
    except Exception:
        pass

    if np.isfinite(mu):
        ax.axvline(mu, color="#333333", linestyle=":", linewidth=1.0)
    thr = st["threshold"]
    if np.isfinite(thr) and lo <= thr <= hi:
        ax.axvline(thr, color="#d62728", linestyle="--", linewidth=1.1)

    ax.set_xlim(lo, hi)
    ax.tick_params(axis="both", which="major", labelsize=7)
    # 關閉科學記號的偏移標註，讓極小均值的 x 軸更好讀
    try:
        ax.xaxis.get_major_formatter().set_useOffset(False)
    except Exception:
        pass


def render_grid_figure(
    title_suffix: str,
    labels: list[str],
    pools: list[np.ndarray],
    stats_list: list[dict],
    out_path: str,
    fig_title: str,
) -> None:
    """4×4 子圖版面。"""
    assert len(labels) == 16, "此函式預設 16 個子圖（4×4）"
    fig, axes = plt.subplots(4, 4, figsize=(15, 13), constrained_layout=False)
    fig.suptitle(fig_title, fontsize=14, y=0.995)

    for ax, lab, vals, st in zip(axes.flat, labels, pools, stats_list):
        _plot_single_hist(ax, vals, st)
        _annotate(ax, lab, st)

    # 軸標
    for r in range(4):
        axes[r, 0].set_ylabel("density", fontsize=9)
    for c in range(4):
        axes[3, c].set_xlabel(r"mean $R^1_{i,j}$ per null", fontsize=9)

    # 偏態等級圖例
    handles = [
        plt.Rectangle((0, 0), 1, 1, fc="white", ec="#2ca02c", lw=1.8, label="|skew|<0.5 near-sym"),
        plt.Rectangle((0, 0), 1, 1, fc="white", ec="#ff7f0e", lw=1.8, label="0.5≤|skew|<1 moderate"),
        plt.Rectangle((0, 0), 1, 1, fc="white", ec="#d62728", lw=1.8, label="|skew|≥1 high"),
        plt.Line2D([0], [0], color="#e45756", lw=1.3, label="KDE"),
        plt.Line2D([0], [0], color="#333333", ls=":", label="mean"),
        plt.Line2D([0], [0], color="#d62728", ls="--", label="mean+2SD"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 0.005))

    plt.tight_layout(rect=(0.0, 0.035, 1.0, 0.975))
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  ➜ 已輸出 {out_path}")


# ──────────────────────────────────────────────────────────────────────
# 統計表寫出
# ──────────────────────────────────────────────────────────────────────
def write_stats_csv(csv_path: str, labels: list[str], stats_list: list[dict]) -> None:
    fields = ["label", "n", "mean", "std", "skew", "excess_kurt",
              "sw_p", "ad_A2", "threshold", "cov_pct", "gap_from_95"]
    with open(csv_path, "w", newline="", encoding="utf-8") as fp:
        w = csv.writer(fp)
        w.writerow(fields)
        for lab, st in zip(labels, stats_list):
            w.writerow([lab] + [f"{st[k]:.6g}" if isinstance(st[k], float) else st[k]
                                for k in fields[1:]])
    print(f"  ➜ 已輸出 {csv_path}")


# ──────────────────────────────────────────────────────────────────────
# 主流程
# ──────────────────────────────────────────────────────────────────────
def process_group(
    group_name: str,
    label_file_pairs: list[tuple[str, str]],
    num_nulls: int,
    Q: int,
    rng: np.random.Generator,
    desc: str,
) -> tuple[list[str], dict[str, list[np.ndarray]], dict[str, list[dict]]]:
    """對整組網路生成 null 並同時計算 min / sqrt 兩版 graph-level 均值與統計量。

    Returns
    -------
    labels : list[str]
    pools_per_metric : {"min": [ndarray(B,), ...], "sqrt": [ndarray(B,), ...]}
    stats_per_metric : {"min": [dict, ...],        "sqrt": [dict, ...]}
    """
    labels: list[str] = []
    pools_per_metric: dict[str, list[np.ndarray]] = {"min": [], "sqrt": []}
    stats_per_metric: dict[str, list[dict]] = {"min": [], "sqrt": []}

    for label, fname in tqdm(label_file_pairs, desc=desc):
        path = os.path.join(DATA_DIR, fname)
        labels.append(label)
        if not os.path.isfile(path):
            print(f"  ! 檔案不存在：{path}")
            for mk in ("min", "sqrt"):
                pools_per_metric[mk].append(np.array([]))
                stats_per_metric[mk].append(compute_stats(np.array([])))
            continue
        G = load_graph(path)
        means = generate_r1_means(G, num_nulls=num_nulls, Q=Q, rng=rng)
        for mk in ("min", "sqrt"):
            pools_per_metric[mk].append(means[mk])
            stats_per_metric[mk].append(compute_stats(means[mk]))
    return labels, pools_per_metric, stats_per_metric


def main():
    ap = argparse.ArgumentParser(
        description="Experiment 5.1 — R¹ graph-mean distribution normality (CLT version)")
    ap.add_argument("--nulls", type=int, default=1000,
                    help="每個網路生成的 null network 數量 B（預設 1000 = 論文 formal）")
    ap.add_argument("--Q", type=int, default=100,
                    help="switching randomization 強度：每個 null 執行 Q × m 次 swap（預設 100 = 論文標準）")
    ap.add_argument("--seed", type=int, default=20260423,
                    help="隨機種子（預設 20260423）")
    args = ap.parse_args()

    t_start = time.time()

    print(f"=== Skewed_exp (Experiment 5.1 — graph-mean CLT version) ===")
    print(f"NUM_NULLS = {args.nulls}")
    print(f"Q         = {args.Q}")
    print(f"SEED      = {args.seed}")
    print(f"DATA_DIR  = {DATA_DIR}")
    print(f"OUT_DIR   = {OUT_DIR}")
    print()

    rng = np.random.default_rng(args.seed)

    # ---------- Group 1: WS ----------
    ws_pairs = [(p, f"ws_swn_{p}.net") for p in P_VALUES]
    ws_labels, ws_pools_dual, ws_stats_dual = process_group(
        "WS", ws_pairs, args.nulls, args.Q, rng, desc="WS 小世界")

    # ---------- Group 2: NW ----------
    nw_pairs = [(p, f"nws_swn_{p}.net") for p in P_VALUES]
    nw_labels, nw_pools_dual, nw_stats_dual = process_group(
        "NW", nw_pairs, args.nulls, args.Q, rng, desc="NW 小世界")

    # ---------- Group 3: 16 Real ----------
    real_labels, real_pools_dual, real_stats_dual = process_group(
        "REAL16", REAL16_NETWORKS, args.nulls, args.Q, rng, desc="16 個實證網路")

    # ---------- 渲染 6 張圖 + 6 張 CSV ----------
    groups = [
        ("ws",     "WS  p=",      ws_labels,    ws_pools_dual,    ws_stats_dual,
         "Figure 5.1 (a) -- WS small-world networks", "(16 p-values)"),
        ("nw",     "NW  p=",      nw_labels,    nw_pools_dual,    nw_stats_dual,
         "Figure 5.1 (b) -- NW small-world networks", "(16 p-values)"),
        ("real16", "",            real_labels,  real_pools_dual,  real_stats_dual,
         "Figure 5.1 (c) -- 16 Real networks",         ""),
    ]

    for g_key, prefix, labels, pools_dual, stats_dual, title_main, title_tail in groups:
        for mk in ("min", "sqrt"):
            display_labels = [f"{prefix}{lab}" if prefix else lab for lab in labels]
            metric_label = METRIC_FUNCS[mk][1]
            fig_title = (f"{title_main}: distribution of graph-mean "
                         f"R^1 over nulls [{metric_label}] {title_tail}").strip()
            fig_title = fig_title.replace("R^1", r"$\overline{R^1_{i,j}}$")
            render_grid_figure(
                title_suffix=f"{g_key}_{mk}",
                labels=display_labels,
                pools=pools_dual[mk],
                stats_list=stats_dual[mk],
                out_path=os.path.join(OUT_DIR, f"fig_{g_key}_skew_{mk}.png"),
                fig_title=fig_title,
            )
            write_stats_csv(
                os.path.join(OUT_DIR, f"stats_{g_key}_{mk}.csv"),
                labels, stats_dual[mk],
            )

    # ---------- run_info.txt ----------
    info_path = os.path.join(OUT_DIR, "run_info.txt")
    t_end = time.time()
    with open(info_path, "w", encoding="utf-8") as fp:
        fp.write("Experiment 5.1 — R¹ Graph-Mean Distribution Normality (CLT version)\n")
        fp.write("=" * 72 + "\n")
        fp.write(f"timestamp    : {datetime.now().isoformat(timespec='seconds')}\n")
        fp.write(f"num_nulls B  : {args.nulls}\n")
        fp.write(f"Q (swap × m) : {args.Q}\n")
        fp.write(f"seed         : {args.seed}\n")
        fp.write(f"runtime(min) : {(t_end - t_start)/60:.2f}\n")
        fp.write(f"data_dir     : {DATA_DIR}\n")
        fp.write(f"out_dir      : {OUT_DIR}\n")
        fp.write("\n")
        fp.write("Sample unit  : per-null graph-mean μ_b = (1/m) · Σ_(u,v) R^1(u,v; G_b)\n")
        fp.write("Sample size  : B means per network (one scalar per null)\n")
        fp.write("\n")
        fp.write("Similarity indices:\n")
        fp.write("  [min]  HPI        : R^1 = |∩| / min(|N(u)\\{v}|, |N(v)\\{u}|)   (論文原版)\n")
        fp.write("  [sqrt] Salton/cos : R^1 = |∩| / sqrt(|N(u)\\{v}| · |N(v)\\{u}|) (cos-HETA)\n")
        fp.write("Null model : NetworkX connected_double_edge_swap（保留度序列）\n")
        fp.write("Stats      : skewness, excess kurtosis, Shapiro-Wilk p, Anderson-Darling A², Mean+2SD coverage\n")
        fp.write("\n")
        fp.write("兩版共用同一批 null networks，以確保比較公平（同一 rng seed 序列）。\n")
        fp.write("\n")

        fp.write("── Summary 偏態等級次數（|skew|<0.5: near-sym; 0.5~1: moderate; ≥1: high） ─\n")
        fp.write(f"{'group':<8} {'metric':<6} {'near-sym':>9} {'moderate':>9} {'high':>5} {'NA':>4}\n")
        for name, sts_dual in [
            ("WS",     ws_stats_dual),
            ("NW",     nw_stats_dual),
            ("REAL16", real_stats_dual),
        ]:
            for mk in ("min", "sqrt"):
                counts = {"near-sym": 0, "moderate": 0, "high": 0, "NA": 0}
                for st in sts_dual[mk]:
                    counts[skew_level(st["skew"])[0]] += 1
                fp.write(
                    f"{name:<8} {mk:<6} "
                    f"{counts['near-sym']:>9} {counts['moderate']:>9} "
                    f"{counts['high']:>5} {counts['NA']:>4}\n"
                )
        fp.write("\n")

        # 偏度 mean / 絕對值 mean 摘要（方便快速對比兩版）
        fp.write("── 平均偏度對照（越接近 0 越接近常態） ────────────────\n")
        fp.write(f"{'group':<8} {'mean(skew)_min':>16} {'mean(skew)_sqrt':>17} {'mean|skew|_min':>16} {'mean|skew|_sqrt':>17}\n")
        for name, sts_dual in [
            ("WS",     ws_stats_dual),
            ("NW",     nw_stats_dual),
            ("REAL16", real_stats_dual),
        ]:
            ms = [st["skew"] for st in sts_dual["min"]]
            ss = [st["skew"] for st in sts_dual["sqrt"]]
            ms_a = [abs(x) for x in ms if not np.isnan(x)]
            ss_a = [abs(x) for x in ss if not np.isnan(x)]
            ms_m = np.nanmean(ms) if len(ms) else float("nan")
            ss_m = np.nanmean(ss) if len(ss) else float("nan")
            ms_abs = np.mean(ms_a) if ms_a else float("nan")
            ss_abs = np.mean(ss_a) if ss_a else float("nan")
            fp.write(f"{name:<8} {ms_m:>16.3f} {ss_m:>17.3f} {ms_abs:>16.3f} {ss_abs:>17.3f}\n")

    print(f"  ➜ 已輸出 {info_path}")
    print(f"\n總耗時：{(t_end - t_start)/60:.2f} 分鐘")


if __name__ == "__main__":
    main()
