# -*- coding: utf-8 -*-
"""
====================================================================
HETA 論文（Huang, Chin, Fu, Tsai, 2019, Physica A 536, 121027）
   §3 Experimental Results 全部 10 項實驗 + Table 1 自動化主控腳本
====================================================================

用途
----
一鍵執行（或選擇性執行）論文 §3 的全部圖表復刻：
    Fig. 4, Fig. 5, Fig. 6, Fig. 7, Fig. 8, Fig. 9,
    Fig. 10, Fig. 11, Fig. 12, Table 2
另外列出 Table 1 的 16 個實證網路基本指標（非實驗，純資料規格）。

執行前置
--------
1. 已安裝 heta 套件：
       cd HETA
       pip install -e .
2. 已安裝相依套件：
       pip install networkx>=3.0 numpy scipy matplotlib tqdm pathos seaborn pandas

使用方式
--------
直接全部跑（預設）：
    python3 auto.py

只列出實驗清單與當前開關狀態：
    python3 auto.py --list

只跑指定的實驗：
    python3 auto.py --only fig4,fig5,table1

跳過某幾個實驗：
    python3 auto.py --skip fig6,fig7,fig8,table2

降低隨機化次數做快速驗收（預設 1000，建議測試時用 50–100）：
    python3 auto.py --times 100

輸出位置
--------
每個實驗的 .py 與輸出（PNG / CSV）皆放在自己的子資料夾，
所有 exp* 與 table* 資料夾統一收納於 HETA/experiment/ 之下：
    HETA/experiment/exp4_smallworld_r1_boxplot/fig4.png
    HETA/experiment/exp5_real16_r1_boxplot/fig5.png
    ...
    HETA/experiment/table2_partition_metrics/table2.csv
共享快取（exp10 重用 exp8、exp11 重用 exp9）仍位於 HETA 根目錄：
    HETA/_cache/

手動開關
--------
直接修改下方 `TOGGLE` 字典即可永久關閉某實驗：
    TOGGLE["fig6"] = False
"""

import argparse
import os
import subprocess
import sys
import time
from datetime import datetime

ROOT = os.path.dirname(os.path.abspath(__file__))
# 所有 exp* / table* 資料夾統一收納於 HETA/experiment/ 之下
EXP_DIR = "experiment"

# ────────────────────────────────────────────────────────────────────
#  實驗清單（key, 子資料夾名, 腳本檔名, 中文說明, 預設啟用, 預估耗時）
# ────────────────────────────────────────────────────────────────────
EXPERIMENTS = [
    ("table1", "table1_dataset_info",
     "table1_dataset_info.py",
     "Table 1：16 個實證網路基本指標（非實驗，列出數值並存 CSV）",
     True,  "< 1 分"),

    ("fig4",   "exp4_smallworld_r1_boxplot",
     "exp4_smallworld_r1_boxplot.py",
     "Fig. 4：WS/NW 32 個小世界網路 R¹ 箱型圖",
     True,  "< 1 分"),

    ("fig5",   "exp5_real16_r1_boxplot",
     "exp5_real16_r1_boxplot.py",
     "Fig. 5：16 個實證網路 R¹ 箱型圖",
     True,  "< 1 分"),

    ("fig6",   "exp6_smallworld_link_identification",
     "exp6_smallworld_link_identification.py",
     "Fig. 6：WS/NW 在 p=0.004/0.064/0.256 的連結識別（圓形佈局）",
     True,  "10–30 分"),

    ("fig7",   "exp7_real16_link_identification",
     "exp7_real16_link_identification.py",
     "Fig. 7：16 個實證網路連結識別（force-directed 佈局）",
     True,  "30–90 分"),

    ("fig8",   "exp8_smallworld_fingerprint",
     "exp8_smallworld_fingerprint.py",
     "Fig. 8：WS/NW 32 個小世界網路指紋堆疊長條",
     True,  "60–120 分（產生 _cache/fp_ws.pkl 與 fp_nw.pkl）"),

    ("fig9",   "exp9_real16_fingerprint",
     "exp9_real16_fingerprint.py",
     "Fig. 9：16 個實證網路指紋堆疊長條",
     True,  "30–60 分（產生 _cache/fp_real16.pkl）"),

    ("fig10",  "exp10_smallworld_correlation",
     "exp10_smallworld_correlation.py",
     "Fig. 10：WS/NW 指紋相關係數矩陣（依賴 _cache/fp_ws.pkl, fp_nw.pkl）",
     True,  "< 1 分（若 fig8 已跑過）；否則 60–120 分"),

    ("fig11",  "exp11_real16_correlation_clustermap",
     "exp11_real16_correlation_clustermap.py",
     "Fig. 11：16 實證網路指紋相關矩陣 + 階層分群（依賴 _cache/fp_real16.pkl）",
     True,  "< 1 分（若 fig9 已跑過）；否則 30–60 分"),

    ("fig12",  "exp12_community_partition",
     "exp12_community_partition.py",
     "Fig. 12：4 個代表網路（rdgam/women/camp92/karate）階層社群分割",
     True,  "1–5 分"),

    ("table2", "table2_partition_metrics",
     "table2_partition_metrics.py",
     "Table 2：16 網路社群分割 6 項量化指標 → table2.csv",
     True,  "30–90 分"),
]

# ────────────────────────────────────────────────────────────────────
#  開關區（手動編輯：True = 跑，False = 略過）
# ────────────────────────────────────────────────────────────────────
TOGGLE = {key: enabled for key, _, _, _, enabled, _ in EXPERIMENTS}

# 範例：永久關閉特別耗時的兩個（取消註解即生效）
# TOGGLE["fig7"] = False
# TOGGLE["fig8"] = False


# ────────────────────────────────────────────────────────────────────
#  內部函式
# ────────────────────────────────────────────────────────────────────
def _print_list():
    """列出所有實驗的開關狀態與預估耗時"""
    print("\n" + "=" * 90)
    print(f"  {'KEY':<8} {'STATUS':<8} {'TIME':<40} 說明")
    print("=" * 90)
    for key, folder, script, desc, _, eta in EXPERIMENTS:
        mark = "[ON]" if TOGGLE.get(key, False) else "[OFF]"
        print(f"  {key:<8} {mark:<8} {eta:<40} {desc}")
    print("=" * 90 + "\n")


def _apply_filters(only, skip):
    """根據 CLI 的 --only / --skip 修改 TOGGLE"""
    if only:
        wanted = {x.strip() for x in only.split(",") if x.strip()}
        for key in TOGGLE:
            TOGGLE[key] = (key in wanted)
        unknown = wanted - set(TOGGLE.keys())
        if unknown:
            print(f"⚠  --only 指定了不存在的實驗 key：{unknown}")
    if skip:
        for s in (x.strip() for x in skip.split(",") if x.strip()):
            if s in TOGGLE:
                TOGGLE[s] = False
            else:
                print(f"⚠  --skip 指定了不存在的實驗 key：{s}")


def _run_one(key, folder, script, desc, times):
    """以 subprocess 執行單一實驗，於該實驗的子資料夾內運作"""
    # 所有實驗資料夾已搬至 HETA/experiment/ 之下
    folder_path = os.path.join(ROOT, EXP_DIR, folder)
    script_path = os.path.join(folder_path, script)
    if not os.path.exists(script_path):
        print(f"❌ 找不到腳本：{script_path}")
        return False, 0.0

    print("\n" + "─" * 90)
    print(f"▶  [{key}] {desc}")
    print(f"   工作目錄：{folder_path}")
    print(f"   執行腳本：{script}")
    print(f"   隨機化次數 HETA_TIMES = {times}")
    print("─" * 90)

    env = os.environ.copy()
    env["HETA_TIMES"] = str(times)
    # 確保子腳本能 import heta（讓 HETA/ 在 sys.path 上）
    py_path = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = ROOT + (os.pathsep + py_path if py_path else "")

    t0 = time.time()
    try:
        result = subprocess.run(
            [sys.executable, script_path],
            cwd=folder_path,
            env=env,
            check=False,
        )
        ok = (result.returncode == 0)
    except Exception as e:
        print(f"❌ 例外：{e}")
        ok = False
    elapsed = time.time() - t0

    status = "✅ 成功" if ok else "❌ 失敗"
    print(f"\n   {status}（耗時 {elapsed:.1f} 秒）")
    return ok, elapsed


def _write_run_log(start_dt, end_dt, grand_elapsed, args, summary):
    """
    將本次執行的參數與時間寫入 HETA/experiment/run_<timestamp>.txt
    每次執行都會產生獨立的 log 檔，不覆蓋歷史。

    參數：
        start_dt       datetime 執行開始時間
        end_dt         datetime 執行結束時間
        grand_elapsed  float    總耗時（秒）
        args           argparse.Namespace  CLI 參數
        summary        list[(key, ok, elapsed)]  各實驗結果
    """
    exp_root = os.path.join(ROOT, EXP_DIR)
    os.makedirs(exp_root, exist_ok=True)

    ts = start_dt.strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(exp_root, f"run_{ts}.txt")

    lines = []
    lines.append("=" * 72)
    lines.append("  HETA 論文實驗自動化 — 本次執行參數與結果")
    lines.append("=" * 72)
    lines.append("")
    lines.append(f"執行開始：{start_dt.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"執行結束：{end_dt.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"總耗時  ：{grand_elapsed:.1f} 秒（{grand_elapsed/60:.2f} 分）")
    lines.append("")

    lines.append("── 執行參數 " + "─" * 60)
    lines.append(f"HETA_TIMES（每網路隨機化次數）：{args.times}")
    lines.append(f"--only：{args.only if args.only else '（未指定，預設全跑）'}")
    lines.append(f"--skip：{args.skip if args.skip else '（未指定）'}")
    lines.append(f"Python ：{sys.executable}")
    lines.append(f"HETA ROOT：{ROOT}")
    lines.append(f"EXP_DIR  ：{os.path.join(ROOT, EXP_DIR)}")
    lines.append("")

    lines.append("── 實驗開關狀態 " + "─" * 56)
    lines.append(f"  {'KEY':<8} {'STATUS':<8} {'TIME':<40} 說明")
    lines.append("  " + "-" * 68)
    for key, folder, script, desc, _, eta in EXPERIMENTS:
        mark = "[ON]" if TOGGLE.get(key, False) else "[OFF]"
        lines.append(f"  {key:<8} {mark:<8} {eta:<40} {desc}")
    lines.append("")

    lines.append("── 執行結果 " + "─" * 60)
    lines.append(f"  {'KEY':<8} {'RESULT':<10} {'TIME (sec)':>12}")
    lines.append("  " + "-" * 40)
    for key, ok, elapsed in summary:
        lines.append(
            f"  {key:<8} {('✅ OK' if ok else '❌ FAIL'):<10} {elapsed:>12.1f}"
        )
    lines.append("")

    n_ok = sum(1 for _, ok, _ in summary if ok)
    n_total = len(summary)
    lines.append("── 最終統計 " + "─" * 60)
    lines.append(f"  通過：{n_ok}/{n_total}")
    lines.append(f"  失敗：{n_total - n_ok}/{n_total}")
    lines.append("=" * 72)
    lines.append("")

    with open(log_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return log_path


# ────────────────────────────────────────────────────────────────────
#  主流程
# ────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="HETA 論文實驗一鍵自動化（10 實驗 + Table 1）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--list", action="store_true",
                        help="只列出實驗與當前開關狀態，不執行")
    parser.add_argument("--only", type=str, default="",
                        help="只跑指定實驗，逗號分隔，例如：--only fig4,fig5")
    parser.add_argument("--skip", type=str, default="",
                        help="跳過指定實驗，逗號分隔，例如：--skip fig6,fig7")
    parser.add_argument("--times", type=int, default=1000,
                        help="randomization 次數（預設 1000；快速驗收建議 50–100）")
    args = parser.parse_args()

    _apply_filters(args.only, args.skip)

    if args.list:
        _print_list()
        return

    print("\n" + "█" * 90)
    print("█" + " HETA 論文實驗自動化主控（auto.py） ".center(88) + "█")
    print("█" * 90)
    _print_list()

    enabled = [(k, f, s, d) for (k, f, s, d, _, _) in EXPERIMENTS if TOGGLE.get(k, False)]
    if not enabled:
        print("⚠  沒有任何實驗被啟用，結束。")
        return

    print(f"準備執行 {len(enabled)} 個實驗，HETA_TIMES = {args.times}\n")

    summary = []
    start_dt = datetime.now()
    grand_t0 = time.time()
    for key, folder, script, desc in enabled:
        ok, elapsed = _run_one(key, folder, script, desc, args.times)
        summary.append((key, ok, elapsed))
    grand_elapsed = time.time() - grand_t0
    end_dt = datetime.now()

    # 最終彙總
    print("\n" + "═" * 90)
    print("  最終彙總")
    print("═" * 90)
    print(f"  {'KEY':<8} {'RESULT':<10} {'TIME (sec)':>12}")
    print("  " + "-" * 80)
    for key, ok, elapsed in summary:
        print(f"  {key:<8} {('✅ OK' if ok else '❌ FAIL'):<10} {elapsed:>12.1f}")
    print("  " + "-" * 80)
    n_ok = sum(1 for _, ok, _ in summary if ok)
    print(f"  通過：{n_ok}/{len(summary)}    總耗時：{grand_elapsed:.1f} 秒"
          f"（{grand_elapsed/60:.1f} 分）")
    print("═" * 90 + "\n")

    # 將本次執行的參數與時間寫入 HETA/experiment/run_<timestamp>.txt
    try:
        log_path = _write_run_log(start_dt, end_dt, grand_elapsed, args, summary)
        print(f"本次執行紀錄已寫入：{log_path}\n")
    except Exception as e:
        print(f"⚠  寫入執行紀錄失敗：{e}\n")


if __name__ == "__main__":
    main()
