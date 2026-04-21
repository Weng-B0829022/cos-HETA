# -*- coding: utf-8 -*-
"""
====================================================================
cos-HETA 自動化主控腳本（cos_auto.py）
   §3 Experimental Results 全部 10 項實驗 + Table 1
   — HETA 幾何平均分母變體（cos_heta 套件）
====================================================================

本腳本為 auto.py 的姊妹版本。行為、CLI 旗標、實驗清單、輸出資料夾結
構、run_info.txt 格式、共享快取位置等全部比照原本；唯一差異為：
改用 `cos_heta/` 內的新演算法（分母由 min 改為幾何平均
sqrt((d_u-1)(d_v-1))）取代原 `heta/` 套件。

為了不修改 src/ 下各實驗腳本（它們統一寫
`from heta import HETA as heta`），本腳本會在執行時自動於
`_cos_shim/heta/` 產生一個轉接（shim）套件，把 `heta` 模組指向
`cos_heta/` 的實作，並把該 shim 路徑插到子程序 PYTHONPATH 最前
方，使 `import heta` 自動改走 cos_heta 版本。

用途
----
一鍵執行（或選擇性執行）論文 §3 的全部圖表復刻（以 cos-HETA 演算法）：
    Fig. 4, Fig. 5, Fig. 6, Fig. 7, Fig. 8, Fig. 9,
    Fig. 10, Fig. 11, Fig. 12, Table 2
另外列出 Table 1 的 16 個實證網路基本指標（非實驗，純資料規格）。

執行前置
--------
1. 已安裝 heta 套件（原版，部分工具與相依仍會查找 heta）：
       cd HETA
       pip install -e .
2. 已安裝相依套件：
       pip install networkx>=3.0 numpy scipy matplotlib tqdm pathos seaborn pandas

使用方式
--------
直接全部跑（預設）：
    python3 cos_auto.py

只列出實驗清單與當前開關狀態：
    python3 cos_auto.py --list

只跑指定的實驗：
    python3 cos_auto.py --only fig4,fig5,table1

跳過某幾個實驗：
    python3 cos_auto.py --skip fig6,fig7,fig8,table2

降低隨機化次數做快速驗收（預設 1000，建議測試時用 50–100）：
    python3 cos_auto.py --times 100

專案結構
--------
    HETA/
    ├── auto.py                                ← 原始 HETA 主控腳本
    ├── cos_auto.py                            ← 本主控腳本（cos-HETA 版）
    ├── heta/                                  ← 原始套件（pip install -e .）
    ├── cos_heta/                              ← 幾何平均分母版本（本腳本使用）
    ├── _cos_shim/heta/__init__.py             ← 自動產生之 heta→cos_heta 轉接
    ├── data/{net,pos}/                        ← 輸入資料
    ├── _cache/                                ← 共用指紋快取
    ├── src/                                   ← 各實驗原始程式碼
    │   ├── exp4_smallworld_r1_boxplot/exp4_smallworld_r1_boxplot.py
    │   ├── exp5_real16_r1_boxplot/exp5_real16_r1_boxplot.py
    │   ├── ...
    │   └── table2_partition_metrics/table2_partition_metrics.py
    └── experiment/                            ← 各次執行結果
        └── experiment_<YYYYMMDD_HHMMSS>/
            ├── fig4.png ~ fig12.png
            ├── table1.csv
            ├── table2.csv, table2_summary.txt
            └── run_info.txt                   ← 本次執行參數與結果

輸出位置
--------
每次執行 cos_auto.py 會建立單一結果資料夾（比照原 auto.py）：
    HETA/experiment/experiment_<時間戳>/
所有啟用實驗的 fig*.png / table*.csv 都寫入其中（不再各自分散）。
執行參數與結果亦寫入同資料夾的 run_info.txt。
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
# 各實驗的原始程式碼資料夾統一收納於 HETA/src/ 之下
SRC_DIR = "src"
# 各次執行的結果（圖、表、run_info.txt）統一輸出到 HETA/experiment/experiment_<時間戳>/
EXP_DIR = "experiment"
# cos-HETA：新演算法所在 Python 套件（分母由 min 改為幾何平均）
COS_HETA_PKG = "cos_heta"
COS_HETA_DIR = os.path.join(ROOT, COS_HETA_PKG)
# cos-HETA：執行階段自動產生的 shim 套件，讓子程序的 `from heta import HETA`
#          被攔截並轉用 cos_heta 的實作。shim 會放在 PYTHONPATH 最前方。
COS_SHIM_DIR = os.path.join(ROOT, "_cos_shim")

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
def _ensure_shim():
    """建立 _cos_shim/heta/ 轉接套件。

    src/ 下各實驗腳本統一以 `from heta import HETA as heta` 取用演算法。
    為了在不修改任何 src/ 腳本的前提下改跑 cos_heta，本函式在
    `_cos_shim/heta/` 產生一個極薄的套件，其 __init__.py 內容會把
    `cos_heta.HETA`（以及可得的 `cos_heta.HETA2`）以 `heta.HETA` / `heta.HETA2`
    的名義註冊到 sys.modules。只要把 `_cos_shim/` 放在子程序 PYTHONPATH
    的最前方，子腳本 `from heta import HETA` 就會先命中這個 shim 並改
    走 cos_heta 版本（原本的 heta/ 套件在搜尋順序中被 shim 遮蔽）。

    回傳：shim 根目錄路徑 COS_SHIM_DIR（放入 PYTHONPATH 的那一級）。
    """
    pkg_dir = os.path.join(COS_SHIM_DIR, "heta")
    os.makedirs(pkg_dir, exist_ok=True)

    init_path = os.path.join(pkg_dir, "__init__.py")
    init_content = (
        "# -*- coding: utf-8 -*-\n"
        "# AUTO-GENERATED by cos_auto.py；請勿手動修改。\n"
        "# 作用：讓 `from heta import HETA as heta` 在 cos_auto.py 啟動的子程序中\n"
        "#       改用 cos_heta 套件的實作（分母由 min 改為幾何平均的 cos-HETA）。\n"
        "import sys as _sys\n"
        "from cos_heta import HETA  # 必備：實驗腳本使用的主模組\n"
        "try:\n"
        "    from cos_heta import HETA2  # 選用：若 cos_heta 有 HETA2 子模組則一併轉接\n"
        "except Exception:\n"
        "    HETA2 = None\n"
        "_sys.modules['heta.HETA'] = HETA\n"
        "if HETA2 is not None:\n"
        "    _sys.modules['heta.HETA2'] = HETA2\n"
    )
    with open(init_path, "w", encoding="utf-8") as f:
        f.write(init_content)

    # 確保 shim 資料夾本身不被誤認為普通套件（放一個說明檔）
    readme_path = os.path.join(COS_SHIM_DIR, "README.txt")
    if not os.path.exists(readme_path):
        with open(readme_path, "w", encoding="utf-8") as f:
            f.write(
                "此資料夾由 cos_auto.py 自動產生，用於在子程序執行階段將\n"
                "`from heta import ...` 轉接到 cos_heta 套件。請勿手動修改。\n"
            )

    return COS_SHIM_DIR


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


def _run_one(key, folder, script, desc, times, output_dir):
    """以 subprocess 執行單一實驗。

    程式碼位於 HETA/src/<folder>/；輸出統一寫到 output_dir（由 main() 一次
    決定並共用給所有實驗）。輸出目的地透過環境變數 HETA_OUTPUT_DIR 傳遞給
    子腳本；子腳本內會 os.chdir(HETA_OUTPUT_DIR) 讓 plt.savefig / to_csv 落點
    一致。
    """
    folder_path = os.path.join(ROOT, SRC_DIR, folder)
    script_path = os.path.join(folder_path, script)
    if not os.path.exists(script_path):
        print(f"❌ 找不到腳本：{script_path}")
        return False, 0.0

    print("\n" + "─" * 90)
    print(f"▶  [{key}] {desc}")
    print(f"   程式碼目錄：{folder_path}")
    print(f"   輸出目錄  ：{output_dir}")
    print(f"   執行腳本  ：{script}")
    print(f"   隨機化次數 HETA_TIMES = {times}")
    print("─" * 90)

    env = os.environ.copy()
    env["HETA_TIMES"] = str(times)
    env["HETA_OUTPUT_DIR"] = output_dir
    # 子腳本寫的是 `from heta import HETA as heta`。為了在不修改子腳本的
    # 前提下改跑 cos_heta，我們把 _cos_shim/ 放在 PYTHONPATH 最前方；
    # 再放 ROOT 使 `import cos_heta` 能找到 cos_heta/ 套件。順序：
    #   _cos_shim/  →  先命中 `heta` shim（重新導向到 cos_heta）
    #   ROOT/       →  提供 cos_heta（供 shim 內部 `from cos_heta import HETA`）
    # 原本的 ROOT/heta/ 仍在，但已被 shim 優先遮蔽。
    py_path = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = os.pathsep.join([COS_SHIM_DIR, ROOT]) + (
        os.pathsep + py_path if py_path else ""
    )

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


def _write_run_log(output_dir, start_dt, end_dt, grand_elapsed, args, summary):
    """
    將本次執行的參數與時間寫入 <output_dir>/run_info.txt
    （output_dir 為 HETA/experiment/experiment_<timestamp>/，由 main() 統一決定）
    run_info.txt 與本次所有 fig*.png / table*.csv 放在同一資料夾，方便歸檔比對。

    參數：
        output_dir     str      本次執行結果資料夾
        start_dt       datetime 執行開始時間
        end_dt         datetime 執行結束時間
        grand_elapsed  float    總耗時（秒）
        args           argparse.Namespace  CLI 參數
        summary        list[(key, ok, elapsed)]  各實驗結果
    """
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, "run_info.txt")

    lines = []
    lines.append("=" * 72)
    lines.append("  cos-HETA 論文實驗自動化（cos_auto.py） — 本次執行參數與結果")
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
    lines.append(f"Python    ：{sys.executable}")
    lines.append(f"HETA ROOT ：{ROOT}")
    lines.append(f"SRC_DIR   ：{os.path.join(ROOT, SRC_DIR)}")
    lines.append(f"演算法套件：{COS_HETA_DIR}  （cos-HETA，分母＝幾何平均）")
    lines.append(f"Shim 路徑 ：{COS_SHIM_DIR}  （heta→cos_heta 轉接；PYTHONPATH 最前）")
    lines.append(f"輸出目錄  ：{output_dir}")
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
        description="cos-HETA 論文實驗一鍵自動化（10 實驗 + Table 1；使用 cos_heta 套件）",
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

    # cos-HETA：執行前先產生（或更新）heta→cos_heta 的 shim 套件。
    shim_dir = _ensure_shim()

    print("\n" + "█" * 90)
    print("█" + " cos-HETA 論文實驗自動化主控（cos_auto.py） ".center(88) + "█")
    print("█" * 90)
    print(f"  使用演算法套件：{COS_HETA_DIR}")
    print(f"  Shim 路徑      ：{shim_dir}  （已注入 PYTHONPATH 最前方）")
    _print_list()

    enabled = [(k, f, s, d) for (k, f, s, d, _, _) in EXPERIMENTS if TOGGLE.get(k, False)]
    if not enabled:
        print("⚠  沒有任何實驗被啟用，結束。")
        return

    # 建立本次執行的結果資料夾 HETA/experiment/experiment_<時間戳>/
    start_dt = datetime.now()
    run_ts = start_dt.strftime("%Y%m%d_%H%M%S")
    run_output_dir = os.path.join(ROOT, EXP_DIR, f"experiment_{run_ts}")
    os.makedirs(run_output_dir, exist_ok=True)

    print(f"準備執行 {len(enabled)} 個實驗，HETA_TIMES = {args.times}")
    print(f"本次結果輸出目錄：{run_output_dir}\n")

    summary = []
    grand_t0 = time.time()
    for key, folder, script, desc in enabled:
        ok, elapsed = _run_one(key, folder, script, desc, args.times, run_output_dir)
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

    # 將本次執行的參數與時間寫入 <run_output_dir>/run_info.txt
    try:
        log_path = _write_run_log(
            run_output_dir, start_dt, end_dt, grand_elapsed, args, summary
        )
        print(f"本次執行紀錄已寫入：{log_path}")
        print(f"本次所有圖表與紀錄：{run_output_dir}\n")
    except Exception as e:
        print(f"⚠  寫入執行紀錄失敗：{e}\n")


if __name__ == "__main__":
    main()
