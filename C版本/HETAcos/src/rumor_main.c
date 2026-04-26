/*
 * rumor_main.c — 謠言擴散模型主程式（C 實作）
 * ================================================
 *
 * 對應 Python main.py 的完整流程：
 *   [1/4] 建立 WS 和 NW 小世界網路
 *   [2/4] HETA 邊分類
 *   [3/4] 選擇初始投放節點（Global Bridge 兩端）
 *   [4/4] 執行謠言擴散模擬
 *
 * 使用方式：
 *   ./rumor_spread [n_random] [n_nodes] [days] [k] [rewire] [seeds]
 *                  [seed_rng] [decay_rate] [zero_day] [willingness]
 *                  [bond] [local] [global] [silk]
 *
 *   所有參數皆為可選，未指定則使用 rumor.h 中的預設值
 *
 * 範例：
 *   ./rumor_spread                              → 全部使用預設值
 *   ./rumor_spread 1 100 500                    → n_random=1, 100節點, 500天
 *   ./rumor_spread 100 100 60 6 0.7 5 40 0.75 5 0.1 0.4 0.2 0.7 1.0
 *
 * 輸出：
 *   result/rumor_spread_result.csv  — 每日感染比例
 *   終端文字摘要
 */

#include "rumor.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

/* ────────────────────────────────────────────────────────────
 * 列印邊類型分布
 * ──────────────────────────────────────────────────────────── */
static void print_edge_distribution(const char *label, const HetaResult *heta) {
    HetaSummary s = heta_summarize(heta);
    fprintf(stdout, "\n  %s 邊類型分布 (總邊數=%d):\n", label, heta->m);
    fprintf(stdout, "    %-8s: %4d 條  (%5.1f%%)\n", "bond",   s.n_bond,   s.pct_bond);
    fprintf(stdout, "    %-8s: %4d 條  (%5.1f%%)\n", "local",  s.n_local,  s.pct_local);
    fprintf(stdout, "    %-8s: %4d 條  (%5.1f%%)\n", "global", s.n_global, s.pct_global);
    fprintf(stdout, "    %-8s: %4d 條  (%5.1f%%)\n", "silk",   s.n_silk,   s.pct_silk);
}

/* ────────────────────────────────────────────────────────────
 * 計時輔助
 * ──────────────────────────────────────────────────────────── */
static double get_time_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

/* ════════════════════════════════════════════════════════════
 * main
 * ════════════════════════════════════════════════════════════ */
int main(int argc, char *argv[]) {
    /* ── 參數設定 ── */
    RumorParams params = rumor_default_params();

    /*
     * 命令列參數順序（全部可選，未指定則使用預設值）：
     *
     *  位置  參數名          範例值   說明
     *  ────  ──────────────  ──────  ─────────────────────────
     *   1    n_random        100     HETA 隨機化網路數量
     *   2    n_nodes         100     網路節點數量
     *   3    days            60      模擬總天數
     *   4    k_neighbors     6       每節點鄰居數（偶數）
     *   5    rewire_prob     0.7     WS 重連機率 / NW 捷徑比例
     *   6    n_seeds         5       種子節點數
     *   7    random_seed     40      隨機種子
     *   8    decay_rate      0.75    每日無傳入信心度保留率
     *   9    decay_zero_day  5       連續無傳入幾天後歸零
     *  10    willingness     0.1     全節點散播意願
     *  11    conf_bond       0.40    bond 邊信心度貢獻率
     *  12    conf_local      0.20    local bridge 信心度貢獻率
     *  13    conf_global     0.70    global bridge 信心度貢獻率
     *  14    conf_silk       1.00    silk 邊信心度貢獻率
     *  15    infected_decay  0.30    感染節點專用衰減率（越低消亡越快）
     */
    if (argc >= 2)  params.n_random_nets = atoi(argv[1]);
    if (argc >= 3)  params.n_nodes       = atoi(argv[2]);
    if (argc >= 4)  params.days          = atoi(argv[3]);
    if (argc >= 5)  params.k_neighbors   = atoi(argv[4]);
    if (argc >= 6)  params.rewire_prob   = atof(argv[5]);
    if (argc >= 7)  params.n_seeds       = atoi(argv[6]);
    if (argc >= 8)  params.random_seed   = atoi(argv[7]);
    if (argc >= 9)  params.decay_rate    = atof(argv[8]);
    if (argc >= 10) params.decay_zero_day= atoi(argv[9]);
    if (argc >= 11) params.willingness   = atof(argv[10]);
    if (argc >= 12) params.conf_bond     = atof(argv[11]);
    if (argc >= 13) params.conf_local    = atof(argv[12]);
    if (argc >= 14) params.conf_global   = atof(argv[13]);
    if (argc >= 15) params.conf_silk     = atof(argv[14]);
    if (argc >= 16) params.infected_decay = atof(argv[15]);

    fprintf(stdout, "==============================================================\n");
    fprintf(stdout, "  謠言擴散模型（HETA 邊類型 + 固定衰減率信心度）— C 版\n");
    fprintf(stdout, "==============================================================\n");
    fprintf(stdout, "  節點數: %d, 鄰居數: %d, 重連率: %.2f\n",
            params.n_nodes, params.k_neighbors, params.rewire_prob);
    fprintf(stdout, "  投放機制: 所有 Global Bridge 端點（自動決定種子數）\n");
    fprintf(stdout, "  天數: %d\n", params.days);
    fprintf(stdout, "  衰減規則: 連續無傳入每天 ×%.2f，%d天後強制歸0\n",
            params.decay_rate, params.decay_zero_day);
    fprintf(stdout, "  有傳入時：累加信心度並重置衰減計時器\n");
    fprintf(stdout, "  HETA n_random: %d\n", params.n_random_nets);
    fprintf(stdout, "\n");

    double t_total_start = get_time_sec();

    /* ═══════════════════════════════════════════════════════════
     * [1/4] 建立 WS 和 NW 小世界網路
     * ═══════════════════════════════════════════════════════════ */
    fprintf(stdout, "[1/4] 建立 WS 和 NW 小世界網路...\n");
    double t0 = get_time_sec();

    Graph *G_ws = build_ws_network(params.n_nodes, params.k_neighbors,
                                    params.rewire_prob,
                                    (uint64_t)params.random_seed);
    Graph *G_nw = build_nw_network(params.n_nodes, params.k_neighbors,
                                    params.rewire_prob,
                                    (uint64_t)params.random_seed);

    double t1 = get_time_sec();
    fprintf(stdout, "  WS: %d 節點, %d 邊\n", G_ws->n, G_ws->m);
    fprintf(stdout, "  NW: %d 節點, %d 邊\n", G_nw->n, G_nw->m);
    fprintf(stdout, "  網路建立耗時: %.3f 秒\n", t1 - t0);

    /* ═══════════════════════════════════════════════════════════
     * [2/4] HETA 邊分類
     * ═══════════════════════════════════════════════════════════ */
    fprintf(stdout, "\n[2/4] WS 網路 HETA 邊分類...\n");
    t0 = get_time_sec();
    HetaResult *heta_ws = heta_run(G_ws, params.n_random_nets);
    t1 = get_time_sec();
    fprintf(stdout, "  WS HETA 耗時: %.3f 秒\n", t1 - t0);

    fprintf(stdout, "\n[2/4] NW 網路 HETA 邊分類...\n");
    t0 = get_time_sec();
    HetaResult *heta_nw = heta_run(G_nw, params.n_random_nets);
    t1 = get_time_sec();
    fprintf(stdout, "  NW HETA 耗時: %.3f 秒\n", t1 - t0);

    print_edge_distribution("WS", heta_ws);
    print_edge_distribution("NW", heta_nw);

    /* ── 邊類型比例摘要 ── */
    {
        HetaSummary ws_s = heta_summarize(heta_ws);
        HetaSummary nw_s = heta_summarize(heta_nw);
        fprintf(stdout, "\n── 邊類型比例摘要 ──\n");
        fprintf(stdout, "%-8s  %8s  %8s\n", "Type", "WS (%)", "NW (%)");
        fprintf(stdout, "------------------------------\n");
        fprintf(stdout, "%-8s  %7.1f%%  %7.1f%%\n", "bond",   ws_s.pct_bond,   nw_s.pct_bond);
        fprintf(stdout, "%-8s  %7.1f%%  %7.1f%%\n", "local",  ws_s.pct_local,  nw_s.pct_local);
        fprintf(stdout, "%-8s  %7.1f%%  %7.1f%%\n", "global", ws_s.pct_global, nw_s.pct_global);
        fprintf(stdout, "%-8s  %7.1f%%  %7.1f%%\n", "silk",   ws_s.pct_silk,   nw_s.pct_silk);
    }

    /* ═══════════════════════════════════════════════════════════
     * [3/4] 選擇初始投放節點（Global Bridge 兩端）
     * ═══════════════════════════════════════════════════════════ */
    fprintf(stdout, "\n[3/4] 選擇初始投放節點（Global Bridge 兩端）...\n");
    int ws_n_seeds, nw_n_seeds;
    int *seeds_ws = select_seeds_from_global(G_ws, heta_ws,
                                              params.n_seeds, &ws_n_seeds);
    int *seeds_nw = select_seeds_from_global(G_nw, heta_nw,
                                              params.n_seeds, &nw_n_seeds);

    fprintf(stdout, "  WS 種子節點:");
    for (int i = 0; i < ws_n_seeds; i++) fprintf(stdout, " %d", seeds_ws[i]);
    fprintf(stdout, "\n");
    fprintf(stdout, "  NW 種子節點:");
    for (int i = 0; i < nw_n_seeds; i++) fprintf(stdout, " %d", seeds_nw[i]);
    fprintf(stdout, "\n");

    /* ═══════════════════════════════════════════════════════════
     * [4/4] 執行謠言擴散模擬
     * ═══════════════════════════════════════════════════════════ */
    fprintf(stdout, "\n[4/4] 執行謠言擴散模擬...\n");
    t0 = get_time_sec();

    SimResult ws_result = simulate_rumor(G_ws, heta_ws, seeds_ws, ws_n_seeds, &params);
    SimResult nw_result = simulate_rumor(G_nw, heta_nw, seeds_nw, nw_n_seeds, &params);

    t1 = get_time_sec();
    fprintf(stdout, "  模擬耗時: %.3f 秒\n", t1 - t0);

    /* ── 列印最終結果 ── */
    double ws_final = ws_result.daily_ratio[ws_result.n_days - 1] * 100.0;
    double nw_final = nw_result.daily_ratio[nw_result.n_days - 1] * 100.0;

    fprintf(stdout, "\n  WS 最終感染率: %.1f%%", ws_final);
    if (ws_final < 1.0) fprintf(stdout, " (Vanished)");
    else if (ws_final > 80.0) fprintf(stdout, " (Epidemic)");
    else fprintf(stdout, " (Stable)");
    fprintf(stdout, "\n");

    fprintf(stdout, "  NW 最終感染率: %.1f%%", nw_final);
    if (nw_final < 1.0) fprintf(stdout, " (Vanished)");
    else if (nw_final > 80.0) fprintf(stdout, " (Epidemic)");
    else fprintf(stdout, " (Stable)");
    fprintf(stdout, "\n");

    /* ── 每 10 天摘要 ── */
    print_daily_summary(&ws_result, &nw_result, 10);

    /* ── 寫入 CSV ── */
#if defined(_WIN32) || defined(_WIN64)
    mkdir("result");
#else
    mkdir("result", 0755);
#endif

    write_result_csv("result/rumor_spread_result.csv",
                     &ws_result, &nw_result, &params);

    /* ── 總耗時 ── */
    double t_total = get_time_sec() - t_total_start;
    fprintf(stdout, "\n✅ 全部完成！總耗時: %.3f 秒\n", t_total);

    /* ── 清理記憶體 ── */
    free(ws_result.daily_ratio);
    free(nw_result.daily_ratio);
    free(seeds_ws);
    free(seeds_nw);
    heta_result_free(heta_ws);
    heta_result_free(heta_nw);
    graph_free(G_ws);
    graph_free(G_nw);

    return EXIT_SUCCESS;
}
