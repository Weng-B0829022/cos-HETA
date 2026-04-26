/*
 * threshold.c — 外部閾值（T_E^k）與內部閾值（T_I^k）計算
 * =========================================================
 *
 * 對應論文：
 *   公式 (5)：T_E^k = Mean_E^k(RG) + 2 × SD_E^k(RG)
 *   公式 (6)：T_I^k = Mean_I^k(candidates) - SD_I^k(candidates)
 *
 * ─────────────────────────────────────────────────────────
 * 重要：本檔案已對齊 Python HETA 的計算邏輯
 * ─────────────────────────────────────────────────────────
 *
 * Python HETA 對 T_E^k 的計算方式：
 *   1. 對每個隨機網路 RG_r，計算其所有邊在每層的 ratio 的 (μ_r^k, σ_r^k)
 *      （population std，ddof=0）
 *   2. T_E^k = (1/N) Σ_r μ_r^k + 2 × (1/N) Σ_r σ_r^k
 *   3. 若 T_E^k > 1.0，clamp 至 1.0
 *
 * Python 在每個圖只做一次 compute_link_property（產生所有邊、所有層的
 * ratio），我們也比照辦理：對每個 random network 只跑一次
 * compute_all_link_ratios，然後對每層分別做 Welford 統計。
 *
 * 主要 API：
 *   - compute_external_thresholds(g, layers, n_random, out_thresholds)
 *       一次 pass 算完所有層的 T_E。建議直接使用此 API。
 *   - compute_external_threshold(g, k, n_random)
 *       單層便利 wrapper（內部仍呼叫 thresholds 版），保留供舊程式碼相容。
 *   - compute_internal_threshold(ratios, n)
 *       T_I^k = mean - sd（population std，clamp ≥ 0）。
 *
 * 仍保留的「不影響邏輯」效能優化：
 *   - OpenMP 並行化 N 個隨機網路（彼此獨立）
 *   - PCG32 RNG（thread-safe，每 thread 獨立 seed）
 *   - 預分配 adj_bs / Bitset workspace（switching_randomize 內部）
 *   - compute_all_link_ratios 內部用 ego table 預先計算 + bitset 交集
 *
 * 已移除（會改動數學結果的優化）：
 *   - 跨網路池化 Welford 統計（pooled mean/std vs. mean-of-per-network-stats）
 *   - BFS 排除 t 整個節點的 ring 計算（與 Python 排除「s-t 直接邊」不同語意）
 */

 #include "heta.h"

 #include <stdlib.h>   /* malloc, free   */
 #include <stdio.h>    /* fprintf        */
 #include <math.h>     /* sqrt           */
 #include <stdint.h>   /* uint64_t       */
 #include <string.h>   /* memset         */

 /* 若有 OpenMP，引入標頭 */
 #ifdef _OPENMP
 #  include <omp.h>
 #endif

 /* ────────────────────────────────────────────────────────────
  * welford_layer — 對某一層的所有邊 ratio 做 Welford 線上統計
  *
  * 等同於 numpy.std(ratios, ddof=0) 加 numpy.mean。
  * 這裡把 mean/std 拆兩個輸出，方便外部累積。
  * ──────────────────────────────────────────────────────────── */
 static void welford_layer(const double *ratios_table, int m, int layers,
                           int k, double *out_mean, double *out_std) {
     long   count = 0;
     double mean  = 0.0;
     double M2    = 0.0;

     for (int e = 0; e < m; e++) {
         double r = ratios_table[(size_t)e * layers + (k - 1)];
         count++;
         double delta  = r - mean;
         mean         += delta / (double)count;
         double delta2 = r - mean;
         M2           += delta * delta2;
     }

     *out_mean = mean;
     *out_std  = (count > 0) ? sqrt(M2 / (double)count) : 0.0;
 }

 /* ════════════════════════════════════════════════════════════
  * compute_external_thresholds — 一次 pass 計算所有層的 T_E^k
  *
  * 對齊 Python HETA。out_thresholds 大小至少 layers+1：
  *   out_thresholds[k] = T_E^k，k = 1..layers
  *   out_thresholds[0] 不使用
  * ════════════════════════════════════════════════════════════ */
 void compute_external_thresholds(const Graph *g, int layers, int n_random,
                                  double *out_thresholds) {
     int n = g->n;

     /* 邊界保護 */
     for (int k = 0; k <= layers; k++) out_thresholds[k] = 0.0;
     if (n_random <= 0 || layers < 1 || g->m == 0) return;

     /* per-network 統計：per_net_mean[r * layers + (k-1)] 等等 */
     size_t stat_count = (size_t)n_random * (size_t)layers;
     double *per_net_mean = (double *)malloc(stat_count * sizeof(double));
     double *per_net_std  = (double *)malloc(stat_count * sizeof(double));
     if (!per_net_mean || !per_net_std) {
         fprintf(stderr, "[threshold] per-network buffer malloc failed\n");
         exit(1);
     }
     for (size_t i = 0; i < stat_count; i++) {
         per_net_mean[i] = 0.0;
         per_net_std[i]  = 0.0;
     }

#ifdef _OPENMP
#pragma omp parallel
     {
         /* thread-local：switching 用的鄰接矩陣 bitset */
         Bitset **adj_bs = adj_bs_alloc(n);
         int tid = omp_get_thread_num();

#pragma omp for schedule(static)
         for (int r = 0; r < n_random; r++) {
             uint64_t seed = ((uint64_t)tid * 1000003ULL)
                           + ((uint64_t)r   * 2654435761ULL);

             /* 建立隨機網路 */
             Graph *gr = graph_copy(g);
             switching_randomize(gr, seed, adj_bs);

             /* 用 Python 風格演算法計算所有邊、所有層的 ratio */
             double *ratios = compute_all_link_ratios(gr, layers);

             if (ratios) {
                 int gr_m = gr->m;
                 for (int k = 1; k <= layers; k++) {
                     double m_r, s_r;
                     welford_layer(ratios, gr_m, layers, k, &m_r, &s_r);
                     per_net_mean[(size_t)r * layers + (k - 1)] = m_r;
                     per_net_std [(size_t)r * layers + (k - 1)] = s_r;
                 }
                 free(ratios);
             }

             graph_free(gr);
         }

         adj_bs_free(adj_bs, n);
     }
#else
     /* ── 無 OpenMP 時的串列版本 ── */
     Bitset **adj_bs = adj_bs_alloc(n);

     for (int r = 0; r < n_random; r++) {
         uint64_t seed = (uint64_t)r * 2654435761ULL ^ 0xDEADBEEFULL;
         Graph *gr = graph_copy(g);
         switching_randomize(gr, seed, adj_bs);

         double *ratios = compute_all_link_ratios(gr, layers);
         if (ratios) {
             int gr_m = gr->m;
             for (int k = 1; k <= layers; k++) {
                 double m_r, s_r;
                 welford_layer(ratios, gr_m, layers, k, &m_r, &s_r);
                 per_net_mean[(size_t)r * layers + (k - 1)] = m_r;
                 per_net_std [(size_t)r * layers + (k - 1)] = s_r;
             }
             free(ratios);
         }
         graph_free(gr);
     }

     adj_bs_free(adj_bs, n);
#endif

     /* ── 對所有隨機網路的 μ_r^k 取平均、對 σ_r^k 取平均 ── */
     for (int k = 1; k <= layers; k++) {
         double sum_means = 0.0;
         double sum_stds  = 0.0;
         for (int r = 0; r < n_random; r++) {
             sum_means += per_net_mean[(size_t)r * layers + (k - 1)];
             sum_stds  += per_net_std [(size_t)r * layers + (k - 1)];
         }
         double avg_mean = sum_means / (double)n_random;
         double avg_std  = sum_stds  / (double)n_random;
         double t_e      = avg_mean + 2.0 * avg_std;
         if (t_e > 1.0) t_e = 1.0;
         out_thresholds[k] = t_e;
     }

     free(per_net_mean);
     free(per_net_std);
 }

 /* ════════════════════════════════════════════════════════════
  * compute_external_threshold — 單層便利 wrapper
  *
  * 內部呼叫 compute_external_thresholds 計算 1..k 全部，回傳第 k 個。
  * 注意：若需要連續多個層的門檻，請直接使用 compute_external_thresholds，
  * 否則每呼叫一次都會重新生成所有隨機網路。
  * ════════════════════════════════════════════════════════════ */
 double compute_external_threshold(const Graph *g, int k, int n_random) {
     if (k < 1) return 0.0;
     double *thresholds = (double *)malloc((size_t)(k + 1) * sizeof(double));
     if (!thresholds) {
         fprintf(stderr, "[threshold] thresholds malloc failed\n");
         exit(1);
     }
     compute_external_thresholds(g, k, n_random, thresholds);
     double result = thresholds[k];
     free(thresholds);
     return result;
 }

 /* ════════════════════════════════════════════════════════════
  * compute_internal_threshold — 計算內部閾值 T_I^k
  *
  * 論文公式 (6)：T_I^k = Mean_I^k - SD_I^k（population std，ddof=0）
  * 對齊 Python HETA：T_I 不可為負，<0 時 clamp 為 0.0
  *
  * 輸入是候選邊的 R^k 值陣列（已在主演算法中計算好）。
  *
  * 參數：
  *   ratios : 候選邊的 R^k 值，長度 n
  *   n      : 候選邊數量
  * ════════════════════════════════════════════════════════════ */
 double compute_internal_threshold(const double *ratios, int n) {
     if (n == 0) return 0.0;

     /* Welford 線上計算（等同 numpy.std(ddof=0)） */
     long   cnt   = 0;
     double mean  = 0.0;
     double M2    = 0.0;
     for (int i = 0; i < n; i++) {
         cnt++;
         double delta  = ratios[i] - mean;
         mean         += delta / (double)cnt;
         double delta2 = ratios[i] - mean;
         M2           += delta * delta2;
     }

     double sd = (cnt > 0) ? sqrt(M2 / (double)cnt) : 0.0;

     /* T_I^k = mean - SD（論文公式 (6)，clamp ≥ 0） */
     double t_i = mean - sd;
     return (t_i < 0.0) ? 0.0 : t_i;
 }
