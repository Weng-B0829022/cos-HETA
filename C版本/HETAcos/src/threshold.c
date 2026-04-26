/*
 * threshold.c — 外部閾值（T_E^k）與內部閾值（T_I^k）計算
 * =========================================================
 *
 * 對應論文：
 *   公式 (5)：T_E^k = Mean_E^k(RG) + 2 × SD_E^k(RG)
 *   公式 (6)：T_I^k = Mean_I^k(candidates) - SD_I^k(candidates)
 *
 * 主要優化：
 *
 * 1. OpenMP 並行化 1000 個隨機網路
 *    - 1000 個網路彼此完全獨立，是天然的並行工作
 *    - 每個 thread 持有自己的：RNG 狀態、BFSBuffer、Bitset workspace
 *    - 使用 Welford online algorithm 累計 mean/variance，每個 thread
 *      維護自己的 (count, mean, M2)，最後用 parallel reduction 合併
 *
 * 2. Welford Online Algorithm 計算均值與變異數
 *    - Python 原版：先把所有 ratio 存進 list（O(|RG|×|E|) 記憶體）
 *      再 np.mean / np.std
 *    - Welford：一次掃過即得 mean 和 variance，記憶體 O(1)
 *    - 數值穩定性優於兩遍演算法（避免大數相減的精度損失）
 *    - 參考：Welford (1962), Technometrics 4(3):419-420
 *
 * 3. 並行 reduction（合併多個 Welford 狀態）
 *    - 每個 thread 各自維護 (n_k, mean_k, M2_k)
 *    - 最後用 Chan's parallel variance 公式合併所有 thread 的結果
 *    - 完全無鎖（lock-free），只在最後合併時需要串列操作
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
  * Welford 線上統計狀態
  * 維護：樣本數 n、目前均值 mean、累計平方偏差 M2
  * 關係：variance = M2 / n（population variance）
  *       std      = sqrt(M2 / n)
  * ──────────────────────────────────────────────────────────── */
 typedef struct {
     long   n;       /* 已處理的樣本數    */
     double mean;    /* 目前的滾動均值    */
     double M2;      /* 累計平方偏差      */
 } Welford;
 
 /* 初始化 Welford 狀態為空 */
 static inline void welford_init(Welford *w) {
     w->n    = 0;
     w->mean = 0.0;
     w->M2   = 0.0;
 }
 
 /* 加入一個新樣本值 x，更新滾動 mean 和 M2 */
 static inline void welford_update(Welford *w, double x) {
     w->n++;
     double delta  = x - w->mean;          /* 新值與目前均值的差 */
     w->mean      += delta / (double)w->n; /* 更新均值           */
     double delta2 = x - w->mean;          /* 更新後的差（注意是新 mean）*/
     w->M2        += delta * delta2;       /* 累計平方偏差       */
 }
 
 /*
  * 合併兩個 Welford 狀態（Chan's parallel algorithm）
  * 用於將各 thread 的區域統計量合併成全域統計量
  *
  * 公式（Chan et al. 1979）：
  *   n_combined = na + nb
  *   delta      = mean_b - mean_a
  *   M2_combined = M2_a + M2_b + delta^2 * na * nb / n_combined
  */
 static void welford_merge(Welford *dst, const Welford *src) {
     if (src->n == 0) return;              /* src 為空，不需合併          */
     if (dst->n == 0) { *dst = *src; return; } /* dst 為空，直接複製    */
 
     long   n_combined = dst->n + src->n;
     double delta      = src->mean - dst->mean;   /* 兩均值之差 */
 
     /* 合併均值（加權平均）*/
     dst->mean = (dst->mean * dst->n + src->mean * src->n) / (double)n_combined;
 
     /* 合併平方偏差（Chan's formula）*/
     dst->M2  += src->M2
               + delta * delta * ((double)dst->n * (double)src->n
                                  / (double)n_combined);
     dst->n    = n_combined;
 }
 
 /* ════════════════════════════════════════════════════════════
  * compute_external_threshold — 計算第 k 層外部閾值 T_E^k
  *
  * 論文公式 (5)：T_E^k = Mean_E^k(RG) + 2 × SD_E^k(RG)
  *
  * 步驟：
  *   1. 生成 n_random 個隨機化網路（OpenMP 並行）
  *   2. 對每個隨機網路的每條邊計算 R^k
  *   3. 用 Welford algorithm 累計 mean 和 std
  *   4. 回傳 mean + 2×std
  *
  * 參數：
  *   g        : 原始圖（不會被修改）
  *   k        : 層數
  *   n_random : 隨機化網路數量（論文建議 1000）
  * ════════════════════════════════════════════════════════════ */
 double compute_external_threshold(const Graph *g, int k, int n_random) {
     int n = g->n;
 
     /* 全域 Welford 狀態（所有 thread 共用，最後合併）*/
     Welford global;
     welford_init(&global);
 
     /*
      * OpenMP parallel for：
      *   - 每個 iteration 處理一個隨機化網路
      *   - reduction 子句自動合併各 thread 的結果（此處手動合併）
      *   - 每個 thread 有自己的 local_wf, buf, bitsets
      */
 #ifdef _OPENMP
 #pragma omp parallel
     {
         /* ── thread-local 資源分配 ── */
 
         /* 每個 thread 自己的 Welford 狀態 */
         Welford local_wf;
         welford_init(&local_wf);
 
         /* thread-local BFS 緩衝區（避免 thread 之間共享記憶體）*/
         BFSBuffer *buf   = bfs_buffer_alloc(n);
 
         /* thread-local Bitset 工作空間（4 個：k 層和 k-1 層各兩個）*/
         Bitset *bs_k_u   = bs_alloc(n);
         Bitset *bs_k_v   = bs_alloc(n);
         Bitset *bs_k1_u  = bs_alloc(n);
         Bitset *bs_k1_v  = bs_alloc(n);
 
         /*
          * 預分配 adj_bs：n 個 Bitset，整個 thread 的所有迭代共用
          * 每次迭代只做 memset 清零，不重複 malloc/free
          * 對大型網路（n=1000）可省去 1000×n 次 malloc，加速約 50 倍
          */
         Bitset **adj_bs  = adj_bs_alloc(n);
 
         /* 用 thread id 和迭代偏移生成唯一 seed，確保各 thread RNG 獨立 */
         int tid = omp_get_thread_num();
 
         /* 靜態排程：每個 thread 分配固定數量的隨機網路 */
 #pragma omp for schedule(static)
         for (int r = 0; r < n_random; r++) {
             /* 為此次隨機化生成唯一 seed（結合 thread id 和迭代號）*/
             uint64_t seed = ((uint64_t)tid * 1000003ULL)
                           + ((uint64_t)r   * 2654435761ULL);
 
             /* 複製原圖（深拷貝，O(n+m)，但這是連續 memcpy，極快）*/
             Graph *gr = graph_copy(g);
 
             /* 對複本進行 switching randomization
              * 傳入預分配的 adj_bs，函式內部只做清零+填入，不 malloc */
             switching_randomize(gr, seed, adj_bs);
 
             /* 計算此隨機網路每條邊的 R^k，用 Welford 累計 */
             /* 遍歷 CSR 鄰接表，只取 u < v 方向（無向圖去重）*/
             for (int u = 0; u < gr->n; u++) {
                 for (int i = gr->offset[u]; i < gr->offset[u + 1]; i++) {
                     int v = gr->adj_list[i];
                     if (u >= v) continue;   /* 只處理 u < v，避免每條邊計算兩次 */
 
                     /* 計算邊 (u,v) 的第 k 層共同鄰居比例 */
                     double ratio = compute_cnr(gr, u, v, k, buf,
                                                bs_k_u, bs_k_v,
                                                bs_k1_u, bs_k1_v);
 
                     /* 更新 thread-local Welford 狀態（線上計算 mean/var）*/
                     welford_update(&local_wf, ratio);
                 }
             }
 
             graph_free(gr);   /* 釋放隨機化圖 */
         }
 
         /* ── 合併 thread-local 結果到全域 Welford（需要 critical section）── */
 #pragma omp critical
         {
             welford_merge(&global, &local_wf);
         }
 
         /* 釋放 thread-local 資源（adj_bs 在此統一釋放，不是在迭代內）*/
         adj_bs_free(adj_bs, n);
         bfs_buffer_free(buf);
         bs_free(bs_k_u);
         bs_free(bs_k_v);
         bs_free(bs_k1_u);
         bs_free(bs_k1_v);
     }
 #else
     /* ── 無 OpenMP 時的串列版本 ── */
     BFSBuffer *buf   = bfs_buffer_alloc(n);
     Bitset    *bs_k_u  = bs_alloc(n);
     Bitset    *bs_k_v  = bs_alloc(n);
     Bitset    *bs_k1_u = bs_alloc(n);
     Bitset    *bs_k1_v = bs_alloc(n);
 
     /* 預分配 adj_bs，1000 次迭代共用，每次只 memset 清零 */
     Bitset   **adj_bs  = adj_bs_alloc(n);
 
     for (int r = 0; r < n_random; r++) {
         uint64_t seed = (uint64_t)r * 2654435761ULL ^ 0xDEADBEEFULL;
         Graph *gr = graph_copy(g);
         switching_randomize(gr, seed, adj_bs);   /* 傳入預分配的 adj_bs */
 
         /* 遍歷隨機圖的所有邊並計算 R^k */
         for (int u = 0; u < gr->n; u++) {
             for (int i = gr->offset[u]; i < gr->offset[u + 1]; i++) {
                 int v = gr->adj_list[i];
                 if (u >= v) continue;   /* 只處理 u < v */
 
                 double ratio = compute_cnr(gr, u, v, k, buf,
                                            bs_k_u, bs_k_v,
                                            bs_k1_u, bs_k1_v);
                 welford_update(&global, ratio);
             }
         }
         graph_free(gr);
     }
 
     adj_bs_free(adj_bs, n);   /* 迴圈結束後統一釋放 */
     bfs_buffer_free(buf);
     bs_free(bs_k_u);
     bs_free(bs_k_v);
     bs_free(bs_k1_u);
     bs_free(bs_k1_v);
 #endif
 
     /* 若沒有任何樣本（空圖），回傳 0 */
     if (global.n == 0) return 0.0;
 
     /* 計算全域標準差：SD = sqrt(M2 / n) */
     double mean = global.mean;
     double sd   = (global.n > 1) ? sqrt(global.M2 / (double)global.n) : 0.0;
 
     /* T_E^k = mean + 2 × SD（論文公式 (5)）*/
     return mean + 2.0 * sd;
 }
 
 /* ════════════════════════════════════════════════════════════
  * compute_internal_threshold — 計算內部閾值 T_I^k
  *
  * 論文公式 (6)：T_I^k = Mean_I^k - SD_I^k
  *
  * 輸入是候選邊的 R^k 值陣列（已在主演算法中計算好）。
  * 使用 Welford 演算法（雖然資料已在陣列中，但保持一致的實作風格）
  *
  * 參數：
  *   ratios : 候選邊的 R^k 值，長度 n
  *   n      : 候選邊數量
  * ════════════════════════════════════════════════════════════ */
 double compute_internal_threshold(const double *ratios, int n) {
     if (n == 0) return 0.0;   /* 無候選邊，閾值為 0 */
 
     Welford w;
     welford_init(&w);
 
     /* 逐一加入每個候選邊的比例值 */
     for (int i = 0; i < n; i++) {
         welford_update(&w, ratios[i]);
     }
 
     /* 計算標準差 */
     double sd = (w.n > 1) ? sqrt(w.M2 / (double)w.n) : 0.0;
 
     /* T_I^k = mean - SD（論文公式 (6)）*/
     return w.mean - sd;
 }