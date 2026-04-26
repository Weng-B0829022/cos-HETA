/*
 * heta.c — HETA 主演算法：邊分類流程
 * =====================================
 *
 * 對應論文 Algorithm 1（三步流程）：
 *
 *   Step 3.1：識別 silk links
 *             → 端點之一 degree == 1 的邊直接標記為 silk
 *
 *   Step 3.2：逐層識別 bond links 與 local bridges
 *             → 對每層 k（1 到 kmax）：
 *               a. 計算未分類邊的 R^k
 *               b. R^k >= T_E^k → bond link
 *               c. 其餘邊計算 T_I^k，R^k > T_I^k → kth-layer local bridge
 *               d. 其餘邊傳入下一層
 *
 *   Step 3.3：識別 global bridges
 *             → 所有層都未能分類的邊 = global bridge
 *
 * 資料流說明：
 *   - 所有邊儲存在 Edge 陣列中，用 active 旗標表示是否尚未分類
 *   - 每層迭代只處理 active == 1 的邊（跳過已分類的邊）
 *   - 外部閾值 T_E^k 預先計算好（最耗時的步驟，OpenMP 並行）
 *   - 內部閾值 T_I^k 在每層迭代內即時計算（O(candidate count)）
 */

 #include "heta.h"

 #include <stdlib.h>   /* malloc, free   */
 #include <stdio.h>    /* fprintf        */
 #include <string.h>   /* memset         */
 
 /* 安全 malloc */
 static void *xmalloc(size_t sz) {
     void *p = malloc(sz);
     if (!p) { fprintf(stderr, "[heta] malloc failed\n"); exit(1); }
     return p;
 }
 
 /* ════════════════════════════════════════════════════════════
  * heta_result_free — 釋放 HetaResult 結構
  * ════════════════════════════════════════════════════════════ */
 void heta_result_free(HetaResult *r) {
     if (!r) return;
     free(r->edges);     /* 釋放邊陣列         */
     free(r->ext_th);    /* 釋放外部閾值陣列   */
     free(r);            /* 釋放結構本體       */
 }
 
 /* ════════════════════════════════════════════════════════════
  * heta_run — 主演算法入口
  *
  * 參數：
  *   g        : CSR 格式圖（不會被修改）
  *   n_random : 隨機化網路數量（論文建議 1000）
  *
  * 回傳：HetaResult*（呼叫者負責 heta_result_free）
  * ════════════════════════════════════════════════════════════ */
 HetaResult *heta_run(const Graph *g, int n_random) {
     int n = g->n;
     int m = g->m;
 
     /* ── 分配結果結構 ── */
     HetaResult *result = (HetaResult *)xmalloc(sizeof(HetaResult));
     result->m          = m;
 
     /* ── 計算 kmax（論文公式 (4)）── */
     result->kmax = graph_kmax(g);
     int kmax     = result->kmax;
 
     fprintf(stdout, "[heta] n=%d, m=%d, kmax=%d, n_random=%d\n",
             n, m, kmax, n_random);
 
     /* ── 建立邊陣列，填入端點與初始狀態 ── */
     result->edges = (Edge *)xmalloc(m * sizeof(Edge));
     int ei = 0;
     for (int u = 0; u < n; u++) {
         for (int i = g->offset[u]; i < g->offset[u + 1]; i++) {
             int v = g->adj_list[i];
             if (u >= v) continue;   /* 只取 u < v 的方向，無向圖不重複 */
             result->edges[ei].u      = u;
             result->edges[ei].v      = v;
             result->edges[ei].type   = LINK_UNCLASSIFIED;  /* 初始未分類 */
             result->edges[ei].layer  = 0;
             result->edges[ei].active = 1;   /* 尚未分類（pass_flag = true）*/
             ei++;
         }
     }
     /* 確認實際邊數與宣告的 m 一致 */
     /* ei 應等於 m */
 
     /* ── Step 3.1：識別 silk links ── */
     /* silk link：端點之一的 degree == 1（對應論文定義）*/
     int n_silk = 0;
     for (int i = 0; i < m; i++) {
         int u = result->edges[i].u;
         int v = result->edges[i].v;
         if (g->degree[u] == 1 || g->degree[v] == 1) {
             result->edges[i].type   = LINK_SILK;
             result->edges[i].active = 0;   /* 標記為已分類（不再處理）*/
             n_silk++;
         }
     }
     fprintf(stdout, "[heta] Step 3.1: %d silk links identified\n", n_silk);
 
     /* ── Step 2：預先計算各層外部閾值（最耗時，OpenMP 並行）── */
     /* ext_th[0] 不使用，ext_th[1..kmax] 對應各層 */
     result->ext_th = (double *)xmalloc((kmax + 1) * sizeof(double));
     result->ext_th[0] = 0.0;   /* 佔位，不使用 */
 
     for (int k = 1; k <= kmax; k++) {
         fprintf(stdout, "[heta] Computing T_E^%d (n_random=%d)...\n",
                 k, n_random);
         result->ext_th[k] = compute_external_threshold(g, k, n_random);
         fprintf(stdout, "[heta] T_E^%d = %.6f\n", k, result->ext_th[k]);
     }
 
     /* ── Step 3.2：逐層識別 bond links 與 local bridges ── */
 
     /* 分配 thread-local 資源（此處為串列主迴圈，只需一套）*/
     BFSBuffer *buf    = bfs_buffer_alloc(n);
     Bitset    *bs_k_u  = bs_alloc(n);
     Bitset    *bs_k_v  = bs_alloc(n);
     Bitset    *bs_k1_u = bs_alloc(n);
     Bitset    *bs_k1_v = bs_alloc(n);
 
     /* 候選邊的 ratio 暫存陣列（最多 m 條）*/
     double *cand_ratios = (double *)xmalloc(m * sizeof(double));
     int    *cand_idx    = (int    *)xmalloc(m * sizeof(int));    /* 候選邊在 edges[] 的 index */
 
     for (int k = 1; k <= kmax; k++) {
         double T_E = result->ext_th[k];   /* 此層的外部閾值 */
         int    n_cand = 0;                 /* 候選邊數量（尚未分類，且 R^k < T_E）*/
         int    n_bond = 0;                 /* 此層識別出的 bond link 數 */
 
         /* ── 遍歷所有尚未分類的邊，計算 R^k ── */
         for (int i = 0; i < m; i++) {
             if (!result->edges[i].active) continue;   /* 已分類，跳過 */
 
             int    u = result->edges[i].u;
             int    v = result->edges[i].v;
 
             /* 計算第 k 層共同鄰居比例（論文公式 (1) 或 (2)）*/
             double ratio = compute_cnr(g, u, v, k, buf,
                                        bs_k_u, bs_k_v,
                                        bs_k1_u, bs_k1_v);
 
             if (ratio >= T_E) {
                 /* ── bond link：R^k 超過外部閾值 ── */
                 result->edges[i].type   = LINK_BOND;
                 result->edges[i].active = 0;   /* 標記已分類 */
                 n_bond++;
             } else {
                 /* ── 候選邊：加入候選列表，等待內部閾值判斷 ── */
                 cand_ratios[n_cand] = ratio;
                 cand_idx[n_cand]    = i;
                 n_cand++;
             }
         }
 
         fprintf(stdout, "[heta] Layer %d: T_E=%.4f, bonds=%d, candidates=%d\n",
                 k, T_E, n_bond, n_cand);
 
         /* ── 若有候選邊，計算內部閾值並識別 local bridges ── */
         if (n_cand > 0) {
             /* 論文公式 (6)：T_I^k = mean(candidates) - SD(candidates) */
             double T_I = compute_internal_threshold(cand_ratios, n_cand);
 
             int n_local = 0;
             for (int ci = 0; ci < n_cand; ci++) {
                 int idx = cand_idx[ci];
 
                 if (cand_ratios[ci] > T_I) {
                     /* ── kth-layer local bridge：R^k > T_I^k ── */
                     result->edges[idx].type   = LINK_LOCAL_BRIDGE;
                     result->edges[idx].layer  = k;   /* 記錄所在層號 */
                     result->edges[idx].active = 0;   /* 標記已分類   */
                     n_local++;
                 }
                 /* 否則：active 保持 1，繼續傳到下一層 */
             }
 
             fprintf(stdout, "[heta] Layer %d: T_I=%.4f, local_bridges=%d\n",
                     k, T_I, n_local);
         }
     }
 
     /* ── Step 3.3：識別 global bridges ── */
     /* 所有層都未能分類（active 仍為 1）的邊 = global bridge */
     int n_global = 0;
     for (int i = 0; i < m; i++) {
         if (result->edges[i].active) {
             result->edges[i].type   = LINK_GLOBAL_BRIDGE;
             result->edges[i].active = 0;
             n_global++;
         }
     }
     fprintf(stdout, "[heta] Step 3.3: %d global bridges identified\n", n_global);
 
     /* ── 清理工作記憶體 ── */
     bfs_buffer_free(buf);
     bs_free(bs_k_u);
     bs_free(bs_k_v);
     bs_free(bs_k1_u);
     bs_free(bs_k1_v);
     free(cand_ratios);
     free(cand_idx);
 
     return result;
 }
 
 /* ════════════════════════════════════════════════════════════
  * heta_summarize — 統計各類型數量與比例
  * ════════════════════════════════════════════════════════════ */
 HetaSummary heta_summarize(const HetaResult *r) {
     HetaSummary s;
     s.n_silk   = 0;
     s.n_bond   = 0;
     s.n_local  = 0;
     s.n_global = 0;
 
     /* 逐邊統計類型 */
     for (int i = 0; i < r->m; i++) {
         switch (r->edges[i].type) {
             case LINK_SILK:          s.n_silk++;   break;
             case LINK_BOND:          s.n_bond++;   break;
             case LINK_LOCAL_BRIDGE:  s.n_local++;  break;
             case LINK_GLOBAL_BRIDGE: s.n_global++; break;
             default: break;   /* LINK_UNCLASSIFIED 不應出現在最終結果 */
         }
     }
 
     /* 計算百分比 */
     double total = (double)r->m;
     s.pct_silk   = (total > 0) ? (s.n_silk   / total * 100.0) : 0.0;
     s.pct_bond   = (total > 0) ? (s.n_bond   / total * 100.0) : 0.0;
     s.pct_local  = (total > 0) ? (s.n_local  / total * 100.0) : 0.0;
     s.pct_global = (total > 0) ? (s.n_global / total * 100.0) : 0.0;
 
     return s;
 }
 
 /* ════════════════════════════════════════════════════════════
  * heta_print_summary — 列印統計摘要到 stdout
  * ════════════════════════════════════════════════════════════ */
 void heta_print_summary(const HetaSummary *s, int total_edges) {
     fprintf(stdout, "\n════════════════════════════════\n");
     fprintf(stdout, " HETA 分類結果摘要（共 %d 條邊）\n", total_edges);
     fprintf(stdout, "════════════════════════════════\n");
     fprintf(stdout, "  %-14s : %4d 條  (%5.1f%%)\n",
             "silk link",     s->n_silk,   s->pct_silk);
     fprintf(stdout, "  %-14s : %4d 條  (%5.1f%%)\n",
             "bond link",     s->n_bond,   s->pct_bond);
     fprintf(stdout, "  %-14s : %4d 條  (%5.1f%%)\n",
             "local bridge",  s->n_local,  s->pct_local);
     fprintf(stdout, "  %-14s : %4d 條  (%5.1f%%)\n",
             "global bridge", s->n_global, s->pct_global);
     fprintf(stdout, "════════════════════════════════\n\n");
 }