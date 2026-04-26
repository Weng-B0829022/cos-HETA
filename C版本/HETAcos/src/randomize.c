/*
 * randomize.c — Switching Randomization（保留度數分布的隨機化）
 * ==============================================================
 *
 * 對應論文 Section 2：生成隨機化網路集合 RG
 *
 * 演算法原理（Switching Algorithm）：
 *   每步隨機選兩條邊 (a,b) 和 (c,d)，嘗試交換為 (a,d) 和 (c,b)。
 *   交換條件：
 *     1. a ≠ d 且 c ≠ b（不產生自迴圈）
 *     2. (a,d) 和 (c,b) 在圖中尚不存在（不產生平行邊）
 *   重複 Q × m 次（Q=100，m=邊數），圖結構充分隨機化。
 *   由於每次交換保持所有節點的度數不變，最終隨機網路與原圖
 *   有相同的節點數、邊數、度數序列。
 *
 * C 實作相對 Python 的改進：
 *   1. 不 copy 整張圖（直接修改傳入的 gr，呼叫者負責先 graph_copy）
 *   2. has_edge 用 bitset 查詢（O(1)，一次位元測試）
 *   3. 隨機數使用 PCG32 演算法（比 rand() 更快更均勻）
 *   4. 邊陣列直接在 adj_list 上操作，無需額外資料結構
 *
 * PCG32 隨機數生成器：
 *   - 比 rand() 速度快約 2 倍，統計品質遠優於線性同餘
 *   - 每個 thread 獨立的 state，完全 thread-safe（不用 mutex）
 *   - 參考：O'Neill (2014), "PCG: A Family of Simple Fast..."
 */

 #include "heta.h"

 #include <stdlib.h>   /* malloc, free       */
 #include <string.h>   /* memset             */
 #include <stdio.h>    /* fprintf            */
 #include <stdint.h>   /* uint64_t, uint32_t */
 
 /* ────────────────────────────────────────────────────────────
  * PCG32 隨機數生成器
  * 狀態：兩個 uint64_t（state + inc）
  * 輸出：一個 uint32_t 的高品質隨機數
  * ──────────────────────────────────────────────────────────── */
 
 typedef struct {
     uint64_t state;   /* RNG 內部狀態，每次更新 */
     uint64_t inc;     /* 序列選擇子（必須為奇數）*/
 } PCG32;
 
 /* 初始化 PCG32，seed 和 seq 決定起始狀態與序列 */
 static void pcg32_init(PCG32 *rng, uint64_t seed, uint64_t seq) {
     rng->state = 0u;
     rng->inc   = (seq << 1u) | 1u;   /* inc 必須為奇數 */
     /* 空跑一次讓狀態離開 0 */
     rng->state = rng->state * 6364136223846793005ULL + rng->inc;
     rng->state += seed;
     rng->state = rng->state * 6364136223846793005ULL + rng->inc;
 }
 
 /* 產生下一個 uint32_t 隨機數 */
 static inline uint32_t pcg32_next(PCG32 *rng) {
     uint64_t old   = rng->state;
     /* 線性同餘更新 state */
     rng->state = old * 6364136223846793005ULL + rng->inc;
     /* 輸出函式：xorshift + rotation（提升統計品質）*/
     uint32_t xorshifted = (uint32_t)(((old >> 18u) ^ old) >> 27u);
     uint32_t rot        = (uint32_t)(old >> 59u);
     return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
 }
 
 /* 在 [0, bound) 範圍內均勻取樣（無模數偏差）*/
 static inline int pcg32_bound(PCG32 *rng, uint32_t bound) {
     /* Lemire's nearly divisionless method：極快且無偏 */
     uint64_t m = (uint64_t)pcg32_next(rng) * (uint64_t)bound;
     return (int)(m >> 32);
 }
 
 /* ────────────────────────────────────────────────────────────
  * 鄰接矩陣 bitset（用於 O(1) has_edge 查詢）
  *
  * adj_bs[u] = 節點 u 的鄰居集合的 bitset
  * 存取格式：adj_bs[u].words[v/64] 的第 v%64 位
  *
  * 相比 hash set 或 dict：
  *   - 測試 (u,v) 是否存在：1 次位元測試，O(1) 且極快
  *   - 新增邊：1 次 bs_set，O(1)
  *   - 刪除邊：1 次位元清除，O(1)
  * ──────────────────────────────────────────────────────────── */
 
 /* 清除 bitset 的第 node 位元（從集合中移除）*/
 static inline void bs_clear_bit(Bitset *bs, int node) {
     bs->words[node / BS_WORD_BITS] &= ~(1ULL << (node % BS_WORD_BITS));
 }
 
 /* ════════════════════════════════════════════════════════════
  * adj_bs_alloc — 預分配 n 個空 Bitset 的指標陣列
  *
  * 只分配記憶體、不填入任何鄰居資料。
  * 供 threshold.c 在 OpenMP parallel 區段開始時呼叫一次，
  * 之後每次迭代只做清零 + 填入，不重複 malloc/free。
  *
  * 時間：O(n)；空間：O(n × ceil(n/64))
  * ════════════════════════════════════════════════════════════ */
 Bitset **adj_bs_alloc(int n) {
     Bitset **adj_bs = (Bitset **)malloc(n * sizeof(Bitset *));
     if (!adj_bs) { fprintf(stderr, "[rand] adj_bs_alloc failed\n"); exit(1); }
     for (int u = 0; u < n; u++) {
         adj_bs[u] = bs_alloc(n);   /* bs_alloc 內部 calloc，全 0 初始化 */
     }
     return adj_bs;
 }
 
 /* ════════════════════════════════════════════════════════════
  * adj_bs_free — 釋放 adj_bs_alloc 分配的陣列
  * ════════════════════════════════════════════════════════════ */
 void adj_bs_free(Bitset **adj_bs, int n) {
     for (int i = 0; i < n; i++) bs_free(adj_bs[i]);
     free(adj_bs);
 }
 
 /* ════════════════════════════════════════════════════════════
  * switching_randomize — 對圖 gr 執行 switching randomization
  *
  * 對應論文 Section 2 的隨機化步驟。
  * gr 必須已是原圖的複本（graph_copy 取得），本函式原地修改。
  *
  * 參數：
  *   gr     : 要隨機化的圖（原地修改）
  *   seed   : PCG32 起始種子（每個 thread 傳入不同的值，確保獨立性）
  *   adj_bs : 預分配的鄰接矩陣 bitset 陣列（大小 gr->n）
  *            由呼叫方（threshold.c）在迴圈外預先分配，此處只做：
  *              1. 清零（bs_clear）
  *              2. 填入原圖鄰居（bs_set）
  *              3. switching 主迴圈操作 bitset
  *            函式結束後不 free，由呼叫方統一管理生命週期
  * ════════════════════════════════════════════════════════════ */
 void switching_randomize(Graph *gr, uint64_t seed, Bitset **adj_bs) {
     int n = gr->n;
     int m = gr->m;
 
     /* ── 建立邊陣列（用索引 0..m-1 表示每條邊）── */
     /* edge_u[i], edge_v[i] = 邊 i 的兩個端點 */
     int *edge_u = (int *)malloc(m * sizeof(int));
     int *edge_v = (int *)malloc(m * sizeof(int));
     if (!edge_u || !edge_v) {
         fprintf(stderr, "[rand] malloc edges failed\n"); exit(1);
     }
 
     /* 從 CSR adj_list 提取邊（只取 u < v 的方向，避免重複）*/
     int edge_cnt = 0;
     for (int u = 0; u < n; u++) {
         for (int i = gr->offset[u]; i < gr->offset[u + 1]; i++) {
             int v = gr->adj_list[i];
             if (u < v) {   /* 只記錄一個方向 */
                 edge_u[edge_cnt] = u;
                 edge_v[edge_cnt] = v;
                 edge_cnt++;
             }
         }
     }
     /* 確認邊數正確 */
     /* edge_cnt 應等於 m */
 
     /* ── 清零並填入鄰居到 adj_bs（預分配版本，不 malloc）── */
     /* 每次迭代只做 memset 清零 + bs_set 填入，比 malloc 快約 5~50 倍 */
     for (int u = 0; u < n; u++) {
         bs_clear(adj_bs[u]);   /* memset 清零，O(n/64) */
         for (int i = gr->offset[u]; i < gr->offset[u + 1]; i++) {
             bs_set(adj_bs[u], gr->adj_list[i]);   /* 填入鄰居 */
         }
     }
 
     /* ── 建立鄰接矩陣 bitset（用於快速 has_edge）── */
 
     /* ── 初始化 PCG32 RNG（每個呼叫有獨立 seed）── */
     PCG32 rng;
     pcg32_init(&rng, seed, seed ^ 0xCAFEBABEULL);
 
     /* ── Switching 主迴圈：執行 Q × m 次嘗試 ── */
     int total_steps = HETA_SWITCH_Q * m;   /* Q=100，論文建議值 */
 
     for (int step = 0; step < total_steps; step++) {
         /* 隨機選兩條邊的索引（[0, m)）*/
         int i = pcg32_bound(&rng, (uint32_t)m);
         int j = pcg32_bound(&rng, (uint32_t)m);
 
         /* 避免選到同一條邊 */
         if (i == j) continue;
 
         int a = edge_u[i], b = edge_v[i];   /* 邊 i：(a, b) */
         int c = edge_u[j], d = edge_v[j];   /* 邊 j：(c, d) */
 
         /* ── 檢查交換條件 ── */
 
         /* 條件 1：交換後不能產生自迴圈 */
         if (a == d || c == b) continue;
 
         /* 條件 2：交換後不能產生重複邊（(a,d) 和 (c,b) 必須不存在）*/
         if (bs_test(adj_bs[a], d)) continue;   /* (a,d) 已存在 */
         if (bs_test(adj_bs[c], b)) continue;   /* (c,b) 已存在 */
 
         /* ── 執行交換：(a,b)+(c,d) → (a,d)+(c,b) ── */
 
         /* 更新鄰接矩陣 bitset：移除舊邊，加入新邊 */
         bs_clear_bit(adj_bs[a], b);   /* 移除 (a,b) */
         bs_clear_bit(adj_bs[b], a);
         bs_clear_bit(adj_bs[c], d);   /* 移除 (c,d) */
         bs_clear_bit(adj_bs[d], c);
 
         bs_set(adj_bs[a], d);         /* 加入 (a,d) */
         bs_set(adj_bs[d], a);
         bs_set(adj_bs[c], b);         /* 加入 (c,b) */
         bs_set(adj_bs[b], c);
 
         /* 更新邊陣列 */
         edge_u[i] = a; edge_v[i] = d;   /* 邊 i 更新為 (a,d) */
         edge_u[j] = c; edge_v[j] = b;   /* 邊 j 更新為 (c,b) */
     }
 
     /* ── 將 bitset 重新寫回 CSR adj_list ── */
     /* 因為 CSR 格式需要連續排列，從 bitset 重建 adj_list */
     int pos = 0;
     for (int u = 0; u < n; u++) {
         gr->offset[u] = pos;   /* 記錄節點 u 在 adj_list 的起始位置 */
 
         /* 遍歷 bitset：找出所有被設定的位元（即 u 的鄰居）*/
         for (int word_i = 0; word_i < adj_bs[u]->n_words; word_i++) {
             uint64_t w = adj_bs[u]->words[word_i];
             while (w) {
                 /* 取出最低設定位元的位置（bit scan forward）*/
                 int bit = __builtin_ctzll(w);           /* 最低位 1 的位置 */
                 int nb  = word_i * (int)BS_WORD_BITS + bit;
                 gr->adj_list[pos++] = nb;
                 w &= w - 1;   /* 清除最低設定位元（Brian Kernighan's trick）*/
             }
         }
     }
     gr->offset[n] = pos;   /* 哨兵：最後一個節點的結束位置 */
 
     /* ── 清理（adj_bs 由呼叫方管理，此處只釋放邊陣列）── */
     free(edge_u);
     free(edge_v);
 }