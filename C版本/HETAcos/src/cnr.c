/*
 * cnr.c — 第 k 層鄰居 BFS 與共同鄰居比例（CNR）計算
 * =====================================================
 *
 * 對應論文公式：
 *   公式 (1)：k=1 時的共同鄰居比例 R^1_{u,v}
 *   公式 (2)：k>1 時的共同鄰居比例 R^k_{u,v}
 *   公式 (3)：第 k 層鄰居集合 V^{k,j}_i 的遞迴定義
 *
 * 關鍵優化（相對 Python 版本）：
 *   1. BFS 使用靜態陣列 queue + visited（無動態 set 分配）
 *   2. 鄰居集合用 Bitset 表示，交集 = bitwise AND
 *   3. BFSBuffer 由呼叫者預分配，thread-local（支援 OpenMP 並行）
 *   4. k>1 時，Vk-1 在計算 Vk 過程中已知，無需重算（快取在 buf 中）
 */

#include "heta.h"

#include <stdlib.h>   /* malloc, free   */
#include <string.h>   /* memset         */
#include <stdio.h>    /* fprintf        */
#include <math.h>     /* sqrt           */

/* 安全 malloc */
static void *xmalloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "[cnr] malloc failed\n"); exit(1); }
    return p;
}

/* ════════════════════════════════════════════════════════════
 * bfs_buffer_alloc — 分配 thread-local BFS 工作緩衝區
 *
 * 每個 OpenMP thread 呼叫一次，之後反覆重用（無需 malloc/free）。
 * 三個陣列大小皆為 n（節點數），確保最大 BFS 深度下不溢位。
 * ════════════════════════════════════════════════════════════ */
BFSBuffer *bfs_buffer_alloc(int n) {
    BFSBuffer *buf  = (BFSBuffer *)xmalloc(sizeof(BFSBuffer));
    buf->queue      = (int *)xmalloc(n * sizeof(int)); /* BFS 隊列         */
    buf->visited    = (int *)xmalloc(n * sizeof(int)); /* 0/1 訪問旗標     */
    buf->layer_buf  = (int *)xmalloc(n * sizeof(int)); /* 節點所屬層號     */
    return buf;
}

/* ════════════════════════════════════════════════════════════
 * bfs_buffer_free — 釋放 BFS 緩衝區
 * ════════════════════════════════════════════════════════════ */
void bfs_buffer_free(BFSBuffer *buf) {
    if (!buf) return;
    free(buf->queue);
    free(buf->visited);
    free(buf->layer_buf);
    free(buf);
}

/* ════════════════════════════════════════════════════════════
 * compute_kth_layer — 計算節點 u（排除 exclude）的第 k 層鄰居
 *
 * 對應論文公式 (3)：
 *   V^{k,j}_i = { v_y | ∀v_x ∈ V^{k-1,j}_i : (v_x, v_y) ∈ E }
 *               \ ( ⋃_{p=1}^{k-1} V^{p,j}_i ) \ {v_i, v_j}
 *
 * 實作方式：
 *   - 做一次 BFS，記錄每個節點第一次被訪問時的層號
 *   - 第 k 層節點 = 層號恰好等於 k 的節點
 *   - 使用 visited 陣列（而非 set）避免重複訪問
 *   - 結果寫入 Bitset out
 *
 * 參數：
 *   g       : CSR 格式圖
 *   u       : 起始節點
 *   exclude : 邊的另一端點（論文要求排除此節點及其前置路徑影響）
 *   k       : 目標層數（1-indexed）
 *   out     : 輸出 bitset（呼叫前必須已 bs_clear）
 *   buf     : thread-local BFS 緩衝區
 * ════════════════════════════════════════════════════════════ */
void compute_kth_layer(const Graph *g, int u, int exclude, int k,
                       Bitset *out, BFSBuffer *buf) {
    int n = g->n;

    /* ── 初始化工作陣列 ── */
    /* 用 memset 清零 visited 和 layer_buf（比 for 迴圈快）*/
    memset(buf->visited,   0, n * sizeof(int));
    memset(buf->layer_buf, -1, n * sizeof(int)); /* -1 代表「未訪問」*/

    /* 清空輸出 bitset */
    bs_clear(out);

    /* ── BFS 初始化：起點 u，排除 exclude ── */
    buf->visited[u]       = 1;   /* 標記 u 已訪問（第 0 層）*/
    buf->visited[exclude] = 1;   /* 標記 exclude 為不可訪問（論文要求排除） */
    buf->layer_buf[u]     = 0;   /* u 本身在第 0 層 */

    int head = 0, tail = 0;
    buf->queue[tail++] = u;      /* 將起點放入隊列 */

    /* ── BFS 主迴圈：展開到第 k 層即停止 ── */
    while (head < tail) {
        int cur       = buf->queue[head++];   /* 取出隊首                */
        int cur_layer = buf->layer_buf[cur];  /* 此節點的層號            */

        /* 若已到達第 k 層，不再往下展開（第 k 層的鄰居在 k+1 層，不需要）*/
        if (cur_layer >= k) continue;

        /* 遍歷 cur 的所有鄰居（CSR 連續記憶體，cache 友善）*/
        for (int i = g->offset[cur]; i < g->offset[cur + 1]; i++) {
            int nb = g->adj_list[i];   /* 鄰居節點 id */

            if (!buf->visited[nb]) {                  /* 尚未訪問 */
                buf->visited[nb]   = 1;               /* 標記已訪問 */
                buf->layer_buf[nb] = cur_layer + 1;   /* 記錄層號   */
                buf->queue[tail++] = nb;              /* 加入隊列   */

                /* 若此鄰居恰好在第 k 層，加入輸出 bitset */
                if (cur_layer + 1 == k) {
                    bs_set(out, nb);
                }
            }
        }
    }
}

/* ════════════════════════════════════════════════════════════
 * compute_cnr — 計算邊 (u,v) 的第 k 層共同鄰居比例 R^k_{u,v}
 *
 * k=1 對應論文公式 (1)（本實作改用幾何平均分母）：
 *   R^1_{u,v} = |V^1_u ∩ V^1_v| / sqrt(|V^1_u| * |V^1_v|)
 *
 * k>1 對應論文公式 (2)（本實作改用幾何平均分母）：
 *   分子 = |Vk_u ∩ Vk_v| + |Vk1_u ∩ Vk_v| + |Vk_u ∩ Vk1_v|
 *   分母 = sqrt(|Vk_u|*|Vk_v|) + sqrt(|Vk1_u|*|Vk_v|) + sqrt(|Vk_u|*|Vk1_v|)
 *   其中 Vk1 = V^{k-1}
 *
 * 註：原始論文公式採用 min(.,.) 作為分母；此版本將 min 改為
 *     sqrt(|A|*|B|)（幾何平均，等價於 Salton/cosine 相似度的歸一化方式），
 *     使指標對兩端集合大小的敏感度更平衡。
 *
 * 若 union 為空（三組交集均為空），回傳 0.0
 *
 * 預分配的 bitset 工作空間由呼叫者提供，避免反覆 malloc：
 *   bs_k_u, bs_k_v   : 第 k 層鄰居
 *   bs_k1_u, bs_k1_v : 第 k-1 層鄰居（k>1 時使用）
 * ════════════════════════════════════════════════════════════ */
double compute_cnr(const Graph *g, int u, int v, int k,
                   BFSBuffer *buf,
                   Bitset *bs_k_u,  Bitset *bs_k_v,
                   Bitset *bs_k1_u, Bitset *bs_k1_v) {

    if (k == 1) {
        /* ── k=1：論文公式 (1)（min → sqrt 版）── */

        /* 計算 V^1_u（排除 v）和 V^1_v（排除 u）*/
        compute_kth_layer(g, u, v, 1, bs_k_u, buf);
        compute_kth_layer(g, v, u, 1, bs_k_v, buf);

        int size_u = bs_popcount(bs_k_u);   /* |V^1_u| */
        int size_v = bs_popcount(bs_k_v);   /* |V^1_v| */

        /* 若任一端點無鄰居（除了對方），比例為 0 */
        if (size_u == 0 || size_v == 0) return 0.0;

        /* 計算交集大小 |V^1_u ∩ V^1_v| */
        int inter = bs_intersect_count(bs_k_u, bs_k_v);

        /* 分母 = sqrt(|V^1_u| * |V^1_v|)（幾何平均，Salton/cosine 形式）*/
        double denom = sqrt((double)size_u * (double)size_v);

        /* size_u、size_v 均 > 0，denom 必 > 0；防禦性檢查仍保留 */
        if (denom <= 0.0) return 0.0;

        return (double)inter / denom;

    } else {
        /* ── k>1：論文公式 (2)（min → sqrt 版）── */

        /* 計算第 k 層鄰居 */
        compute_kth_layer(g, u, v, k,     bs_k_u,  buf);
        compute_kth_layer(g, v, u, k,     bs_k_v,  buf);

        /* 計算第 k-1 層鄰居 */
        compute_kth_layer(g, u, v, k - 1, bs_k1_u, buf);
        compute_kth_layer(g, v, u, k - 1, bs_k1_v, buf);

        int sz_k_u  = bs_popcount(bs_k_u);    /* |Vk_u|  */
        int sz_k_v  = bs_popcount(bs_k_v);    /* |Vk_v|  */
        int sz_k1_u = bs_popcount(bs_k1_u);   /* |Vk1_u| */
        int sz_k1_v = bs_popcount(bs_k1_v);   /* |Vk1_v| */

        /* ── 計算三組交集大小（論文公式 (2) 分子各項）── */

        /* 項 1：|Vk_u ∩ Vk_v|   — 兩邊都在第 k 層的共同節點 */
        int inter_kk   = bs_intersect_count(bs_k_u,  bs_k_v);

        /* 項 2：|Vk1_u ∩ Vk_v|  — u 的 k-1 層 與 v 的 k 層 的共同節點 */
        int inter_k1k  = bs_intersect_count(bs_k1_u, bs_k_v);

        /* 項 3：|Vk_u ∩ Vk1_v|  — u 的 k 層 與 v 的 k-1 層 的共同節點 */
        int inter_kk1  = bs_intersect_count(bs_k_u,  bs_k1_v);

        /* 若三組交集的聯集為空（無任何共同鄰居），比例為 0 */
        if (inter_kk == 0 && inter_k1k == 0 && inter_kk1 == 0) return 0.0;

        /* ── 計算分子（論文公式 (2) numerator）── */
        int num = inter_kk + inter_k1k + inter_kk1;

        /* ── 計算分母（公式 (2) denominator，min → sqrt 版）── */
        /* 每項取兩個集合大小的幾何平均 sqrt(|A|*|B|)（Salton/cosine 形式）*/
        double sq_kk  = sqrt((double)sz_k_u  * (double)sz_k_v);   /* sqrt(|Vk_u|*|Vk_v|)   */
        double sq_k1k = sqrt((double)sz_k1_u * (double)sz_k_v);   /* sqrt(|Vk1_u|*|Vk_v|)  */
        double sq_kk1 = sqrt((double)sz_k_u  * (double)sz_k1_v);  /* sqrt(|Vk_u|*|Vk1_v|)  */

        double denom = sq_kk + sq_k1k + sq_kk1;

        /* 防止除以 0（任一集合為空會讓對應項變 0；三項皆 0 時整體為 0）*/
        if (denom <= 0.0) return 0.0;

        return (double)num / denom;
    }
}
