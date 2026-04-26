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
 * k=1 對應論文公式 (1)：
 *   R^1_{u,v} = |V^1_u ∩ V^1_v| / min(|V^1_u|, |V^1_v|)
 *
 * k>1 對應論文公式 (2)：
 *   分子 = |Vk_u ∩ Vk_v| + |Vk1_u ∩ Vk_v| + |Vk_u ∩ Vk1_v|
 *   分母 = min(|Vk_u|,|Vk_v|) + min(|Vk1_u|,|Vk_v|) + min(|Vk_u|,|Vk1_v|)
 *   其中 Vk1 = V^{k-1}
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
        /* ── k=1：論文公式 (1) ── */

        /* 計算 V^1_u（排除 v）和 V^1_v（排除 u）*/
        compute_kth_layer(g, u, v, 1, bs_k_u, buf);
        compute_kth_layer(g, v, u, 1, bs_k_v, buf);

        int size_u = bs_popcount(bs_k_u);   /* |V^1_u| */
        int size_v = bs_popcount(bs_k_v);   /* |V^1_v| */

        /* 若任一端點無鄰居（除了對方），比例為 0 */
        if (size_u == 0 || size_v == 0) return 0.0;

        /* 計算交集大小 |V^1_u ∩ V^1_v| */
        int inter = bs_intersect_count(bs_k_u, bs_k_v);

        /* 分母 = min(|V^1_u|, |V^1_v|)（不可超過較小集合的大小）*/
        int denom = (size_u < size_v) ? size_u : size_v;

        return (double)inter / (double)denom;

    } else {
        /* ── k>1：論文公式 (2) ── */

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

        /* ── 計算分母（論文公式 (2) denominator）── */
        /* 每項取兩個集合大小的 min，確保比例不超過 1 */
        int min_kk   = (sz_k_u  < sz_k_v)  ? sz_k_u  : sz_k_v;   /* min(|Vk_u|, |Vk_v|)   */
        int min_k1k  = (sz_k1_u < sz_k_v)  ? sz_k1_u : sz_k_v;   /* min(|Vk1_u|, |Vk_v|)  */
        int min_kk1  = (sz_k_u  < sz_k1_v) ? sz_k_u  : sz_k1_v;  /* min(|Vk_u|, |Vk1_v|)  */

        int denom = min_kk + min_k1k + min_kk1;

        /* 防止除以 0（理論上 union 非空時 denom > 0，但防禦性處理）*/
        if (denom == 0) return 0.0;

        return (double)num / (double)denom;
    }
}

/* ════════════════════════════════════════════════════════════
 * Ego Table — 對齊 Python HETA 的 generate_ego_graph()
 *
 * 預先計算所有節點的 0~(layers-1) 階「不受限」鄰域 bitset：
 *   ego(u, 0) = {u}
 *   ego(u, r) = ego(u, r-1) ∪ ⋃_{ng ∈ N(u)} ego(ng, r-1)
 *
 * 注意：「不受限」意指這個 ego 不排除任何節點。當邊 (s,t) 之分析需要
 * 排除「直接 s-t 邊」時，是在 compute_all_link_ratios 內處理（透過
 * 「只 union s 的非 t 鄰居的 ego」與「結果扣除 {s, t}」），而不是改
 * 動 ego table 本身。這正好對應 Python HETA 的設計。
 * ════════════════════════════════════════════════════════════ */
EgoTable *ego_table_alloc(const Graph *g, int layers) {
    if (layers < 1) {
        fprintf(stderr, "[cnr] ego_table_alloc: layers must be >= 1\n");
        exit(1);
    }
    int n = g->n;

    EgoTable *et = (EgoTable *)xmalloc(sizeof(EgoTable));
    et->n        = n;
    et->n_layers = layers;
    et->bitsets  = (Bitset **)xmalloc((size_t)n * (size_t)layers * sizeof(Bitset *));
    for (int u = 0; u < n; u++) {
        for (int r = 0; r < layers; r++) {
            et->bitsets[(size_t)u * layers + r] = bs_alloc(n);
        }
    }

    /* ── Layer 0：每個節點的 0-hop ego = {自己}  ── */
    for (int u = 0; u < n; u++) {
        bs_set(et->bitsets[(size_t)u * layers + 0], u);
    }

    /* ── Layer r ≥ 1：ego(u, r) = ego(u, r-1) ∪ ⋃_{ng ∈ N(u)} ego(ng, r-1) ── */
    for (int r = 1; r < layers; r++) {
        for (int u = 0; u < n; u++) {
            Bitset *cur = et->bitsets[(size_t)u * layers + r];
            /* 起點：複製 ego(u, r-1) */
            bs_copy(cur, et->bitsets[(size_t)u * layers + r - 1]);
            /* 對每個鄰居 ng，把 ego(ng, r-1) 聯集進來 */
            for (int i = g->offset[u]; i < g->offset[u + 1]; i++) {
                int ng = g->adj_list[i];
                bs_or_into(cur, et->bitsets[(size_t)ng * layers + r - 1]);
            }
        }
    }

    return et;
}

void ego_table_free(EgoTable *et) {
    if (!et) return;
    if (et->bitsets) {
        for (int u = 0; u < et->n; u++) {
            for (int r = 0; r < et->n_layers; r++) {
                bs_free(et->bitsets[(size_t)u * et->n_layers + r]);
            }
        }
        free(et->bitsets);
    }
    free(et);
}

/* ────────────────────────────────────────────────────────────
 * RingWorkspace — 每邊計算所需的 thread-local bitset 工作空間
 *
 * 因為每條邊的 ring 計算彼此獨立，所以 workspace 可在迴圈外只分配一次。
 * 在 OpenMP 並行裡每個 thread 各持一份。
 * ──────────────────────────────────────────────────────────── */
typedef struct {
    int      n;
    int      layers;
    Bitset **ring_s;     /* ring_s[0..layers]，ring_s[0] 永遠空（佔位）   */
    Bitset **ring_t;
    Bitset  *acc_s;      /* 累積已涵蓋節點（每邊重置）                    */
    Bitset  *acc_t;
    Bitset  *get_ego_s;  /* get_ego_graph(s, t, l) 暫存                   */
    Bitset  *get_ego_t;
    Bitset  *common;     /* 三組交集的聯集                                */
    Bitset  *temp;       /* 計算交集時的暫存                              */
} RingWorkspace;

static RingWorkspace *ring_workspace_alloc(int n, int layers) {
    RingWorkspace *w = (RingWorkspace *)xmalloc(sizeof(RingWorkspace));
    w->n      = n;
    w->layers = layers;
    w->ring_s = (Bitset **)xmalloc((size_t)(layers + 1) * sizeof(Bitset *));
    w->ring_t = (Bitset **)xmalloc((size_t)(layers + 1) * sizeof(Bitset *));
    for (int l = 0; l <= layers; l++) {
        w->ring_s[l] = bs_alloc(n);
        w->ring_t[l] = bs_alloc(n);
    }
    w->acc_s     = bs_alloc(n);
    w->acc_t     = bs_alloc(n);
    w->get_ego_s = bs_alloc(n);
    w->get_ego_t = bs_alloc(n);
    w->common    = bs_alloc(n);
    w->temp      = bs_alloc(n);
    return w;
}

static void ring_workspace_free(RingWorkspace *w) {
    if (!w) return;
    if (w->ring_s) {
        for (int l = 0; l <= w->layers; l++) bs_free(w->ring_s[l]);
        free(w->ring_s);
    }
    if (w->ring_t) {
        for (int l = 0; l <= w->layers; l++) bs_free(w->ring_t[l]);
        free(w->ring_t);
    }
    bs_free(w->acc_s);
    bs_free(w->acc_t);
    bs_free(w->get_ego_s);
    bs_free(w->get_ego_t);
    bs_free(w->common);
    bs_free(w->temp);
    free(w);
}

/* ────────────────────────────────────────────────────────────
 * 對單一邊 (s, t) 計算 layer 1..layers 的 ratio
 *
 * 邏輯完全對應 Python HETA 的 compute_link_property（見論文 / heta.core.py）。
 *
 * 寫入：out_ratios[l-1] = 邊在層 l 的 ratio，l = 1..layers
 * ──────────────────────────────────────────────────────────── */
static void compute_edge_ratios(const Graph *g, int s, int t, int layers,
                                const EgoTable *et, RingWorkspace *w,
                                double *out_ratios) {
    /* 重設 accumulator 與 ring[0] 為空（每邊起點皆為空集合）*/
    bs_clear(w->acc_s);
    bs_clear(w->acc_t);
    bs_clear(w->ring_s[0]);   /* l-1=0 時的「ring」就是空（accumulator 起始狀態）*/
    bs_clear(w->ring_t[0]);

    for (int l = 1; l <= layers; l++) {
        /* ── Step 1：get_ego_graph(s, t, l)
         *   = ⋃_{ng ∈ N(s), ng ≠ t} ego(ng, l-1)  −  {s}
         */
        bs_clear(w->get_ego_s);
        for (int i = g->offset[s]; i < g->offset[s + 1]; i++) {
            int ng = g->adj_list[i];
            if (ng == t) continue;
            bs_or_into(w->get_ego_s,
                       et->bitsets[(size_t)ng * et->n_layers + (l - 1)]);
        }
        bs_clear_bit(w->get_ego_s, s);

        bs_clear(w->get_ego_t);
        for (int i = g->offset[t]; i < g->offset[t + 1]; i++) {
            int ng = g->adj_list[i];
            if (ng == s) continue;
            bs_or_into(w->get_ego_t,
                       et->bitsets[(size_t)ng * et->n_layers + (l - 1)]);
        }
        bs_clear_bit(w->get_ego_t, t);

        /* ── Step 2：ring_s[l] = get_ego_s − acc_s − {s, t} ── */
        bs_copy(w->ring_s[l], w->get_ego_s);
        bs_subtract_into(w->ring_s[l], w->acc_s);
        bs_clear_bit(w->ring_s[l], s);
        bs_clear_bit(w->ring_s[l], t);

        bs_copy(w->ring_t[l], w->get_ego_t);
        bs_subtract_into(w->ring_t[l], w->acc_t);
        bs_clear_bit(w->ring_t[l], s);
        bs_clear_bit(w->ring_t[l], t);

        /* ── Step 3：common = (s_l ∩ t_l) ∪ (s_l ∩ t_{l-1}) ∪ (s_{l-1} ∩ t_l)
         *
         * 不能直接用三個 popcount 加總（會把同一節點重複算），
         * 必須先聯集再 popcount。
         */
        bs_clear(w->common);
        bs_intersect_into(w->ring_s[l],     w->ring_t[l],     w->temp);
        bs_or_into(w->common, w->temp);
        bs_intersect_into(w->ring_s[l],     w->ring_t[l - 1], w->temp);
        bs_or_into(w->common, w->temp);
        bs_intersect_into(w->ring_s[l - 1], w->ring_t[l],     w->temp);
        bs_or_into(w->common, w->temp);

        int common_count = bs_popcount(w->common);

        if (common_count == 0) {
            out_ratios[l - 1] = 0.0;
        } else {
            int sz_s_l   = bs_popcount(w->ring_s[l]);
            int sz_t_l   = bs_popcount(w->ring_t[l]);
            int sz_s_lm1 = bs_popcount(w->ring_s[l - 1]);
            int sz_t_lm1 = bs_popcount(w->ring_t[l - 1]);

            int min_a = (sz_s_l   < sz_t_l)   ? sz_s_l   : sz_t_l;
            int min_b = (sz_s_l   < sz_t_lm1) ? sz_s_l   : sz_t_lm1;
            int min_c = (sz_s_lm1 < sz_t_l)   ? sz_s_lm1 : sz_t_l;
            int denom = min_a + min_b + min_c;

            out_ratios[l - 1] = (denom > 0)
                ? (double)common_count / (double)denom
                : 0.0;
        }

        /* ── Step 4：擴展 accumulator
         *   acc_s ∪= ring_s[l]; acc_t ∪= ring_t[l]
         */
        bs_or_into(w->acc_s, w->ring_s[l]);
        bs_or_into(w->acc_t, w->ring_t[l]);
    }
}

/* ════════════════════════════════════════════════════════════
/* ════════════════════════════════════════════════════════════
 * compute_all_link_ratios — 對所有邊計算 layer 1..layers 的 ratio
 * Python HETA 風格：先建 ego table，再對每邊跑環狀分層演算法。
 * 邊以 CSR (u<v) 順序枚舉，回傳 out[edge_idx*layers + (l-1)]。
 * 空圖回傳 NULL；呼叫者負責 free。
 * ════════════════════════════════════════════════════════════ */
double *compute_all_link_ratios(const Graph *g, int layers) {
    int m = g->m;
    if (m == 0 || layers < 1) return NULL;

    double *out = (double *)xmalloc((size_t)m * (size_t)layers * sizeof(double));
    for (size_t i = 0; i < (size_t)m * (size_t)layers; i++) out[i] = 0.0;

    EgoTable *et = ego_table_alloc(g, layers);
    RingWorkspace *w = ring_workspace_alloc(g->n, layers);

    int edge_idx = 0;
    for (int u = 0; u < g->n; u++) {
        for (int i = g->offset[u]; i < g->offset[u + 1]; i++) {
            int v = g->adj_list[i];
            if (u >= v) continue;
            compute_edge_ratios(g, u, v, layers, et, w,
                                &out[(size_t)edge_idx * layers]);
            edge_idx++;
        }
    }

    ring_workspace_free(w);
    ego_table_free(et);

    return out;
}
