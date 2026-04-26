/*
 * network.c — WS 與 NW 小世界網路生成
 * ======================================
 *
 * 對應 Python main.py 的 build_ws() 和 build_nw() 函式
 *
 * WS (Watts-Strogatz)：
 *   1. 建立 n 節點環形格子，每節點連接左右各 k/2 個鄰居
 *   2. 對每條邊以機率 p 重連到隨機節點
 *
 * NW (Newman-Watts)：
 *   1. 建立 n 節點環形格子（同 WS 步驟 1）
 *   2. 不刪除原有邊，額外新增 p * m 條隨機捷徑
 *
 * 使用 PCG32 作為隨機數生成器，與 randomize.c 一致
 */

#include "rumor.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

/* ────────────────────────────────────────────────────────────
 * PCG32 隨機數生成器（與 randomize.c 相同實作）
 * ──────────────────────────────────────────────────────────── */
typedef struct {
    uint64_t state;
    uint64_t inc;
} NetPCG32;

static void net_pcg32_init(NetPCG32 *rng, uint64_t seed, uint64_t seq) {
    rng->state = 0u;
    rng->inc   = (seq << 1u) | 1u;
    rng->state = rng->state * 6364136223846793005ULL + rng->inc;
    rng->state += seed;
    rng->state = rng->state * 6364136223846793005ULL + rng->inc;
}

static inline uint32_t net_pcg32_next(NetPCG32 *rng) {
    uint64_t old = rng->state;
    rng->state = old * 6364136223846793005ULL + rng->inc;
    uint32_t xorshifted = (uint32_t)(((old >> 18u) ^ old) >> 27u);
    uint32_t rot        = (uint32_t)(old >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
}

static inline int net_pcg32_bound(NetPCG32 *rng, uint32_t bound) {
    uint64_t m = (uint64_t)net_pcg32_next(rng) * (uint64_t)bound;
    return (int)(m >> 32);
}

static inline double net_pcg32_double(NetPCG32 *rng) {
    return (double)net_pcg32_next(rng) / 4294967296.0;
}

/* 安全 malloc */
static void *xmalloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "[network] malloc failed\n"); exit(1); }
    return p;
}

/* ────────────────────────────────────────────────────────────
 * 邊列表動態陣列（用於建構網路）
 * ──────────────────────────────────────────────────────────── */
typedef struct {
    int *us;
    int *vs;
    int  m;
    int  capacity;
} NetEdgeList;

static void nel_init(NetEdgeList *el, int cap) {
    el->capacity = cap;
    el->m = 0;
    el->us = (int *)xmalloc(cap * sizeof(int));
    el->vs = (int *)xmalloc(cap * sizeof(int));
}

static void nel_push(NetEdgeList *el, int u, int v) {
    if (el->m >= el->capacity) {
        el->capacity *= 2;
        el->us = (int *)realloc(el->us, el->capacity * sizeof(int));
        el->vs = (int *)realloc(el->vs, el->capacity * sizeof(int));
        if (!el->us || !el->vs) {
            fprintf(stderr, "[network] realloc failed\n"); exit(1);
        }
    }
    /* 正規化方向：u < v */
    if (u > v) { int t = u; u = v; v = t; }
    el->us[el->m] = u;
    el->vs[el->m] = v;
    el->m++;
}

static void nel_free(NetEdgeList *el) {
    free(el->us);
    free(el->vs);
}

/* ────────────────────────────────────────────────────────────
 * 鄰接矩陣（用於快速 has_edge 查詢）
 * 使用 bitset 陣列，與 randomize.c 一致
 * ──────────────────────────────────────────────────────────── */
static int has_edge_matrix(Bitset **adj_bs, int u, int v) {
    return bs_test(adj_bs[u], v);
}

static void add_edge_matrix(Bitset **adj_bs, int u, int v) {
    bs_set(adj_bs[u], v);
    bs_set(adj_bs[v], u);
}

/* ════════════════════════════════════════════════════════════
 * build_ws_network — 建立 Watts-Strogatz 小世界網路
 *
 * 對應 Python: nx.watts_strogatz_graph(n, k, p, seed=...)
 *
 * 步驟：
 *   1. 建立環形格子：節點 i 連接 i+1, i+2, ..., i+k/2（mod n）
 *   2. 對每個節點 i 的每條右側邊 (i, i+j)，以機率 p 重連
 *      重連目標：隨機選擇一個非自身、非已存在鄰居的節點
 * ════════════════════════════════════════════════════════════ */
Graph *build_ws_network(int n, int k, double p, uint64_t seed) {
    int half_k = k / 2;  /* 每側鄰居數 */

    /* 分配鄰接矩陣 bitset 用於快速查詢 */
    Bitset **adj_bs = adj_bs_alloc(n);

    /* 建立邊列表 */
    NetEdgeList el;
    nel_init(&el, n * half_k);

    /* 步驟 1：建立環形格子 */
    for (int i = 0; i < n; i++) {
        for (int j = 1; j <= half_k; j++) {
            int target = (i + j) % n;
            if (!has_edge_matrix(adj_bs, i, target)) {
                nel_push(&el, i, target);
                add_edge_matrix(adj_bs, i, target);
            }
        }
    }

    /* 步驟 2：重連 */
    NetPCG32 rng;
    net_pcg32_init(&rng, seed, seed ^ 0xABCD1234ULL);

    for (int i = 0; i < n; i++) {
        for (int j = 1; j <= half_k; j++) {
            if (net_pcg32_double(&rng) < p) {
                int old_target = (i + j) % n;

                /* 嘗試找到新的目標節點 */
                int attempts = 0;
                int new_target;
                do {
                    new_target = net_pcg32_bound(&rng, (uint32_t)n);
                    attempts++;
                } while ((new_target == i ||
                          has_edge_matrix(adj_bs, i, new_target)) &&
                         attempts < n * 2);

                if (attempts >= n * 2) continue;  /* 放棄此次重連 */

                /* 移除舊邊的 bitset 記錄 */
                adj_bs[i]->words[old_target / BS_WORD_BITS]
                    &= ~(1ULL << (old_target % BS_WORD_BITS));
                adj_bs[old_target]->words[i / BS_WORD_BITS]
                    &= ~(1ULL << (i % BS_WORD_BITS));

                /* 在邊列表中替換（找到舊邊替換）*/
                int u_norm = (i < old_target) ? i : old_target;
                int v_norm = (i < old_target) ? old_target : i;
                for (int e = 0; e < el.m; e++) {
                    if (el.us[e] == u_norm && el.vs[e] == v_norm) {
                        int nu = (i < new_target) ? i : new_target;
                        int nv = (i < new_target) ? new_target : i;
                        el.us[e] = nu;
                        el.vs[e] = nv;
                        break;
                    }
                }

                /* 加入新邊的 bitset 記錄 */
                add_edge_matrix(adj_bs, i, new_target);
            }
        }
    }

    /* 去重：使用 bitset 重建乾淨的邊列表 */
    NetEdgeList clean;
    nel_init(&clean, el.m);
    for (int u = 0; u < n; u++) {
        for (int word_i = 0; word_i < adj_bs[u]->n_words; word_i++) {
            uint64_t w = adj_bs[u]->words[word_i];
            while (w) {
                int bit = __builtin_ctzll(w);
                int v = word_i * (int)BS_WORD_BITS + bit;
                if (u < v) {
                    nel_push(&clean, u, v);
                }
                w &= w - 1;
            }
        }
    }

    /* 建立 CSR 圖 */
    Graph *g = graph_build(n, clean.m, clean.us, clean.vs);

    /* 清理 */
    adj_bs_free(adj_bs, n);
    nel_free(&el);
    nel_free(&clean);

    return g;
}

/* ════════════════════════════════════════════════════════════
 * build_nw_network — 建立 Newman-Watts 小世界網路
 *
 * 對應 Python main.py 的 build_nw() 函式
 *
 * 步驟：
 *   1. 建立環形格子（p=0 的 WS）
 *   2. 在原有邊基礎上新增 p * m 條隨機捷徑（不刪除原有邊）
 * ════════════════════════════════════════════════════════════ */
Graph *build_nw_network(int n, int k, double p, uint64_t seed) {
    (void)seed;  /* NW 使用固定 seed=42 以匹配 Python 版本 */
    int half_k = k / 2;

    /* 分配鄰接矩陣 bitset */
    Bitset **adj_bs = adj_bs_alloc(n);

    /* 建立邊列表 */
    NetEdgeList el;
    nel_init(&el, n * half_k * 2);

    /* 步驟 1：建立環形格子（無重連） */
    for (int i = 0; i < n; i++) {
        for (int j = 1; j <= half_k; j++) {
            int target = (i + j) % n;
            if (!has_edge_matrix(adj_bs, i, target)) {
                nel_push(&el, i, target);
                add_edge_matrix(adj_bs, i, target);
            }
        }
    }

    int base_edges = el.m;
    int shortcuts = (int)(p * base_edges);

    /* 步驟 2：新增隨機捷徑 */
    NetPCG32 rng;
    net_pcg32_init(&rng, 42, 42 ^ 0xFEDCBA98ULL);  /* 對應 Python rng seed=42 */

    int added = 0;
    int attempts = 0;
    int max_attempts = shortcuts * 100;

    while (added < shortcuts && attempts < max_attempts) {
        int u = net_pcg32_bound(&rng, (uint32_t)n);
        int v = net_pcg32_bound(&rng, (uint32_t)n);
        if (u == v) { attempts++; continue; }

        if (!has_edge_matrix(adj_bs, u, v)) {
            nel_push(&el, u, v);
            add_edge_matrix(adj_bs, u, v);
            added++;
        }
        attempts++;
    }

    /* 建立 CSR 圖 */
    Graph *g = graph_build(n, el.m, el.us, el.vs);

    /* 清理 */
    adj_bs_free(adj_bs, n);
    nel_free(&el);

    return g;
}
