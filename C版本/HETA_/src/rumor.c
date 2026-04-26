/*
 * rumor.c — 謠言擴散模擬（信心度累積 + 衰減）
 * ================================================
 *
 * 對應 Python main.py 的：
 *   - select_seed_nodes_from_global()
 *   - simulate_rumor()
 *
 * 【信心度規則】
 * - 每天計算來自感染鄰居的傳入信心度總和
 * - 若當天有傳入：conf += 傳入量，重置 decay_timer = 0
 * - 若當天無傳入：decay_timer += 1
 *                 conf *= DECAY_RATE
 *                 若 decay_timer >= DECAY_ZERO_DAY：conf = 0
 * - conf >= 1 的節點視為「感染中」，會對外傳播
 */

#include "rumor.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* 安全 malloc */
static void *xmalloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "[rumor] malloc failed\n"); exit(1); }
    return p;
}

static void *xcalloc(size_t n, size_t sz) {
    void *p = calloc(n, sz);
    if (!p) { fprintf(stderr, "[rumor] calloc failed\n"); exit(1); }
    return p;
}

/* ════════════════════════════════════════════════════════════
 * get_confidence_rate — 根據邊類型回傳信心度貢獻率
 * ════════════════════════════════════════════════════════════ */
double get_confidence_rate(LinkType type, const RumorParams *params) {
    switch (type) {
        case LINK_BOND:          return params->conf_bond;
        case LINK_LOCAL_BRIDGE:  return params->conf_local;
        case LINK_GLOBAL_BRIDGE: return params->conf_global;
        case LINK_SILK:          return params->conf_silk;
        default:                 return params->conf_global;  /* fallback */
    }
}

/* ════════════════════════════════════════════════════════════
 * 建立邊類型查找表
 *
 * 為了在模擬中快速查詢邊 (u,v) 的類型，建立一個
 * 從 CSR adj_list index 映射到 LinkType 的陣列。
 *
 * edge_type_map[offset[u] + i] = 邊 (u, adj_list[offset[u]+i]) 的類型
 * ════════════════════════════════════════════════════════════ */
static LinkType *build_edge_type_map(const Graph *g, const HetaResult *heta) {
    int n = g->n;

    /* 分配映射陣列（大小 = 2m，因為無向圖每條邊存兩次）*/
    LinkType *map = (LinkType *)xcalloc(g->offset[n], sizeof(LinkType));

    /* 預設所有邊為 global bridge */
    for (int i = 0; i < g->offset[n]; i++) {
        map[i] = LINK_GLOBAL_BRIDGE;
    }

    /* 填入 HETA 分類結果 */
    for (int e = 0; e < heta->m; e++) {
        int u = heta->edges[e].u;
        int v = heta->edges[e].v;
        LinkType type = heta->edges[e].type;

        /* 找到 u->v 在 adj_list 中的位置 */
        for (int i = g->offset[u]; i < g->offset[u + 1]; i++) {
            if (g->adj_list[i] == v) {
                map[i] = type;
                break;
            }
        }
        /* 找到 v->u 在 adj_list 中的位置 */
        for (int i = g->offset[v]; i < g->offset[v + 1]; i++) {
            if (g->adj_list[i] == u) {
                map[i] = type;
                break;
            }
        }
    }

    return map;
}

/* ════════════════════════════════════════════════════════════
 * select_seeds_from_global — 選取所有 global bridge 端點作為種子
 *
 * 蒐集所有 global bridge 邊的兩端節點（去重），全部投放。
 * 不需要指定投放數量，有幾個 global bridge 端點就投幾個。
 * ════════════════════════════════════════════════════════════ */
int *select_seeds_from_global(const Graph *g, const HetaResult *heta,
                              int n_seeds, int *out_count) {
    (void)n_seeds;   /* 不再使用，保留參數以維持介面相容 */
    int n = g->n;

    /* 用 flag 陣列標記所有 global bridge 端點（自動去重）*/
    int *is_seed = (int *)xcalloc(n, sizeof(int));
    int count = 0;

    for (int e = 0; e < heta->m; e++) {
        if (heta->edges[e].type == LINK_GLOBAL_BRIDGE) {
            int u = heta->edges[e].u;
            int v = heta->edges[e].v;
            if (!is_seed[u]) { is_seed[u] = 1; count++; }
            if (!is_seed[v]) { is_seed[v] = 1; count++; }
        }
    }

    /* 蒐集種子節點 */
    int *seeds = (int *)xmalloc(count * sizeof(int));
    int idx = 0;
    for (int v = 0; v < n; v++) {
        if (is_seed[v]) seeds[idx++] = v;
    }
    *out_count = count;

    fprintf(stdout, "  Global bridge 端點全部投放：共 %d 個種子節點\n", count);
    if (count > 0) {
        fprintf(stdout, "  種子節點：");
        int show = (count < 20) ? count : 20;
        for (int i = 0; i < show; i++) fprintf(stdout, " %d", seeds[i]);
        if (count > 20) fprintf(stdout, " ... (共%d個)", count);
        fprintf(stdout, "\n");
    }

    free(is_seed);
    return seeds;
}

/* ════════════════════════════════════════════════════════════
 * simulate_rumor — 執行謠言擴散模擬
 *
 * 每天流程：
 *   Step A：記錄當天 conf >= 1 的節點比例
 *   Step B：統計傳入 — conf>=1 的節點對每條出邊的鄰居貢獻信心度
 *   Step C：同步更新
 *       有傳入 > 0：conf += 傳入量；decay_timer = 0
 *       無傳入    ：decay_timer += 1
 *                  conf *= DECAY_RATE
 *                  decay_timer >= DECAY_ZERO_DAY → conf = 0
 * ════════════════════════════════════════════════════════════ */
SimResult simulate_rumor(const Graph *g, const HetaResult *heta,
                         const int *seed_nodes, int n_seeds,
                         const RumorParams *params) {
    int n    = g->n;
    int days = params->days;

    /* 建立邊類型查找表 */
    LinkType *edge_type_map = build_edge_type_map(g, heta);

    /* 節點狀態陣列 */
    double *conf        = (double *)xcalloc(n, sizeof(double));
    int    *decay_timer = (int    *)xcalloc(n, sizeof(int));
    double *incoming    = (double *)xmalloc(n * sizeof(double));
    double *new_conf    = (double *)xmalloc(n * sizeof(double));
    int    *new_timer   = (int    *)xmalloc(n * sizeof(int));

    /* 初始化種子節點 */
    for (int i = 0; i < n_seeds; i++) {
        conf[seed_nodes[i]] = 1.0;
        decay_timer[seed_nodes[i]] = 0;
    }

    /* 結果陣列（days + 1 天，包含最後一天的記錄）*/
    SimResult result;
    result.n_days      = days + 1;
    result.daily_ratio = (double *)xmalloc(result.n_days * sizeof(double));

    for (int day = 0; day < days; day++) {
        /* ── Step A：記錄當天感染比例 ── */
        int infected = 0;
        for (int v = 0; v < n; v++) {
            if (conf[v] >= 1.0) infected++;
        }
        result.daily_ratio[day] = (double)infected / (double)n;

        /* ── Step B：計算各節點傳入信心度 ── */
        memset(incoming, 0, n * sizeof(double));

        for (int u = 0; u < n; u++) {
            if (conf[u] < 1.0) continue;  /* 未感染，不傳播 */

            /* 遍歷 u 的所有鄰居 */
            for (int i = g->offset[u]; i < g->offset[u + 1]; i++) {
                int v = g->adj_list[i];
                LinkType etype = edge_type_map[i];
                double rate = get_confidence_rate(etype, params);
                incoming[v] += rate * params->willingness;
            }
        }

        /* ── Step C：同步更新 ── */
        /* 感染節點(conf>=1)使用 infected_decay（更強衰減），
         * 未感染節點使用 decay_rate（正常衰減）。
         * 數學：穩態 conf* = inc / (1 - d)
         *   未感染 d=0.75 → 需 inc≥0.25 才能感染
         *   已感染 d=0.30 → 需 inc≥0.70 才能維持
         * → 感染節點需要更多鄰居支持才能存活，自然產生部分感染均衡 */
        double d_infected = (params->infected_decay > 0.0)
                            ? params->infected_decay
                            : params->decay_rate;

        for (int v = 0; v < n; v++) {
            double inc = incoming[v];

            if (inc > 0.0) {
                /* 有傳入：先衰減再累加，上限 cap=1.0 */
                double d = (conf[v] >= 1.0) ? d_infected : params->decay_rate;
                double updated = conf[v] * d + inc;
                if (updated > 1.0) updated = 1.0;
                new_conf[v]  = updated;
                new_timer[v] = 0;
            } else {
                /* 無傳入：推進計時器，衰減信心度 */
                double d = (conf[v] >= 1.0) ? d_infected : params->decay_rate;
                int t = decay_timer[v] + 1;
                new_timer[v] = t;
                if (t >= params->decay_zero_day) {
                    new_conf[v] = 0.0;           /* 超過上限，歸零 */
                } else {
                    new_conf[v] = conf[v] * d;
                }
            }
        }

        /* 寫回狀態 */
        memcpy(conf,        new_conf,  n * sizeof(double));
        memcpy(decay_timer, new_timer, n * sizeof(int));
    }

    /* 最後一天記錄 */
    {
        int infected = 0;
        for (int v = 0; v < n; v++) {
            if (conf[v] >= 1.0) infected++;
        }
        result.daily_ratio[days] = (double)infected / (double)n;
    }

    /* 清理 */
    free(edge_type_map);
    free(conf);
    free(decay_timer);
    free(incoming);
    free(new_conf);
    free(new_timer);

    return result;
}

/* ════════════════════════════════════════════════════════════
 * print_daily_summary — 每 interval 天列印摘要
 * ════════════════════════════════════════════════════════════ */
void print_daily_summary(const SimResult *ws, const SimResult *nw, int interval) {
    fprintf(stdout, "\n── 每%d天數據摘要 ──\n", interval);
    fprintf(stdout, "%5s | %8s | %8s\n", "Day", "WS (%)", "NW (%)");
    fprintf(stdout, "------------------------------\n");

    int max_days = (ws->n_days > nw->n_days) ? ws->n_days : nw->n_days;
    for (int d = 0; d < max_days; d += interval) {
        double ws_pct = (d < ws->n_days) ? ws->daily_ratio[d] * 100.0 : 0.0;
        double nw_pct = (d < nw->n_days) ? nw->daily_ratio[d] * 100.0 : 0.0;
        fprintf(stdout, "%5d | %7.1f%% | %7.1f%%\n", d, ws_pct, nw_pct);
    }
}

/* ════════════════════════════════════════════════════════════
 * write_result_csv — 將每日結果寫入 CSV 檔
 * ════════════════════════════════════════════════════════════ */
void write_result_csv(const char *path,
                      const SimResult *ws, const SimResult *nw,
                      const RumorParams *params) {
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "[rumor] 無法建立 CSV 檔案：%s\n", path);
        return;
    }

    /* 寫入參數資訊（注釋行）*/
    fprintf(fp, "# Rumor Spreading Model (HETA Edge Types + Fixed Decay Rate)\n");
    fprintf(fp, "# N=%d, k=%d, p=%.2f, seeds=%d, decay_rate=%.2f, "
                "zero_at=%d days\n",
            params->n_nodes, params->k_neighbors, params->rewire_prob,
            params->n_seeds, params->decay_rate, params->decay_zero_day);
    fprintf(fp, "# bond=%.2f, local=%.2f, global=%.2f, silk=%.2f, "
                "willingness=%.2f\n",
            params->conf_bond, params->conf_local, params->conf_global,
            params->conf_silk, params->willingness);

    /* CSV 標頭 */
    fprintf(fp, "day,ws_ratio,nw_ratio,ws_pct,nw_pct\n");

    int max_days = (ws->n_days > nw->n_days) ? ws->n_days : nw->n_days;
    for (int d = 0; d < max_days; d++) {
        double ws_r = (d < ws->n_days) ? ws->daily_ratio[d] : 0.0;
        double nw_r = (d < nw->n_days) ? nw->daily_ratio[d] : 0.0;
        fprintf(fp, "%d,%.6f,%.6f,%.2f,%.2f\n",
                d, ws_r, nw_r, ws_r * 100.0, nw_r * 100.0);
    }

    fclose(fp);
    fprintf(stdout, "[rumor] CSV 已寫入：%s\n", path);
}
