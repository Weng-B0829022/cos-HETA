/*
 * auto.c — HETA 批次自動執行（含節點座標與 R^1 輸出）
 * =====================================================
 *
 * 與 main.c 差異：
 *   - main.c：命令列指定單一 edge_file，輸入是 "u v" 邊列表
 *   - auto.c：無參數，掃整個 IN_DIR 資料夾，輸入是 Pajek (.net) 格式
 *
 * 結果檔額外輸出（為了 Fig. 4–7、12 視覺化）：
 *   - 節點座標表（從 Pajek vertices 抽出 x,y）
 *   - 逐邊 R^1（第一層共同鄰居比，論文公式 (1)）
 *   - 原始邊列表（即逐邊分類結果中的 u,v 兩欄）
 *
 * 編譯：
 *   gcc -std=c11 -O3 -fopenmp -Iinclude src/{graph,bitset,cnr,randomize,
 *       threshold,heta,auto}.c -o auto -lm
 */

#define _GNU_SOURCE
#include "heta.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>

/* ════════════════════════════════════════════════════════════
 * 設定常數（修改後需重新編譯）
 * ════════════════════════════════════════════════════════════ */
#define IN_DIR    "data/net"
#define OUT_DIR   "result"
#define N_RANDOM  1000
#define FILE_EXT  ".net"
#define MAX_PATH  1024
#define MAX_FILES 1024

/* ════════════════════════════════════════════════════════════
 * 邊列表動態陣列
 * ════════════════════════════════════════════════════════════ */
typedef struct {
    int *us, *vs;
    int  m, capacity, max_node;
} EdgeList;

static void el_init(EdgeList *e) {
    e->capacity = 1024; e->m = 0; e->max_node = -1;
    e->us = (int *)malloc(e->capacity * sizeof(int));
    e->vs = (int *)malloc(e->capacity * sizeof(int));
}
static void el_push(EdgeList *e, int u, int v) {
    if (e->m >= e->capacity) {
        e->capacity *= 2;
        e->us = (int *)realloc(e->us, e->capacity * sizeof(int));
        e->vs = (int *)realloc(e->vs, e->capacity * sizeof(int));
    }
    e->us[e->m] = u; e->vs[e->m] = v; e->m++;
    if (u > e->max_node) e->max_node = u;
    if (v > e->max_node) e->max_node = v;
}
static void el_free(EdgeList *e) { free(e->us); free(e->vs); }

/* ════════════════════════════════════════════════════════════
 * 節點座標暫存（Pajek vertices 區段抽出）
 * ════════════════════════════════════════════════════════════ */
typedef struct {
    int    *vid;       /* 原始 Pajek id（0 或 1 起算）       */
    double *x;
    double *y;
    int     n;
    int     capacity;
} VertexList;

static void vl_init(VertexList *v) {
    v->capacity = 1024; v->n = 0;
    v->vid = (int *)malloc(v->capacity * sizeof(int));
    v->x   = (double *)malloc(v->capacity * sizeof(double));
    v->y   = (double *)malloc(v->capacity * sizeof(double));
}
static void vl_push(VertexList *v, int id, double x, double y) {
    if (v->n >= v->capacity) {
        v->capacity *= 2;
        v->vid = (int *)realloc(v->vid, v->capacity * sizeof(int));
        v->x   = (double *)realloc(v->x,   v->capacity * sizeof(double));
        v->y   = (double *)realloc(v->y,   v->capacity * sizeof(double));
    }
    v->vid[v->n] = id; v->x[v->n] = x; v->y[v->n] = y; v->n++;
}
static void vl_free(VertexList *v) { free(v->vid); free(v->x); free(v->y); }

/* ════════════════════════════════════════════════════════════
 * 帶座標的網路資料
 * ════════════════════════════════════════════════════════════ */
typedef struct {
    Graph  *g;
    double *x;          /* 0-indexed 後的座標，大小 g->n */
    double *y;
    int     has_coord;  /* 是否有「非全為 0」的座標      */
    int     offset;     /* 0 或 1（原 Pajek 起始 id）    */
} NetData;

static int starts_with_i(const char *s, const char *prefix) {
    while (*prefix) {
        if (tolower((unsigned char)*s) != tolower((unsigned char)*prefix)) return 0;
        s++; prefix++;
    }
    return 1;
}

/* ════════════════════════════════════════════════════════════
 * load_pajek — 讀取 .net 檔，回傳 NetData*
 *
 * Pajek 節點行格式：
 *   <vid> "<label>" <x> <y> [shape] [...]
 *   或     <vid> <label> <x> <y> ...
 *
 * 自動偵測 0/1-indexed
 * ════════════════════════════════════════════════════════════ */
static NetData *load_pajek(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "[auto] 無法開啟 %s\n", path);
        return NULL;
    }

    EdgeList raw_edges; el_init(&raw_edges);
    VertexList raw_verts; vl_init(&raw_verts);
    int min_vid = INT32_MAX;
    int section = 0;

    char line[2048];
    while (fgets(line, sizeof(line), fp)) {
        char *p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '\0' || *p == '\n' || *p == '\r' || *p == '#') continue;

        if (*p == '*') {
            if (starts_with_i(p, "*vertices"))      section = 1;
            else if (starts_with_i(p, "*edges") ||
                     starts_with_i(p, "*arcs"))     section = 2;
            else                                    section = 0;
            continue;
        }

        if (section == 1) {
            /* 節點行：vid [label] x y ... */
            int vid;
            double x = 0.0, y = 0.0;
            char rest[2048];
            if (sscanf(p, "%d %2047[^\n]", &vid, rest) >= 1) {
                /* 嘗試從 rest 中找出 x, y（跳過 quoted label）*/
                char *q = rest;
                /* 跳過 label：可能是 "..." 或 一個 token */
                while (*q && isspace((unsigned char)*q)) q++;
                if (*q == '"') {
                    q++;
                    while (*q && *q != '"') q++;
                    if (*q == '"') q++;
                } else {
                    while (*q && !isspace((unsigned char)*q)) q++;
                }
                /* 從 label 之後嘗試讀兩個 double */
                sscanf(q, "%lf %lf", &x, &y);
                vl_push(&raw_verts, vid, x, y);
                if (vid < min_vid) min_vid = vid;
            }
        } else if (section == 2) {
            int u, v;
            if (sscanf(p, "%d %d", &u, &v) == 2) {
                el_push(&raw_edges, u, v);
            }
        }
    }
    fclose(fp);

    /* 偵測 offset (0 or 1) */
    int edge_min = INT32_MAX;
    for (int i = 0; i < raw_edges.m; i++) {
        if (raw_edges.us[i] < edge_min) edge_min = raw_edges.us[i];
        if (raw_edges.vs[i] < edge_min) edge_min = raw_edges.vs[i];
    }
    int offset = min_vid;
    if (edge_min < offset) offset = edge_min;
    if (offset != 0 && offset != 1) offset = 0;

    /* 平移 + 去自迴圈 */
    EdgeList el; el_init(&el);
    for (int i = 0; i < raw_edges.m; i++) {
        int u = raw_edges.us[i] - offset;
        int v = raw_edges.vs[i] - offset;
        if (u < 0 || v < 0 || u == v) continue;
        if (u > v) { int t = u; u = v; v = t; }
        el_push(&el, u, v);
    }
    el_free(&raw_edges);

    /* n = max(vertex 區段最大 id, 邊端點最大 id) + 1 */
    int n = el.max_node + 1;
    if (raw_verts.n > 0) {
        for (int i = 0; i < raw_verts.n; i++) {
            int idx = raw_verts.vid[i] - offset;
            if (idx + 1 > n) n = idx + 1;
        }
    }
    if (n <= 0) n = 1;

    /* 建立座標陣列（無對應 vertex 行的節點留 0,0）*/
    double *xs = (double *)calloc(n, sizeof(double));
    double *ys = (double *)calloc(n, sizeof(double));
    int has_coord = 0;
    for (int i = 0; i < raw_verts.n; i++) {
        int idx = raw_verts.vid[i] - offset;
        if (idx >= 0 && idx < n) {
            xs[idx] = raw_verts.x[i];
            ys[idx] = raw_verts.y[i];
            if (raw_verts.x[i] != 0.0 || raw_verts.y[i] != 0.0) has_coord = 1;
        }
    }
    vl_free(&raw_verts);

    fprintf(stdout, "[auto] %s  →  n=%d edges=%d offset=%d coord=%s\n",
            path, n, el.m, offset, has_coord ? "yes" : "no(全 0)");

    Graph *g = graph_build(n, el.m, el.us, el.vs);
    el_free(&el);

    NetData *nd = (NetData *)malloc(sizeof(NetData));
    nd->g = g; nd->x = xs; nd->y = ys;
    nd->has_coord = has_coord; nd->offset = offset;
    return nd;
}

static void net_data_free(NetData *nd) {
    if (!nd) return;
    graph_free(nd->g);
    free(nd->x); free(nd->y);
    free(nd);
}

/* ════════════════════════════════════════════════════════════
 * compute_r1_all — 對每條邊算 R^1（論文公式 (1)）
 *
 * 回傳 double 陣列，大小 = result->m，需要呼叫者 free
 * ════════════════════════════════════════════════════════════ */
static double *compute_r1_all(const Graph *g, const HetaResult *result) {
    int m = result->m;
    double *r1 = (double *)malloc(m * sizeof(double));

    BFSBuffer *buf = bfs_buffer_alloc(g->n);
    Bitset *bs1 = bs_alloc(g->n);
    Bitset *bs2 = bs_alloc(g->n);
    Bitset *bs3 = bs_alloc(g->n);   /* k=1 用不到 bs_k1，但 API 需要傳 */
    Bitset *bs4 = bs_alloc(g->n);

    for (int i = 0; i < m; i++) {
        const Edge *e = &result->edges[i];
        r1[i] = compute_cnr(g, e->u, e->v, 1, buf, bs1, bs2, bs3, bs4);
    }

    bs_free(bs1); bs_free(bs2); bs_free(bs3); bs_free(bs4);
    bfs_buffer_free(buf);
    return r1;
}

/* ════════════════════════════════════════════════════════════
 * 比較函式（給 qsort 用）
 * ════════════════════════════════════════════════════════════ */
static int cmp_d(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}
static int cmp_str(const void *a, const void *b) {
    return strcmp(*(const char **)a, *(const char **)b);
}

/* ════════════════════════════════════════════════════════════
 * 輸出：類型字串
 * ════════════════════════════════════════════════════════════ */
static const char *type_str(LinkType t, int layer) {
    static char buf[32];
    switch (t) {
        case LINK_SILK:          return "silk";
        case LINK_BOND:          return "bond";
        case LINK_LOCAL_BRIDGE:
            snprintf(buf, sizeof(buf), "local_bridge_k%d", layer);
            return buf;
        case LINK_GLOBAL_BRIDGE: return "global_bridge";
        default:                 return "unclassified";
    }
}

/* ════════════════════════════════════════════════════════════
 * write_result — 輸出完整結果 txt（含座標表 + R^1 欄）
 * ════════════════════════════════════════════════════════════ */
static void write_result(const char *out_path, const NetData *nd,
                         const HetaResult *result, const HetaSummary *s,
                         const double *r1, double elapsed) {
    FILE *fp = fopen(out_path, "w");
    if (!fp) { fprintf(stderr, "[auto] 無法寫入 %s\n", out_path); return; }

    const Graph *g = nd->g;

    /* ── 標頭 ── */
    fprintf(fp, "HETA 分析結果\n");
    fprintf(fp, "================================================================================\n");
    fprintf(fp, "節點數     : %d\n", g->n);
    fprintf(fp, "邊數       : %d\n", g->m);
    fprintf(fp, "kmax       : %d\n", result->kmax);
    fprintf(fp, "n_random   : %d\n", N_RANDOM);
    fprintf(fp, "執行時間   : %.3f 秒\n", elapsed);
    fprintf(fp, "原始 indexing : %d-indexed\n", nd->offset);
    fprintf(fp, "節點座標   : %s\n", nd->has_coord ? "有（從 Pajek 抽出）" : "無（全 0）");
    fprintf(fp, "================================================================================\n\n");

    /* ── 統計摘要 ── */
    fprintf(fp, "分類統計\n");
    fprintf(fp, "--------------------------------------------------------------------------------\n");
    fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "silk link",     s->n_silk,   s->pct_silk);
    fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "bond link",     s->n_bond,   s->pct_bond);
    fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "local bridge",  s->n_local,  s->pct_local);
    fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "global bridge", s->n_global, s->pct_global);
    fprintf(fp, "--------------------------------------------------------------------------------\n\n");

    /* ── R^1 統計（min/max/mean/median）── */
    int m = result->m;
    if (m > 0) {
        double rmin = r1[0], rmax = r1[0], rsum = 0.0;
        for (int i = 0; i < m; i++) {
            if (r1[i] < rmin) rmin = r1[i];
            if (r1[i] > rmax) rmax = r1[i];
            rsum += r1[i];
        }
        double rmean = rsum / m;
        /* 中位數需要排序，先複製再排序 */
        double *sorted = (double *)malloc(m * sizeof(double));
        memcpy(sorted, r1, m * sizeof(double));
        qsort(sorted, m, sizeof(double), cmp_d);
        double rmed = (m % 2 == 1) ? sorted[m/2]
                                    : 0.5 * (sorted[m/2 - 1] + sorted[m/2]);
        free(sorted);

        fprintf(fp, "R^1 (第一層共同鄰居比) 統計\n");
        fprintf(fp, "--------------------------------------------------------------------------------\n");
        fprintf(fp, "  min    : %.6f\n", rmin);
        fprintf(fp, "  max    : %.6f\n", rmax);
        fprintf(fp, "  mean   : %.6f\n", rmean);
        fprintf(fp, "  median : %.6f\n", rmed);
        fprintf(fp, "--------------------------------------------------------------------------------\n\n");
    }

    /* ── 節點座標表 ── */
    fprintf(fp, "節點座標 (Vertices)\n");
    fprintf(fp, "--------------------------------------------------------------------------------\n");
    fprintf(fp, "  %-6s  %-12s  %-12s\n", "id", "x", "y");
    fprintf(fp, "  %-6s  %-12s  %-12s\n", "------", "------------", "------------");
    for (int i = 0; i < g->n; i++) {
        fprintf(fp, "  %-6d  %12.6f  %12.6f\n", i, nd->x[i], nd->y[i]);
    }
    fprintf(fp, "--------------------------------------------------------------------------------\n\n");

    /* ── 逐邊分類結果 + R^1 ── */
    fprintf(fp, "逐邊分類結果 (含 R^1)\n");
    fprintf(fp, "--------------------------------------------------------------------------------\n");
    fprintf(fp, "  %-6s  %-6s  %-22s  %-5s  %-10s\n",
            "u", "v", "type", "layer", "R^1");
    fprintf(fp, "  %-6s  %-6s  %-22s  %-5s  %-10s\n",
            "------", "------", "----------------------", "-----", "----------");
    for (int i = 0; i < result->m; i++) {
        const Edge *e = &result->edges[i];
        fprintf(fp, "  %-6d  %-6d  %-22s  %-5d  %.6f\n",
                e->u, e->v, type_str(e->type, e->layer), e->layer, r1[i]);
    }
    fprintf(fp, "--------------------------------------------------------------------------------\n");

    fclose(fp);
}

/* ════════════════════════════════════════════════════════════
 * 工具：判副檔名
 * ════════════════════════════════════════════════════════════ */
static int has_ext(const char *name, const char *ext) {
    size_t n = strlen(name), e = strlen(ext);
    if (n <= e) return 0;
    return strcmp(name + (n - e), ext) == 0;
}

/* ════════════════════════════════════════════════════════════
 * main
 * ════════════════════════════════════════════════════════════ */
int main(void) {
#if defined(_WIN32) || defined(_WIN64)
    mkdir(OUT_DIR);
#else
    mkdir(OUT_DIR, 0755);
#endif

    DIR *dir = opendir(IN_DIR);
    if (!dir) {
        fprintf(stderr, "[auto] 無法開啟資料夾 %s\n", IN_DIR);
        return EXIT_FAILURE;
    }
    char *files[MAX_FILES];
    int   n_files = 0;
    struct dirent *de;
    while ((de = readdir(dir)) != NULL) {
        if (!has_ext(de->d_name, FILE_EXT)) continue;
        if (n_files >= MAX_FILES) break;
        files[n_files++] = strdup(de->d_name);
    }
    closedir(dir);
    qsort(files, n_files, sizeof(char *), cmp_str);

    if (n_files == 0) {
        fprintf(stderr, "[auto] 在 %s 找不到 *%s 檔\n", IN_DIR, FILE_EXT);
        return EXIT_FAILURE;
    }

    fprintf(stdout, "════════════════════════════════════════════════════════\n");
    fprintf(stdout, " HETA 批次執行：%d 個檔案 (n_random=%d)\n", n_files, N_RANDOM);
    fprintf(stdout, " 來源：%s/\n 輸出：%s/\n", IN_DIR, OUT_DIR);
    fprintf(stdout, "════════════════════════════════════════════════════════\n");

    typedef struct {
        char  name[128];
        int   n, m, kmax;
        int   silk, bond, lb, gb;
        double sec, r1_mean;
        int   ok, has_coord;
    } Row;
    Row *rows = (Row *)calloc(n_files, sizeof(Row));

    double total_sec = 0.0;
    int n_ok = 0;

    for (int i = 0; i < n_files; i++) {
        char in_path[MAX_PATH], out_path[MAX_PATH], stem[128];
        snprintf(in_path, sizeof(in_path), "%s/%s", IN_DIR, files[i]);
        strncpy(stem, files[i], sizeof(stem) - 1);
        stem[sizeof(stem) - 1] = '\0';
        char *dot = strrchr(stem, '.');
        if (dot) *dot = '\0';
        snprintf(out_path, sizeof(out_path), "%s/%s_result.txt", OUT_DIR, stem);

        fprintf(stdout, "\n[%d/%d] %s\n", i + 1, n_files, files[i]);

        NetData *nd = load_pajek(in_path);
        if (!nd) {
            strncpy(rows[i].name, stem, sizeof(rows[i].name) - 1);
            rows[i].ok = 0;
            continue;
        }

        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        HetaResult *res = heta_run(nd->g, N_RANDOM);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double sec = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

        HetaSummary s = heta_summarize(res);

        /* 算 R^1 */
        double *r1 = compute_r1_all(nd->g, res);

        write_result(out_path, nd, res, &s, r1, sec);

        /* 紀錄 */
        double r1_mean = 0.0;
        for (int j = 0; j < res->m; j++) r1_mean += r1[j];
        if (res->m > 0) r1_mean /= res->m;

        strncpy(rows[i].name, stem, sizeof(rows[i].name) - 1);
        rows[i].n = nd->g->n; rows[i].m = nd->g->m; rows[i].kmax = res->kmax;
        rows[i].silk = s.n_silk; rows[i].bond = s.n_bond;
        rows[i].lb = s.n_local; rows[i].gb = s.n_global;
        rows[i].sec = sec; rows[i].r1_mean = r1_mean;
        rows[i].ok = 1; rows[i].has_coord = nd->has_coord;

        total_sec += sec; n_ok++;

        fprintf(stdout,
                "    n=%d m=%d kmax=%d silk=%d bond=%d LB=%d GB=%d  R^1_mean=%.4f  (%.2fs) → %s\n",
                nd->g->n, nd->g->m, res->kmax,
                s.n_silk, s.n_bond, s.n_local, s.n_global,
                r1_mean, sec, out_path);

        free(r1);
        heta_result_free(res);
        net_data_free(nd);
        free(files[i]);
    }

    /* ── summary CSV ── */
    char sum_csv[MAX_PATH];
    snprintf(sum_csv, sizeof(sum_csv), "%s/summary.csv", OUT_DIR);
    FILE *fc = fopen(sum_csv, "w");
    if (fc) {
        fprintf(fc, "network,nodes,edges,kmax,n_random,runtime_sec,has_coord,"
                    "silk,silk_pct,bond,bond_pct,"
                    "local_bridge,local_bridge_pct,"
                    "global_bridge,global_bridge_pct,r1_mean\n");
        for (int i = 0; i < n_files; i++) {
            if (!rows[i].ok) continue;
            int total = rows[i].silk + rows[i].bond + rows[i].lb + rows[i].gb;
            double inv = total > 0 ? 100.0 / total : 0.0;
            fprintf(fc, "%s,%d,%d,%d,%d,%.3f,%d,%d,%.2f,%d,%.2f,%d,%.2f,%d,%.2f,%.6f\n",
                    rows[i].name, rows[i].n, rows[i].m, rows[i].kmax,
                    N_RANDOM, rows[i].sec, rows[i].has_coord,
                    rows[i].silk, rows[i].silk * inv,
                    rows[i].bond, rows[i].bond * inv,
                    rows[i].lb,   rows[i].lb   * inv,
                    rows[i].gb,   rows[i].gb   * inv,
                    rows[i].r1_mean);
        }
        fclose(fc);
    }

    fprintf(stdout, "\n════════════════════════════════════════════════════════\n");
    fprintf(stdout, " 總彙整 (%d / %d 成功，總耗時 %.2f 秒)\n", n_ok, n_files, total_sec);
    fprintf(stdout, "════════════════════════════════════════════════════════\n");
    fprintf(stdout, " %-25s %5s %6s %5s %5s %5s %5s %5s %8s %7s\n",
            "network", "n", "m", "kmax", "silk", "bond", "LB", "GB", "R^1_avg", "sec");
    fprintf(stdout, " %-25s %5s %6s %5s %5s %5s %5s %5s %8s %7s\n",
            "-------------------------", "-----", "------", "-----",
            "-----", "-----", "-----", "-----", "--------", "-------");
    for (int i = 0; i < n_files; i++) {
        if (!rows[i].ok) {
            fprintf(stdout, " %-25s   FAILED\n", rows[i].name);
            continue;
        }
        fprintf(stdout, " %-25s %5d %6d %5d %5d %5d %5d %5d %8.4f %7.2f\n",
                rows[i].name, rows[i].n, rows[i].m, rows[i].kmax,
                rows[i].silk, rows[i].bond, rows[i].lb, rows[i].gb,
                rows[i].r1_mean, rows[i].sec);
    }
    fprintf(stdout, "════════════════════════════════════════════════════════\n");
    fprintf(stdout, " summary CSV → %s\n", sum_csv);

    free(rows);
    return EXIT_SUCCESS;
}
