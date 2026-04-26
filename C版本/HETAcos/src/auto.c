#define _GNU_SOURCE
/*
 * auto.c — HETA 批次自動執行入口
 * ================================
 *
 * 用途：
 *   一鍵掃描 IN_DIR 下所有 *.net（Pajek 格式）檔案，
 *   逐一執行 HETA 演算法，將每個結果寫入 OUT_DIR/<name>_result.txt，
 *   最後印出總彙整表（每個網路一列）。
 *
 * 與 main.c 差異：
 *   - main.c：命令列指定單一 edge_file，輸入是 "u v" 邊列表
 *   - auto.c：無參數，掃整個目錄，輸入是 Pajek (.net) 格式
 *
 * 編譯：
 *   gcc -std=c11 -O3 -fopenmp -Iinclude src/{graph,bitset,cnr,randomize,
 *       threshold,heta,auto}.c -o auto -lm
 *
 * 設定（修改下方 #define 後重編）：
 *   IN_DIR     輸入資料夾（含 *.net）
 *   OUT_DIR    輸出資料夾
 *   N_RANDOM   隨機化網路數量
 *   FILE_EXT   篩選副檔名
 */

#include "heta.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <dirent.h>
#include <sys/stat.h>

/* ════════════════════════════════════════════════════════════
 * 設定常數（修改後需重新編譯）
 * ════════════════════════════════════════════════════════════ */
#define IN_DIR    "data/net"      /* Pajek .net 檔來源 */
#define OUT_DIR   "result"        /* 結果輸出資料夾    */
#define N_RANDOM  1000            /* 論文標準值        */
#define FILE_EXT  ".net"          /* 只處理 .net 檔    */
#define MAX_PATH  1024
#define MAX_FILES 1024            /* 最多處理檔案數    */

/* ════════════════════════════════════════════════════════════
 * 邊列表動態陣列（與 main.c 共用設計）
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

/* 工具：忽略大小寫比對開頭 */
static int starts_with_i(const char *s, const char *prefix) {
    while (*prefix) {
        if (tolower((unsigned char)*s) != tolower((unsigned char)*prefix)) return 0;
        s++; prefix++;
    }
    return 1;
}

/* ════════════════════════════════════════════════════════════
 * load_pajek — 讀取 Pajek .net 檔，回傳 Graph*
 *
 * Pajek 格式概述：
 *   *vertices N            ← 節點區段
 *   1 "label" ...          ← 每行第一欄是節點 id（通常 1-indexed，亦可 0）
 *   ...
 *   *edges                 ← 邊區段（無向）
 *   u v [weight ...]       ← 每行前兩欄是端點 id
 *   *arcs                  ← 邊區段（有向，本程式視為無向）
 *
 * 自動偵測 0 / 1-indexed（用 *vertices 區段或邊端點的最小值決定 offset）
 * 自動去自迴圈與重複邊（用 (min,max) 排序後比對）
 * ════════════════════════════════════════════════════════════ */
static Graph *load_pajek(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "[auto] 無法開啟 %s\n", path);
        return NULL;
    }

    /* 第一遍：抓 *vertices 區段所有 id 與 *edges/*arcs 區段所有邊 */
    EdgeList raw;
    el_init(&raw);
    int min_vid = INT32_MAX;     /* 用來判斷 offset 是 0 還 1 */
    int section = 0;              /* 0=none, 1=vertices, 2=edges */

    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        /* 跳過空白行 */
        char *p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '\0' || *p == '\n' || *p == '\r') continue;
        if (*p == '#') continue;

        /* 判斷區段切換 */
        if (*p == '*') {
            if (starts_with_i(p, "*vertices"))      section = 1;
            else if (starts_with_i(p, "*edges") ||
                     starts_with_i(p, "*arcs"))     section = 2;
            else                                    section = 0;
            continue;
        }

        if (section == 1) {
            /* 節點行：取第一個整數作為 vertex id */
            int vid;
            if (sscanf(p, "%d", &vid) == 1) {
                if (vid < min_vid) min_vid = vid;
            }
        } else if (section == 2) {
            /* 邊行：前兩個整數 */
            int u, v;
            if (sscanf(p, "%d %d", &u, &v) == 2) {
                el_push(&raw, u, v);
            }
        }
    }
    fclose(fp);

    /* 計算 offset：min(節點最小id, 邊端點最小id) */
    int edge_min = INT32_MAX;
    for (int i = 0; i < raw.m; i++) {
        if (raw.us[i] < edge_min) edge_min = raw.us[i];
        if (raw.vs[i] < edge_min) edge_min = raw.vs[i];
    }
    int offset = min_vid;
    if (edge_min < offset) offset = edge_min;
    if (offset != 0 && offset != 1) offset = 0;  /* 安全保險 */

    /* 把 1-indexed 平移到 0-indexed，去自迴圈與重複邊 */
    EdgeList el;
    el_init(&el);

    /* 簡易去重：因為原始邊不一定排序，用 hash 太重，
     * 改用「正規化方向 (u<v) 後排序」的方式留待 graph_build 容錯，
     * 不過 Pajek 通常本來就無重複，這邊只去自迴圈 */
    for (int i = 0; i < raw.m; i++) {
        int u = raw.us[i] - offset;
        int v = raw.vs[i] - offset;
        if (u < 0 || v < 0) continue;   /* 安全略過 */
        if (u == v) continue;            /* 去自迴圈 */
        if (u > v) { int t = u; u = v; v = t; }
        el_push(&el, u, v);
    }
    el_free(&raw);

    fprintf(stdout, "[auto] %s  →  edges=%d  offset=%d\n",
            path, el.m, offset);

    Graph *g = graph_build(el.max_node + 1, el.m, el.us, el.vs);
    el_free(&el);
    return g;
}

/* ════════════════════════════════════════════════════════════
 * 輸出格式（精簡複用 main.c 的版面）
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

static void write_result(const char *out_path, const Graph *g,
                         const HetaResult *result, const HetaSummary *s,
                         double elapsed) {
    FILE *fp = fopen(out_path, "w");
    if (!fp) { fprintf(stderr, "[auto] 無法寫入 %s\n", out_path); return; }

    fprintf(fp, "HETA 分析結果\n");
    fprintf(fp, "================================================================================\n");
    fprintf(fp, "節點數     : %d\n", g->n);
    fprintf(fp, "邊數       : %d\n", g->m);
    fprintf(fp, "kmax       : %d\n", result->kmax);
    fprintf(fp, "n_random   : %d\n", N_RANDOM);
    fprintf(fp, "執行時間   : %.3f 秒\n", elapsed);
    fprintf(fp, "================================================================================\n\n");

    fprintf(fp, "分類統計\n");
    fprintf(fp, "--------------------------------------------------------------------------------\n");
    fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "silk link",     s->n_silk,   s->pct_silk);
    fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "bond link",     s->n_bond,   s->pct_bond);
    fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "local bridge",  s->n_local,  s->pct_local);
    fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "global bridge", s->n_global, s->pct_global);
    fprintf(fp, "--------------------------------------------------------------------------------\n\n");

    fprintf(fp, "逐邊分類結果\n");
    fprintf(fp, "--------------------------------------------------------------------------------\n");
    fprintf(fp, "  %-6s  %-6s  %-22s  %s\n", "u", "v", "type", "layer");
    fprintf(fp, "  %-6s  %-6s  %-22s  %s\n", "------", "------", "----------------------", "-----");
    for (int i = 0; i < result->m; i++) {
        const Edge *e = &result->edges[i];
        fprintf(fp, "  %-6d  %-6d  %-22s  %d\n",
                e->u, e->v, type_str(e->type, e->layer), e->layer);
    }
    fprintf(fp, "--------------------------------------------------------------------------------\n");

    fclose(fp);
}

/* ════════════════════════════════════════════════════════════
 * 篩出符合 FILE_EXT 的檔名 → 字典序排序
 * ════════════════════════════════════════════════════════════ */
static int has_ext(const char *name, const char *ext) {
    size_t n = strlen(name), e = strlen(ext);
    if (n <= e) return 0;
    return strcmp(name + (n - e), ext) == 0;
}
static int cmp_str(const void *a, const void *b) {
    return strcmp(*(const char **)a, *(const char **)b);
}

/* ════════════════════════════════════════════════════════════
 * main — 批次入口
 * ════════════════════════════════════════════════════════════ */
int main(void) {
    /* 建立輸出資料夾 */
#if defined(_WIN32) || defined(_WIN64)
    mkdir(OUT_DIR);
#else
    mkdir(OUT_DIR, 0755);
#endif

    /* 開資料夾 */
    DIR *dir = opendir(IN_DIR);
    if (!dir) {
        fprintf(stderr, "[auto] 無法開啟資料夾 %s\n", IN_DIR);
        return EXIT_FAILURE;
    }

    /* 收集所有 .net 檔名 */
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

    /* 暫存每筆統計，最後印總彙整 */
    typedef struct {
        char  name[128];
        int   n, m, kmax;
        int   silk, bond, lb, gb;
        double sec;
        int   ok;
    } Row;
    Row *rows = (Row *)calloc(n_files, sizeof(Row));

    double total_sec = 0.0;
    int    n_ok = 0;

    for (int i = 0; i < n_files; i++) {
        char in_path[MAX_PATH], out_path[MAX_PATH], stem[128];
        snprintf(in_path,  sizeof(in_path),  "%s/%s", IN_DIR, files[i]);

        /* 取 stem（去副檔名）*/
        strncpy(stem, files[i], sizeof(stem) - 1);
        stem[sizeof(stem) - 1] = '\0';
        char *dot = strrchr(stem, '.');
        if (dot) *dot = '\0';
        snprintf(out_path, sizeof(out_path), "%s/%s_result.txt", OUT_DIR, stem);

        fprintf(stdout, "\n[%d/%d] %s\n", i + 1, n_files, files[i]);

        Graph *g = load_pajek(in_path);
        if (!g) {
            strncpy(rows[i].name, stem, sizeof(rows[i].name) - 1);
            rows[i].ok = 0;
            continue;
        }

        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        HetaResult *res = heta_run(g, N_RANDOM);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double sec = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

        HetaSummary s = heta_summarize(res);
        write_result(out_path, g, res, &s, sec);

        /* 紀錄該列 */
        strncpy(rows[i].name, stem, sizeof(rows[i].name) - 1);
        rows[i].n    = g->n;
        rows[i].m    = g->m;
        rows[i].kmax = res->kmax;
        rows[i].silk = s.n_silk;
        rows[i].bond = s.n_bond;
        rows[i].lb   = s.n_local;
        rows[i].gb   = s.n_global;
        rows[i].sec  = sec;
        rows[i].ok   = 1;

        total_sec += sec;
        n_ok++;

        fprintf(stdout, "    n=%d  m=%d  kmax=%d  silk=%d  bond=%d  LB=%d  GB=%d  (%.2fs) → %s\n",
                g->n, g->m, res->kmax, s.n_silk, s.n_bond, s.n_local, s.n_global,
                sec, out_path);

        heta_result_free(res);
        graph_free(g);
        free(files[i]);
    }

    /* ── 總彙整：CSV + 文字表 ── */
    char sum_csv[MAX_PATH];
    snprintf(sum_csv, sizeof(sum_csv), "%s/summary.csv", OUT_DIR);
    FILE *fc = fopen(sum_csv, "w");
    if (fc) {
        fprintf(fc, "network,nodes,edges,kmax,n_random,runtime_sec,"
                    "silk,silk_pct,bond,bond_pct,"
                    "local_bridge,local_bridge_pct,"
                    "global_bridge,global_bridge_pct\n");
        for (int i = 0; i < n_files; i++) {
            if (!rows[i].ok) continue;
            int total = rows[i].silk + rows[i].bond + rows[i].lb + rows[i].gb;
            double inv = total > 0 ? 100.0 / total : 0.0;
            fprintf(fc, "%s,%d,%d,%d,%d,%.3f,%d,%.2f,%d,%.2f,%d,%.2f,%d,%.2f\n",
                    rows[i].name, rows[i].n, rows[i].m, rows[i].kmax,
                    N_RANDOM, rows[i].sec,
                    rows[i].silk, rows[i].silk * inv,
                    rows[i].bond, rows[i].bond * inv,
                    rows[i].lb,   rows[i].lb   * inv,
                    rows[i].gb,   rows[i].gb   * inv);
        }
        fclose(fc);
    }

    fprintf(stdout, "\n════════════════════════════════════════════════════════\n");
    fprintf(stdout, " 總彙整 (%d / %d 成功，總耗時 %.2f 秒)\n", n_ok, n_files, total_sec);
    fprintf(stdout, "════════════════════════════════════════════════════════\n");
    fprintf(stdout, " %-25s %5s %6s %5s %5s %5s %5s %5s %7s\n",
            "network", "n", "m", "kmax", "silk", "bond", "LB", "GB", "sec");
    fprintf(stdout, " %-25s %5s %6s %5s %5s %5s %5s %5s %7s\n",
            "-------------------------", "-----", "------", "-----", "-----", "-----", "-----", "-----", "-------");
    for (int i = 0; i < n_files; i++) {
        if (!rows[i].ok) {
            fprintf(stdout, " %-25s   FAILED\n", rows[i].name);
            continue;
        }
        fprintf(stdout, " %-25s %5d %6d %5d %5d %5d %5d %5d %7.2f\n",
                rows[i].name, rows[i].n, rows[i].m, rows[i].kmax,
                rows[i].silk, rows[i].bond, rows[i].lb, rows[i].gb,
                rows[i].sec);
    }
    fprintf(stdout, "════════════════════════════════════════════════════════\n");
    fprintf(stdout, " summary CSV → %s\n", sum_csv);

    free(rows);
    return EXIT_SUCCESS;
}
