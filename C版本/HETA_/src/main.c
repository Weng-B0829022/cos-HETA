/*
 * main.c — HETA 程式入口
 * =======================
 *
 * 使用方式：
 *   ./heta <edge_file> [n_random]
 *
 * 邊檔案格式（純文字，每行一條邊）：
 *   <節點u id> <節點v id>
 *   節點 id 從 0 開始（0-indexed）
 *   每行以空白或 TAB 分隔
 *   '#' 開頭的行視為注釋，直接跳過
 *
 * 範例輸入檔（karate.txt）：
 *   # Zachary's Karate Club Network
 *   0 1
 *   0 2
 *   0 3
 *   ...
 *
 * 輸出：
 *   每條邊的分類結果（邊 u v type layer）+ 統計摘要
 *   若指定 --csv，輸出 CSV 格式到 stdout
 *
 * 參數：
 *   <edge_file> : 邊列表檔案路徑
 *   [n_random]  : 隨機化網路數量（預設 1000，測試時建議 100）
 */

 #include "heta.h"

 #include <stdio.h>    /* FILE, fopen, fclose, fprintf, fgets  */
 #include <stdlib.h>   /* malloc, realloc, free, atoi, exit    */
 #include <string.h>   /* strcmp, strchr, strrchr, strcat      */
 #include <time.h>     /* clock_gettime, CLOCK_MONOTONIC       */
 #include <sys/stat.h> /* mkdir（建立 result 資料夾）           */
 
 /* ────────────────────────────────────────────────────────────
  * 邊列表讀取結構（動態陣列，初始容量 1024，不足時倍增）
  * ──────────────────────────────────────────────────────────── */
 typedef struct {
     int  *us;        /* 邊起點陣列 */
     int  *vs;        /* 邊終點陣列 */
     int   m;         /* 已讀入的邊數 */
     int   capacity;  /* 目前分配的容量 */
     int   max_node;  /* 出現過的最大節點 id（+1 即節點總數）*/
 } EdgeList;
 
 /* 初始化 EdgeList */
 static void edgelist_init(EdgeList *el) {
     el->capacity = 1024;
     el->m        = 0;
     el->max_node = -1;
     el->us = (int *)malloc(el->capacity * sizeof(int));
     el->vs = (int *)malloc(el->capacity * sizeof(int));
     if (!el->us || !el->vs) {
         fprintf(stderr, "[main] EdgeList malloc failed\n");
         exit(EXIT_FAILURE);
     }
 }
 
 /* 加入一條邊到 EdgeList（必要時倍增容量）*/
 static void edgelist_push(EdgeList *el, int u, int v) {
     /* 容量不足時倍增（動態陣列標準策略）*/
     if (el->m >= el->capacity) {
         el->capacity *= 2;
         el->us = (int *)realloc(el->us, el->capacity * sizeof(int));
         el->vs = (int *)realloc(el->vs, el->capacity * sizeof(int));
         if (!el->us || !el->vs) {
             fprintf(stderr, "[main] EdgeList realloc failed\n");
             exit(EXIT_FAILURE);
         }
     }
     el->us[el->m]  = u;
     el->vs[el->m]  = v;
     el->m++;
     /* 更新最大節點 id */
     if (u > el->max_node) el->max_node = u;
     if (v > el->max_node) el->max_node = v;
 }
 
 /* 釋放 EdgeList */
 static void edgelist_free(EdgeList *el) {
     free(el->us);
     free(el->vs);
 }
 
 /* ════════════════════════════════════════════════════════════
  * load_edge_file — 從文字檔讀取邊列表
  *
  * 支援格式：
  *   - 每行兩個整數（節點 id），空白或 TAB 分隔
  *   - '#' 開頭的行為注釋，跳過
  *   - 空行跳過
  *
  * 自動去除重複邊與自迴圈（u == v 的邊）
  * ════════════════════════════════════════════════════════════ */
 static Graph *load_edge_file(const char *path) {
     FILE *fp = fopen(path, "r");
     if (!fp) {
         fprintf(stderr, "[main] Cannot open '%s'\n", path);
         exit(EXIT_FAILURE);
     }
 
     EdgeList el;
     edgelist_init(&el);
 
     char line[256];   /* 每行最多 255 個字元 */
     int  lineno = 0;
 
     while (fgets(line, sizeof(line), fp)) {
         lineno++;
 
         /* 跳過注釋行（'#' 開頭）*/
         if (line[0] == '#') continue;
 
         /* 跳過空行 */
         if (line[0] == '\n' || line[0] == '\r') continue;
 
         int u, v;
         /* sscanf 解析兩個整數 */
         if (sscanf(line, "%d %d", &u, &v) != 2) {
             /* 解析失敗（可能是 TAB 分隔，嘗試 %d\t%d）*/
             if (sscanf(line, "%d\t%d", &u, &v) != 2) {
                 fprintf(stderr, "[main] Warning: cannot parse line %d: %s",
                         lineno, line);
                 continue;
             }
         }
 
         /* 跳過自迴圈（u == v）*/
         if (u == v) continue;
 
         /* 正規化方向：確保 u < v（無向圖）*/
         if (u > v) { int tmp = u; u = v; v = tmp; }
 
         edgelist_push(&el, u, v);
     }
     fclose(fp);
 
     fprintf(stdout, "[main] Loaded %d edges from '%s'\n", el.m, path);
     fprintf(stdout, "[main] Max node id = %d (n = %d)\n",
             el.max_node, el.max_node + 1);
 
     /* 建立 CSR 圖 */
     Graph *g = graph_build(el.max_node + 1, el.m, el.us, el.vs);
 
     edgelist_free(&el);
     return g;
 }
 
 
 /* ════════════════════════════════════════════════════════════
  * print_usage — 顯示使用說明
  * ════════════════════════════════════════════════════════════ */
 static void print_usage(const char *prog) {
     fprintf(stderr,
         "用法：%s <edge_file> [n_random]\n\n"
         "  edge_file : 邊列表檔案（每行 \"u v\"，# 開頭為注釋）\n"
         "  n_random  : 隨機化網路數量（預設 %d，測試用 100）\n\n"
         "範例：\n"
         "  %s karate.txt 1000\n"
         "  %s karate.txt 100   # 快速測試\n",
         prog, HETA_N_RANDOM_DEFAULT, prog, prog);
 }
 
 /* ════════════════════════════════════════════════════════════
  * build_output_path — 從輸入檔名產生輸出路徑
  *
  * 規則：
  *   輸入  "karate_test.txt"   → "result/karate_test_result.txt"
  *   輸入  "data/my_graph.txt" → "result/my_graph_result.txt"
  *   （去掉目錄部分，只取純檔名，去掉副檔名後加 _result.txt）
  * ════════════════════════════════════════════════════════════ */
 static void build_output_path(const char *edge_file, char *out, size_t out_sz) {
     /* 取得純檔名（去掉路徑）*/
     const char *base = strrchr(edge_file, '/');
     base = base ? base + 1 : edge_file;   /* 若無 '/' 直接用原字串 */
 
     /* 複製純檔名到暫存緩衝區 */
     char name[256];
     strncpy(name, base, sizeof(name) - 1);
     name[sizeof(name) - 1] = '\0';
 
     /* 去掉副檔名（找最後一個 '.'）*/
     char *dot = strrchr(name, '.');
     if (dot) *dot = '\0';
 
     /* 組合輸出路徑：result/<name>_result.txt */
     snprintf(out, out_sz, "result/%s_result.txt", name);
 }
 
 /* ════════════════════════════════════════════════════════════
  * write_result — 將統計摘要與逐邊分類結果寫入 txt 檔
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
                          double elapsed, int n_random) {
     FILE *fp = fopen(out_path, "w");
     if (!fp) {
         fprintf(stderr, "[main] 無法建立輸出檔案：%s\n", out_path);
         return;
     }
 
     /* ── 標頭資訊 ── */
     fprintf(fp, "HETA 分析結果\n");
     fprintf(fp, "================================================================================\n");
     fprintf(fp, "節點數     : %d\n", g->n);
     fprintf(fp, "邊數       : %d\n", g->m);
     fprintf(fp, "kmax       : %d\n", result->kmax);
     fprintf(fp, "n_random   : %d\n", n_random);
     fprintf(fp, "執行時間   : %.3f 秒\n", elapsed);
     fprintf(fp, "================================================================================\n\n");
 
     /* ── 統計摘要 ── */
     fprintf(fp, "分類統計\n");
     fprintf(fp, "--------------------------------------------------------------------------------\n");
     fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "silk link",     s->n_silk,   s->pct_silk);
     fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "bond link",     s->n_bond,   s->pct_bond);
     fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "local bridge",  s->n_local,  s->pct_local);
     fprintf(fp, "  %-14s : %4d 條  (%5.1f%%)\n", "global bridge", s->n_global, s->pct_global);
     fprintf(fp, "--------------------------------------------------------------------------------\n\n");
 
     /* ── 逐邊分類結果 ── */
     fprintf(fp, "逐邊分類結果\n");
     fprintf(fp, "--------------------------------------------------------------------------------\n");
     fprintf(fp, "  %-6s  %-6s  %-22s  %s\n", "u", "v", "type", "layer");
     fprintf(fp, "  %-6s  %-6s  %-22s  %s\n", "──────", "──────", "──────────────────────", "─────");
     for (int i = 0; i < result->m; i++) {
         const Edge *e = &result->edges[i];
         fprintf(fp, "  %-6d  %-6d  %-22s  %d\n",
                 e->u, e->v, type_str(e->type, e->layer), e->layer);
     }
     fprintf(fp, "--------------------------------------------------------------------------------\n");
 
     fclose(fp);
     fprintf(stdout, "[main] 結果已寫入：%s\n", out_path);
 }
 
 /* ════════════════════════════════════════════════════════════
  * main — 程式入口
  * ════════════════════════════════════════════════════════════ */
 int main(int argc, char *argv[]) {
     /* ── 參數解析 ── */
     if (argc < 2) {
         print_usage(argv[0]);
         return EXIT_FAILURE;
     }
 
     const char *edge_file = argv[1];
 
     /* n_random 預設為 HETA_N_RANDOM_DEFAULT（1000），可由命令列覆寫 */
     int n_random = HETA_N_RANDOM_DEFAULT;
     if (argc >= 3) {
         n_random = atoi(argv[2]);
         if (n_random <= 0) {
             fprintf(stderr, "[main] n_random 必須為正整數，收到 '%s'\n", argv[2]);
             return EXIT_FAILURE;
         }
     }
 
     /* ── 建立 result 資料夾（若已存在則忽略）── */
     /* 0755：owner 可讀寫執行，group/other 可讀執行 */
 #if defined(_WIN32) || defined(_WIN64)
     mkdir("result");          /* Windows：mkdir 不需要權限參數 */
 #else
     mkdir("result", 0755);    /* macOS / Linux */
 #endif
 
     /* ── 讀取邊列表，建立圖 ── */
     fprintf(stdout, "[main] Loading edge file: %s\n", edge_file);
     Graph *g = load_edge_file(edge_file);
 
     fprintf(stdout, "[main] Graph: n=%d, m=%d\n", g->n, g->m);
 
     /* ── 計時開始（牆鐘時間，反映實際等待時間）── */
     struct timespec t_start, t_end;
     clock_gettime(CLOCK_MONOTONIC, &t_start);
 
     /* ── 執行 HETA ── */
     fprintf(stdout, "[main] Running HETA (n_random=%d)...\n", n_random);
     HetaResult *result = heta_run(g, n_random);
 
     /* ── 計時結束 ── */
     clock_gettime(CLOCK_MONOTONIC, &t_end);
     double elapsed = (double)(t_end.tv_sec  - t_start.tv_sec)
                    + (double)(t_end.tv_nsec - t_start.tv_nsec) / 1e9;
     fprintf(stdout, "[main] HETA completed in %.3f seconds\n", elapsed);
 
     /* ── 列印統計摘要到終端 ── */
     HetaSummary summary = heta_summarize(result);
     heta_print_summary(&summary, result->m);
 
     /* ── 建立輸出路徑並寫入 result 資料夾 ── */
     char out_path[512];
     build_output_path(edge_file, out_path, sizeof(out_path));
     write_result(out_path, g, result, &summary, elapsed, n_random);
 
     /* ── 釋放記憶體 ── */
     heta_result_free(result);
     graph_free(g);
 
     return EXIT_SUCCESS;
 }