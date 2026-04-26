/*
 * heta.h — HETA 演算法核心資料結構與函式宣告
 * ================================================
 * 論文來源：Huang et al. (2019), Physica A 536, 121027
 * "Beyond bond links in complex networks:
 *  Local bridges, global bridges and silk links"
 *
 * 設計原則：
 *   1. 圖用 CSR（Compressed Sparse Row）格式儲存，確保 cache 友善
 *   2. 鄰居集合用 bitset（uint64_t 陣列）表示，交集變成 bitwise AND
 *   3. BFS 使用靜態陣列 queue + bool visited，避免動態記憶體分配
 *   4. 1000 個隨機網路的計算透過 OpenMP 並行處理
 */

#ifndef HETA_H
#define HETA_H

#include <stdint.h>   /* uint64_t, uint32_t */
#include <stddef.h>   /* size_t             */

/* ────────────────────────────────────────────────────────────
 * 編譯期常數
 * ──────────────────────────────────────────────────────────── */

/* 論文建議的隨機化網路數量（預設 1000，測試時可縮小） */
#define HETA_N_RANDOM_DEFAULT  1000

/* switching randomization 每條邊的 swap 次數（論文建議 Q=100） */
#define HETA_SWITCH_Q          100

/* bitset 每個 word 的位元數 */
#define BS_WORD_BITS           64u

/* 計算存放 n 個節點所需的 word 數量（無號，向上取整） */
#define BS_WORDS(n)            (((n) + BS_WORD_BITS - 1u) / BS_WORD_BITS)

/* ────────────────────────────────────────────────────────────
 * 邊類型枚舉
 * ──────────────────────────────────────────────────────────── */
typedef enum {
    LINK_UNCLASSIFIED  = 0,  /* 尚未分類（初始狀態）                    */
    LINK_SILK          = 1,  /* silk link：端點之一 degree == 1          */
    LINK_BOND          = 2,  /* bond link：R^k >= T_E^k                  */
    LINK_LOCAL_BRIDGE  = 3,  /* kth-layer local bridge：T_I^k < R^k < T_E^k */
    LINK_GLOBAL_BRIDGE = 4   /* global bridge：所有層皆未能分類          */
} LinkType;

/* ────────────────────────────────────────────────────────────
 * 邊的完整資訊
 * ──────────────────────────────────────────────────────────── */
typedef struct {
    int      u;          /* 邊的起點節點 id（0-indexed）                */
    int      v;          /* 邊的終點節點 id（0-indexed）                */
    LinkType type;       /* 最終分類結果                                 */
    int      layer;      /* local bridge 所在層（其餘類型為 0）          */
    int      active;     /* 1 = 尚未分類（pass_flag），0 = 已分類        */
} Edge;

/* ────────────────────────────────────────────────────────────
 * CSR 格式圖結構
 *
 * 鄰接表佈局（以 4 節點、5 邊的圖為例）：
 *   adj_list  = [1,2, 0,2,3, 0,1, 1]   ← 所有節點的鄰居連續存放
 *   offset[0] = 0   → 節點 0 的鄰居是 adj_list[0..1] = {1,2}
 *   offset[1] = 2   → 節點 1 的鄰居是 adj_list[2..4] = {0,2,3}
 *   offset[2] = 5   → 節點 2 的鄰居是 adj_list[5..6] = {0,1}
 *   offset[3] = 7   → 節點 3 的鄰居是 adj_list[7..7] = {1}
 *   offset[4] = 8   ← 哨兵：節點 n-1 結束位置
 *
 * 優點：連續記憶體存取，cache miss 極少
 * ──────────────────────────────────────────────────────────── */
typedef struct {
    int  n;          /* 節點數                                           */
    int  m;          /* 邊數                                             */
    int *offset;     /* 大小 n+1，offset[u]..offset[u+1]-1 是節點 u 的鄰居範圍 */
    int *adj_list;   /* 大小 2m（無向圖每條邊存兩次）                    */
    int *degree;     /* 大小 n，degree[u] = 節點 u 的度數                */
} Graph;

/* ────────────────────────────────────────────────────────────
 * Bitset 結構
 *
 * 用途：表示一個節點的「kth 層鄰居集合」
 *       兩個 bitset 的交集 = 逐 word AND，比 hash set 快 ~64 倍
 *
 * 使用方式：
 *   Bitset *bs = bs_alloc(n);   // 建立可容納 n 個節點的 bitset
 *   bs_set(bs, node_id);        // 將 node_id 加入集合
 *   int cnt = bs_intersect_count(a, b, n_words); // |A ∩ B|
 * ──────────────────────────────────────────────────────────── */
typedef struct {
    uint64_t *words;   /* bitset 資料陣列，words[i] 的第 j 位代表節點 i*64+j */
    int       n_words; /* 陣列長度 = ceil(n / 64)                         */
} Bitset;

/* ────────────────────────────────────────────────────────────
 * BFS 工作緩衝區（thread-local，避免 malloc）
 *
 * 在並行區段內每個 thread 持有一份，確保 thread safety
 * ──────────────────────────────────────────────────────────── */
typedef struct {
    int  *queue;     /* BFS 隊列，大小 n                                 */
    int  *visited;   /* visited 旗標陣列，大小 n（0/1）                  */
    int  *layer_buf; /* 暫存每個節點所屬層號（-1 = 未訪問），大小 n      */
} BFSBuffer;

/* ────────────────────────────────────────────────────────────
 * HETA 主結果結構
 * ──────────────────────────────────────────────────────────── */
typedef struct {
    Edge  *edges;     /* 邊陣列，大小 m                                   */
    int    m;         /* 邊數                                             */
    int    kmax;      /* 最大層數 = floor(avg_spl / 2)                    */
    double *ext_th;   /* 外部閾值陣列，ext_th[k]（k 從 1 開始，大小 kmax+1）*/
} HetaResult;

/* ────────────────────────────────────────────────────────────
 * 統計摘要
 * ──────────────────────────────────────────────────────────── */
typedef struct {
    int    n_silk;         /* silk link 數量   */
    int    n_bond;         /* bond link 數量   */
    int    n_local;        /* local bridge 總數 */
    int    n_global;       /* global bridge 數量 */
    double pct_silk;       /* silk link 百分比  */
    double pct_bond;       /* bond link 百分比  */
    double pct_local;      /* local bridge 百分比 */
    double pct_global;     /* global bridge 百分比 */
} HetaSummary;

/* ════════════════════════════════════════════════════════════
 * 函式宣告 — graph.c
 * ════════════════════════════════════════════════════════════ */

/* 從邊列表建立 CSR 圖；呼叫者負責之後呼叫 graph_free() */
Graph *graph_build(int n, int m, const int *us, const int *vs);

/* 釋放 CSR 圖佔用的所有記憶體 */
void   graph_free(Graph *g);

/* 複製圖（用於 switching randomization）*/
Graph *graph_copy(const Graph *g);

/* 計算圖的平均最短路徑長度（BFS，不連通圖取最大連通子圖）*/
double graph_avg_spl(const Graph *g);

/* 計算 kmax = floor(avg_spl / 2)，最小值 1 */
int    graph_kmax(const Graph *g);

/* ════════════════════════════════════════════════════════════
 * 函式宣告 — bitset.c
 * ════════════════════════════════════════════════════════════ */

/* 分配可容納 n 個節點的 bitset（全 0 初始化）*/
Bitset *bs_alloc(int n);

/* 釋放 bitset */
void    bs_free(Bitset *bs);

/* 將 bitset 全部清零 */
void    bs_clear(Bitset *bs);

/* 複製 src → dst（兩者 n_words 必須相同）*/
void    bs_copy(Bitset *dst, const Bitset *src);

/* 設定第 node 位元 */
static inline void bs_set(Bitset *bs, int node) {
    bs->words[node / BS_WORD_BITS] |= (1ULL << (node % BS_WORD_BITS));
}

/* 測試第 node 位元是否已設定 */
static inline int bs_test(const Bitset *bs, int node) {
    return (bs->words[node / BS_WORD_BITS] >> (node % BS_WORD_BITS)) & 1;
}

/* 計算兩個 bitset 的交集大小 |A ∩ B| */
int bs_intersect_count(const Bitset *a, const Bitset *b);

/* 計算兩個 bitset 交集的 popcount，並將結果存入 out（out = A AND B）*/
int bs_intersect_into(const Bitset *a, const Bitset *b, Bitset *out);

/* 計算 bitset 中被設定的位元總數（popcount）*/
int bs_popcount(const Bitset *bs);

/* ════════════════════════════════════════════════════════════
 * 函式宣告 — cnr.c
 * ════════════════════════════════════════════════════════════ */

/* 分配 BFS 緩衝區（每個 thread 呼叫一次）*/
BFSBuffer *bfs_buffer_alloc(int n);

/* 釋放 BFS 緩衝區 */
void       bfs_buffer_free(BFSBuffer *buf);

/*
 * 計算節點 u 排除 exclude 後的第 k 層鄰居，結果存入 bitset out
 * buf 是 thread-local 的 BFS 緩衝區
 */
void compute_kth_layer(const Graph *g, int u, int exclude, int k,
                       Bitset *out, BFSBuffer *buf);

/*
 * 計算邊 (u,v) 的第 k 層共同鄰居比例 R^k_{u,v}
 * 論文公式 (1)（k=1）和公式 (2)（k>1）
 * buf 是 thread-local BFS 緩衝區
 * bs_k_u, bs_k_v, bs_k1_u, bs_k1_v 是預分配的 bitset 工作空間
 */
double compute_cnr(const Graph *g, int u, int v, int k,
                   BFSBuffer *buf,
                   Bitset *bs_k_u, Bitset *bs_k_v,
                   Bitset *bs_k1_u, Bitset *bs_k1_v);

/* ════════════════════════════════════════════════════════════
 * 函式宣告 — randomize.c
 * ════════════════════════════════════════════════════════════ */

/*
 * 對圖 g 執行一次 switching randomization（原地修改 gr）
 * gr 必須是 g 的複本（graph_copy 取得）
 * seed 是此次隨機化的 RNG 種子（每個 thread 不同）
 */
void switching_randomize(Graph *gr, uint64_t seed, Bitset **adj_bs);
/* 分配 n 個空 Bitset 的指標陣列（預分配用，不填入鄰居資料）*/
Bitset **adj_bs_alloc(int n);
/* 釋放 adj_bs_alloc 分配的陣列 */
void     adj_bs_free(Bitset **adj_bs, int n);
/* ════════════════════════════════════════════════════════════
 * 函式宣告 — threshold.c
 * ════════════════════════════════════════════════════════════ */

/*
 * 計算第 k 層的外部閾值 T_E^k
 * 論文公式 (5)：T_E^k = Mean_E^k(RG) + 2 * SD_E^k(RG)
 * n_random：隨機化網路數量（論文建議 1000）
 * 使用 OpenMP 並行生成 n_random 個隨機網路
 */
double compute_external_threshold(const Graph *g, int k, int n_random);

/*
 * 計算內部閾值 T_I^k
 * 論文公式 (6)：T_I^k = Mean_I^k - SD_I^k
 * ratios：候選邊的 R^k 值陣列；n：候選邊數量
 */
double compute_internal_threshold(const double *ratios, int n);

/* ════════════════════════════════════════════════════════════
 * 函式宣告 — heta.c
 * ════════════════════════════════════════════════════════════ */

/*
 * 主演算法入口：對圖 g 執行完整 HETA
 * n_random：隨機化網路數量（論文建議 1000）
 * 回傳：HetaResult*，呼叫者負責呼叫 heta_result_free()
 */
HetaResult *heta_run(const Graph *g, int n_random);

/* 釋放 HetaResult */
void heta_result_free(HetaResult *r);

/* 統計各類型數量與比例 */
HetaSummary heta_summarize(const HetaResult *r);

/* 列印統計摘要到 stdout */
void heta_print_summary(const HetaSummary *s, int total_edges);

#endif /* HETA_H */
