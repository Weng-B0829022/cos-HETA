/*
 * rumor.h — 謠言擴散模型資料結構與函式宣告
 * =============================================
 * 對應 Python main.py 的完整動力傳播模型
 *
 * 包含：
 *   1. WS / NW 小世界網路生成
 *   2. 基於 HETA 邊類型的種子節點選取
 *   3. 信心度累積 + 衰減的謠言擴散模擬
 */

#ifndef RUMOR_H
#define RUMOR_H

#include "heta.h"

/* ────────────────────────────────────────────────────────────
 * 模型參數結構（對應 Python 全域參數）
 * ──────────────────────────────────────────────────────────── */
typedef struct {
    int    n_nodes;          /* 節點數                        */
    int    k_neighbors;      /* 每節點鄰居數（WS/NW 環形基底）*/
    double rewire_prob;      /* 重連 / 捷徑機率               */
    int    days;             /* 模擬天數                      */
    int    n_seeds;          /* 種子節點數                    */
    int    n_random_nets;    /* HETA 隨機化網路數量           */
    int    q_switch;         /* switching randomization Q 值  */
    int    random_seed;      /* 全域隨機種子                  */

    double decay_rate;       /* 每日無傳入信心度保留率        */
    int    decay_zero_day;   /* 連續無傳入幾天後歸零          */
    double willingness;      /* 全節點固定散播意願            */

    /* 各邊類型信心度貢獻率 */
    double conf_bond;
    double conf_local;
    double conf_global;
    double conf_silk;

    /* 感染節點加速衰減 */
    double infected_decay;   /* 感染節點(conf>=1)專用衰減率，比 decay_rate 更低 → 更快消亡 */
} RumorParams;

/* 預設參數（對應 Python main.py 設定）*/
static inline RumorParams rumor_default_params(void) {
    RumorParams p;

    /* ── 網路拓撲參數 ── */
    p.n_nodes        = 100;    /* 網路節點數量                                */
    p.k_neighbors    = 6;      /* WS/NW 環形基底中每節點的鄰居數（必須為偶數）*/
    p.rewire_prob    = 0.7;    /* WS 重連機率 / NW 捷徑新增比例              */

    /* ── 模擬控制參數 ── */
    p.days           = 60;     /* 模擬總天數                                  */
    p.n_seeds        = 5;      /* 初始感染種子節點數（從 global bridge 端點選取）*/
    p.n_random_nets  = 100;    /* HETA 演算法產生的隨機化網路數量（越大越精確但越慢）*/
    p.q_switch       = 100;    /* switching randomization 每條邊的交換次數    */
    p.random_seed    = 40;     /* 全域隨機種子（控制網路生成的可重現性）       */

    /* ── 信心度衰減參數（Ebbinghaus 遺忘曲線）── */
    p.decay_rate     = 0.75;   /* 每日無傳入時信心度保留率（0~1，越小遺忘越快）*/
    p.decay_zero_day = 5;     /* 連續無傳入超過此天數後信心度強制歸零        */
    p.willingness    = 0.1;   /* 全節點固定散播意願（0~1，越高越願意傳播）   */

    /* ── 各邊類型信心度貢獻率（每天感染節點對鄰居傳遞的信心度量）── */
    p.conf_bond      = 0.40;   /* bond（強連結）：傳得深但困在群組內          */
    p.conf_local     = 0.20;   /* local bridge（局部橋接）：中等傳遞力        */
    p.conf_global    = 0.70;   /* global bridge（全域橋接）：跨群但信任有限   */
    p.conf_silk      = 1.00;   /* silk（末梢節點）：無抵抗力，幾乎必然接受    */

    /* ── 感染節點加速衰減 ── */
    p.infected_decay = 0.0;    /* 感染節點專用衰減率（0=使用 decay_rate，<decay_rate 則更快消亡）*/

    return p;
}

/* ────────────────────────────────────────────────────────────
 * 模擬結果
 * ──────────────────────────────────────────────────────────── */
typedef struct {
    double *daily_ratio;     /* 每日感染比例陣列，大小 days+1 */
    int     n_days;          /* 陣列長度（= days + 1）        */
} SimResult;

/* ════════════════════════════════════════════════════════════
 * 函式宣告 — network.c
 * ════════════════════════════════════════════════════════════ */

/*
 * 建立 Watts-Strogatz 小世界網路
 * 參數：n 節點數, k 每側鄰居數, p 重連機率, seed 隨機種子
 * 回傳：CSR 格式 Graph*（呼叫者負責 graph_free）
 */
Graph *build_ws_network(int n, int k, double p, uint64_t seed);

/*
 * 建立 Newman-Watts 小世界網路（在環形上新增捷徑，不刪除原有邊）
 * 參數：n 節點數, k 每側鄰居數, p 新增捷徑比例, seed 隨機種子
 * 回傳：CSR 格式 Graph*
 */
Graph *build_nw_network(int n, int k, double p, uint64_t seed);

/* ════════════════════════════════════════════════════════════
 * 函式宣告 — rumor.c
 * ════════════════════════════════════════════════════════════ */

/*
 * 從 HETA 分類結果中選取所有 global bridge 端點作為種子節點
 * n_seeds 參數已不使用（保留介面相容），投放數量由 global bridge 數量決定
 * 回傳：動態分配的種子節點陣列（呼叫者負責 free）
 * out_count：實際選取的種子數量
 */
int *select_seeds_from_global(const Graph *g, const HetaResult *heta,
                              int n_seeds, int *out_count);

/*
 * 執行謠言擴散模擬
 * 回傳：SimResult（呼叫者負責 free(result.daily_ratio)）
 */
SimResult simulate_rumor(const Graph *g, const HetaResult *heta,
                         const int *seed_nodes, int n_seeds,
                         const RumorParams *params);

/*
 * 根據邊類型回傳對應的信心度貢獻率
 */
double get_confidence_rate(LinkType type, const RumorParams *params);

/*
 * 列印每 10 天摘要
 */
void print_daily_summary(const SimResult *ws, const SimResult *nw, int interval);

/*
 * 將結果寫入 CSV 檔
 */
void write_result_csv(const char *path,
                      const SimResult *ws, const SimResult *nw,
                      const RumorParams *params);

#endif /* RUMOR_H */
