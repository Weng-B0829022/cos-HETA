# HETA — Hierarchical Edge Type Analysis（C 高效能實作）

論文來源：Huang et al. (2019), Physica A 536, 121027  
"Beyond bond links in complex networks: Local bridges, global bridges and silk links"

---

## 檔案結構

```
heta_c/
├── include/
│   └── heta.h          # 所有資料結構與函式宣告
├── src/
│   ├── graph.c         # CSR 圖結構、BFS、最短路徑
│   ├── bitset.c        # Bitset 操作（含 AVX2 SIMD 支援）
│   ├── cnr.c           # k 層鄰居 BFS 與共同鄰居比例計算
│   ├── randomize.c     # Switching randomization（PCG32 RNG）
│   ├── threshold.c     # 外部/內部閾值計算（OpenMP + Welford）
│   ├── heta.c          # 主演算法流程（三步分類）
│   └── main.c          # 程式入口、邊表讀取、CSV 輸出
└── Makefile
```

---

## 編譯

```bash
# 標準編譯（含 OpenMP 並行）
make

# AVX2 加速版（Haswell / Zen 以上 CPU）
make avx2

# 偵錯版（AddressSanitizer + UBSan）
make debug

# 快速功能測試（karate club, n_random=50）
make test
```

**依賴**：GCC 7+（或 Clang 6+）、libm、OpenMP

---

## 使用方式

```bash
./heta <edge_file> [n_random]
```

- `edge_file`：純文字邊列表，每行 `u v`，`#` 開頭為注釋
- `n_random`：隨機化網路數量，預設 1000（論文標準值），測試用 100

**範例：**
```bash
./heta karate_test.txt 1000
./heta my_network.txt 100
```

---

## 輸入格式

```
# 注釋行（# 開頭）
0 1
0 2
1 3
...
```

- 節點 id 從 0 開始（0-indexed）
- 無向圖：每條邊只需輸入一次
- 自動去除自迴圈（u == v 的邊）

---

## 輸出格式

程式輸出兩部分：

1. **統計摘要**（到 stderr）：
```
════════════════════════════════
 HETA 分類結果摘要（共 78 條邊）
════════════════════════════════
  silk link      :    0 條  (  0.0%)
  bond link      :   45 條  ( 57.7%)
  local bridge   :   18 條  ( 23.1%)
  global bridge  :   15 條  ( 19.2%)
════════════════════════════════
```

2. **逐邊分類結果**（CSV，到 stdout）：
```csv
u,v,type,layer
0,1,bond,0
0,2,bond,0
0,5,local_bridge_k1,1
...
```

可重定向 stdout 儲存 CSV：
```bash
./heta karate_test.txt 1000 > results.csv
```

---

## 優化技術說明

| 模組 | Python 原版瓶頸 | C 優化方式 | 預期加速 |
|------|----------------|-----------|---------|
| 圖複製（×1000） | dict deepcopy | memcpy CSR 陣列 | ~100× |
| has_edge 查詢 | dict lookup | bitset 位元測試 O(1) | ~20× |
| 鄰居集合交集 | Python set & | bitwise AND + POPCNT | ~64× |
| BFS visited | Python set | bool 陣列 + memset | ~10× |
| 1000 隨機網路 | 串列 | OpenMP 多核並行 | ~N_core× |
| 隨機數生成 | Python random | PCG32 演算法 | ~5× |
| mean/std 計算 | 先存陣列再 np | Welford 線上演算法 | 記憶體 O(1) |
| **綜合估計** | 基準 | | **~500×~2000×** |

---

## 演算法對應

| C 函式 | 對應論文 |
|--------|--------|
| `graph_kmax()` | 公式 (4)：kmax = floor(avg_spl / 2) |
| `compute_cnr()` k=1 | 公式 (1)：R^1_{u,v} |
| `compute_cnr()` k>1 | 公式 (2)：R^k_{u,v} |
| `compute_kth_layer()` | 公式 (3)：V^{k,j}_i |
| `switching_randomize()` | Section 2：Switching Algorithm |
| `compute_external_threshold()` | 公式 (5)：T_E^k |
| `compute_internal_threshold()` | 公式 (6)：T_I^k |
| `heta_run()` Step 3.1 | Algorithm 1 Step 3.1：silk links |
| `heta_run()` Step 3.2 | Algorithm 1 Step 3.2：bond / local bridge |
| `heta_run()` Step 3.3 | Algorithm 1 Step 3.3：global bridge |
