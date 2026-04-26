# ════════════════════════════════════════════════════════════
# Makefile — HETA C 跨平台建置規則
# ════════════════════════════════════════════════════════════
#
# 支援平台：
#   - macOS   (Apple Silicon arm64 / Intel x86_64)
#   - Windows (透過 MSYS2/MinGW-w64 或 Cygwin)
#   - Linux   (x86_64 / aarch64)
#
# 使用方式：
#   make              - 自動偵測平台並編譯（含 OpenMP，若可用）
#   make no_omp       - 不使用 OpenMP（任何平台通用，最簡單）
#   make avx2         - 啟用 AVX2（x86 平台，Haswell/Zen 以上）
#   make debug        - 偵錯版本（停用優化）
#   make test         - 快速功能測試（karate club, n_random=50）
#   make bench        - 效能測試（n_random=1000）
#   make clean        - 清除編譯產物
#   make info         - 顯示目前偵測到的平台與編譯設定
#
# ────────────────────────────────────────────────────────────
# macOS 安裝步驟（OpenMP 需要額外安裝）：
#   方案 A（推薦，最簡單）：
#     brew install gcc
#     make
#
#   方案 B（使用系統 clang）：
#     brew install libomp
#     make
#
#   方案 C（不想裝任何東西，串列模式）：
#     make no_omp
#
# Windows 安裝步驟（建議使用 MSYS2）：
#   1. 安裝 MSYS2：https://www.msys2.org/
#   2. 開啟「MSYS2 MINGW64」終端，執行：
#      pacman -S mingw-w64-x86_64-gcc
#   3. 在同一終端執行：make
#
# Linux：
#   sudo apt install gcc    # Debian/Ubuntu
#   make
# ════════════════════════════════════════════════════════════

# ────────────────────────────────────────────────────────────
# 平台自動偵測
#   uname -s 回傳值：
#     Linux              → "Linux"
#     macOS              → "Darwin"
#     MSYS2/MinGW        → "MINGW64_NT-10.0-..." 等
#     Cygwin             → "CYGWIN_NT-10.0-..."
#   uname -m 回傳值：
#     x86_64, amd64      → 64-bit Intel/AMD
#     arm64, aarch64     → ARM 64-bit（Apple Silicon 等）
# ──────────────────────────────────────────────────────────────
UNAME_S := $(shell uname -s 2>/dev/null || echo Windows)
UNAME_M := $(shell uname -m 2>/dev/null || echo x86_64)

# ────────────────────────────────────────────────────────────
# 依平台設定執行檔副檔名、清除指令
# ──────────────────────────────────────────────────────────────

ifneq (,$(findstring MINGW,$(UNAME_S)))
    PLATFORM := Windows
    EXE      := .exe       # Windows 執行檔需要 .exe
    RM       := rm -rf
    MKDIR    := mkdir -p
else ifneq (,$(findstring CYGWIN,$(UNAME_S)))
    PLATFORM := Cygwin
    EXE      := .exe
    RM       := rm -rf
    MKDIR    := mkdir -p
else ifeq ($(UNAME_S),Darwin)
    PLATFORM := macOS
    EXE      :=
    RM       := rm -rf
    MKDIR    := mkdir -p
else
    PLATFORM := Linux
    EXE      :=
    RM       := rm -rf
    MKDIR    := mkdir -p
endif

# ────────────────────────────────────────────────────────────
# 編譯器選擇與 OpenMP 旗標（依平台不同）
# ──────────────────────────────────────────────────────────────

ifeq ($(PLATFORM),macOS)
    # ── macOS 策略：優先找 brew gcc，次選 clang + libomp ──

    # 嘗試找到 brew 安裝路徑
    BREW_PREFIX := $(shell brew --prefix 2>/dev/null || echo /opt/homebrew)

    # 找 brew gcc 的最新版本號（例：gcc-14、gcc-13）
    GCC_VER := $(shell ls $(BREW_PREFIX)/bin/gcc-* 2>/dev/null \
                | grep -oE 'gcc-[0-9]+' | sort -t- -k2 -rn | head -1 \
                | grep -oE '[0-9]+')

    ifneq ($(GCC_VER),)
        # ── 找到 brew gcc → 使用它（OpenMP 開箱即用）──
        CC          := $(BREW_PREFIX)/bin/gcc-$(GCC_VER)
        OMP_CFLAGS  := -fopenmp
        OMP_LDFLAGS := -fopenmp
        HAS_OMP     := yes
        CC_DESC     := brew gcc-$(GCC_VER) (OpenMP built-in)
    else
        # brew gcc 不存在，改用系統 clang
        CC := clang

        # 檢查 brew libomp 是否安裝
        LIBOMP_DIR := $(shell brew --prefix libomp 2>/dev/null)
        LIBOMP_OK  := $(shell test -f "$(LIBOMP_DIR)/lib/libomp.dylib" \
                       2>/dev/null && echo yes || echo no)

        ifeq ($(LIBOMP_OK),yes)
            # ── clang + brew libomp ──
            # clang 需要 -Xpreprocessor 才能識別 -fopenmp
            OMP_CFLAGS  := -Xpreprocessor -fopenmp -I$(LIBOMP_DIR)/include
            OMP_LDFLAGS := -L$(LIBOMP_DIR)/lib -lomp
            HAS_OMP     := yes
            CC_DESC     := Apple clang + brew libomp
        else
            # ── 純 clang，無 OpenMP（功能完整，只有 threshold 較慢）──
            OMP_CFLAGS  :=
            OMP_LDFLAGS :=
            HAS_OMP     := no
            CC_DESC     := Apple clang (no OpenMP)
        endif
    endif

    # macOS 數學函式庫內建，不需要 -lm
    MATH_LIB := 

    # 依 CPU 架構選擇優化旗標
    ifeq ($(UNAME_M),arm64)
        # Apple Silicon：-mcpu=apple-m1 啟用 NEON 向量化
        # 不支援 -mavx2，所以不加 -march=native
        ARCH_FLAGS := -mcpu=apple-m1
    else
        # Intel Mac：-march=native 啟用 SSE4.2 + POPCNT
        ARCH_FLAGS := -march=native
    endif

else ifeq ($(PLATFORM),Windows)
    # ── Windows MSYS2 MinGW：gcc 含完整 OpenMP ──
    CC          := gcc
    OMP_CFLAGS  := -fopenmp
    OMP_LDFLAGS := -fopenmp
    HAS_OMP     := yes
    CC_DESC     := MinGW gcc (OpenMP built-in)
    MATH_LIB    := -lm
    ARCH_FLAGS  := -march=native

else ifeq ($(PLATFORM),Cygwin)
    # ── Cygwin：gcc 含 OpenMP ──
    CC          := gcc
    OMP_CFLAGS  := -fopenmp
    OMP_LDFLAGS := -fopenmp
    HAS_OMP     := yes
    CC_DESC     := Cygwin gcc
    MATH_LIB    := -lm
    ARCH_FLAGS  := -march=native

else
    # ── Linux：gcc + libgomp（OpenMP 標準支援）──
    CC          := gcc
    OMP_CFLAGS  := -fopenmp
    OMP_LDFLAGS := -fopenmp
    HAS_OMP     := yes
    CC_DESC     := gcc (OpenMP built-in)
    MATH_LIB    := -lm
    ARCH_FLAGS  := -march=native
endif

# ────────────────────────────────────────────────────────────
# 目錄與目標檔名
# ──────────────────────────────────────────────────────────────
SRCDIR := src
INCDIR := include
OBJDIR := obj
TARGET := heta$(EXE)
RUMOR_TARGET := rumor_spread$(EXE)

# ────────────────────────────────────────────────────────────
# 編譯旗標組合
#
#   -std=c11      : C11 標準（所有平台一致）
#   -Wall -Wextra : 完整警告（捕捉潛在 bug）
#   -O3           : 最高優化（迴圈展開、函式內聯、向量化）
#   -funroll-loops: 手動提示展開（bitset 迴圈效果顯著）
#   -ffast-math   : 允許浮點重排序（Welford 加速 ~10%）
#                   注意：若需嚴格 IEEE 754 相容，請移除
#   $(ARCH_FLAGS) : 依 CPU 架構（native / apple-m1）
#   $(OMP_CFLAGS) : OpenMP（依平台差異由上方自動設定）
# ──────────────────────────────────────────────────────────────
BASE_CFLAGS := -std=c11 -Wall -Wextra -O3 \
               -funroll-loops -ffast-math \
               $(ARCH_FLAGS) \
               -I$(INCDIR)

CFLAGS  := $(BASE_CFLAGS) $(OMP_CFLAGS)
LDFLAGS := $(OMP_LDFLAGS) $(MATH_LIB)

# ────────────────────────────────────────────────────────────
# 原始檔列表
# ──────────────────────────────────────────────────────────────
SRCS := $(SRCDIR)/graph.c      \
        $(SRCDIR)/bitset.c     \
        $(SRCDIR)/cnr.c        \
        $(SRCDIR)/randomize.c  \
        $(SRCDIR)/threshold.c  \
        $(SRCDIR)/heta.c       \
        $(SRCDIR)/main.c

OBJS := $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCS))

# 謠言擴散模型的原始檔（共用 HETA 核心 + 新增 network.c, rumor.c, rumor_main.c）
RUMOR_SRCS := $(SRCDIR)/graph.c      \
              $(SRCDIR)/bitset.c     \
              $(SRCDIR)/cnr.c        \
              $(SRCDIR)/randomize.c  \
              $(SRCDIR)/threshold.c  \
              $(SRCDIR)/heta.c       \
              $(SRCDIR)/network.c    \
              $(SRCDIR)/rumor.c      \
              $(SRCDIR)/rumor_main.c

RUMOR_OBJS := $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(RUMOR_SRCS))

# ════════════════════════════════════════════════════════════
# 主要建置目標
# ════════════════════════════════════════════════════════════

.PHONY: all
all: clean _show_platform $(OBJDIR) $(TARGET) $(RUMOR_TARGET)
	@echo ""
	@echo "  建置完成：$(TARGET) 與 $(RUMOR_TARGET)"
	@echo "  HETA 使用方式：./$(TARGET) <edge_file> [n_random]"
	@echo "  謠言模型使用方式：./$(RUMOR_TARGET) [n_random] [n_nodes] [days]"

# 建置前顯示平台資訊
.PHONY: _show_platform
_show_platform:
	@echo "  平台：$(PLATFORM) / $(UNAME_M) / $(CC_DESC) / OpenMP=$(HAS_OMP)"

# 建立 obj 目錄
$(OBJDIR):
	$(MKDIR) $(OBJDIR)

# 連結 HETA 主程式
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 連結謠言擴散模型主程式
$(RUMOR_TARGET): $(RUMOR_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 編譯規則：.c → .o
# $< = 第一個依賴（.c 檔）；$@ = 目標（.o 檔）
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(INCDIR)/heta.h
	$(CC) $(CFLAGS) -c $< -o $@

# ════════════════════════════════════════════════════════════
# no_omp — 強制停用 OpenMP（最通用，任何平台均可編譯）
#
# 適合：
#   - macOS 不想裝 brew（make no_omp 直接用系統 clang）
#   - 排查 OpenMP 相關問題
#   - CI 環境不確定是否有 OpenMP
# ════════════════════════════════════════════════════════════
.PHONY: no_omp
no_omp: CFLAGS  := $(BASE_CFLAGS)
no_omp: LDFLAGS := $(MATH_LIB)
no_omp: clean $(OBJDIR) $(TARGET)
	@echo ""
	@echo "  串列模式（無 OpenMP）"
	@echo "  threshold 計算為單核，其餘優化不變"

# ════════════════════════════════════════════════════════════
# avx2 — 啟用 AVX2（僅 x86_64）
# Apple Silicon 請勿使用
# ════════════════════════════════════════════════════════════
.PHONY: avx2
avx2:
ifeq ($(UNAME_M),arm64)
	@echo "  錯誤：AVX2 不支援 Apple Silicon (arm64)"
	@echo "  arm64 已透過 -mcpu=apple-m1 自動啟用 NEON 向量化"
	@echo "  請直接使用 make 或 make no_omp"
	@exit 1
else
	$(MAKE) CFLAGS="$(BASE_CFLAGS) $(OMP_CFLAGS) -mavx2" \
	        LDFLAGS="$(LDFLAGS)"
	@echo "  AVX2 版本建置完成"
endif

# ════════════════════════════════════════════════════════════
# debug — 偵錯版本（停用優化與 OpenMP）
# Windows MinGW 的 ASAN 支援有限，自動降級為無 ASAN 版本
# ════════════════════════════════════════════════════════════
.PHONY: debug
ifeq ($(PLATFORM),Windows)
debug: CFLAGS  := -std=c11 -Wall -Wextra -g -O0 -I$(INCDIR)
debug: LDFLAGS := $(MATH_LIB)
debug: clean $(OBJDIR) $(TARGET)
	@echo "  偵錯版本（Windows，無 ASAN）"
else
debug: CFLAGS  := -std=c11 -Wall -Wextra -g -O0 \
                  -fsanitize=address,undefined \
                  -fno-omit-frame-pointer -I$(INCDIR)
debug: LDFLAGS := -fsanitize=address,undefined $(MATH_LIB)
debug: clean $(OBJDIR) $(TARGET)
	@echo "  偵錯版本（含 AddressSanitizer + UBSan）"
endif

# ════════════════════════════════════════════════════════════
# test — 快速功能測試（karate club, n_random=50）
# ════════════════════════════════════════════════════════════
.PHONY: test
test: all karate_test.txt
	@echo "── 執行 HETA 測試（karate club, n_random=50）──"
	./$(TARGET) karate_test.txt 50

# 生成 Zachary's Karate Club 邊列表（34 節點、78 條邊）
# 使用 printf（比 echo 跨平台相容性更好）
karate_test.txt:
	@printf "# Zachary's Karate Club Network (34 nodes, 78 edges)\n" > $@
	@printf "# Source: Zachary (1977), J. Anthropol. Res. 33\n" >> $@
	@printf "0 1\n0 2\n0 3\n0 4\n0 5\n0 6\n0 7\n0 8\n0 10\n0 11\n" >> $@
	@printf "0 12\n0 13\n0 17\n0 19\n0 21\n0 31\n1 2\n1 3\n1 7\n1 13\n" >> $@
	@printf "1 17\n1 19\n1 21\n1 30\n2 3\n2 7\n2 8\n2 9\n2 13\n2 27\n" >> $@
	@printf "2 28\n2 32\n3 7\n3 12\n3 13\n4 6\n4 10\n5 6\n5 10\n5 16\n" >> $@
	@printf "6 16\n8 30\n8 32\n8 33\n9 33\n13 33\n14 32\n14 33\n15 32\n" >> $@
	@printf "15 33\n18 32\n18 33\n19 33\n20 32\n20 33\n22 32\n22 33\n" >> $@
	@printf "23 25\n23 27\n23 29\n23 32\n23 33\n24 25\n24 27\n24 31\n" >> $@
	@printf "25 31\n26 29\n26 33\n27 33\n28 31\n28 33\n29 32\n29 33\n" >> $@
	@printf "30 32\n30 33\n31 32\n31 33\n32 33\n" >> $@

# ════════════════════════════════════════════════════════════
# bench — 效能基準測試
# ════════════════════════════════════════════════════════════
.PHONY: bench
bench: all karate_test.txt
	./$(TARGET) karate_test.txt 1000

# ════════════════════════════════════════════════════════════
# clean — 清除編譯產物
# ════════════════════════════════════════════════════════════
# ════════════════════════════════════════════════════════════
# rumor — 只建置謠言擴散模型
# ════════════════════════════════════════════════════════════
.PHONY: rumor
rumor: clean _show_platform $(OBJDIR) $(RUMOR_TARGET)
	@echo ""
	@echo "  謠言擴散模型建置完成：$(RUMOR_TARGET)"
	@echo "  使用方式：./$(RUMOR_TARGET) [n_random] [n_nodes] [days]"

# ════════════════════════════════════════════════════════════
# rumor_test — 快速測試謠言擴散模型
# ════════════════════════════════════════════════════════════
.PHONY: rumor_test
rumor_test: rumor
	@echo "── 執行謠言擴散模型測試（n_random=50）──"
	./$(RUMOR_TARGET) 50

.PHONY: clean
clean:
	$(RM) $(OBJDIR) $(TARGET) $(RUMOR_TARGET)
	@echo "  清除完成"

# ════════════════════════════════════════════════════════════
# info — 顯示平台偵測結果
# ════════════════════════════════════════════════════════════
.PHONY: info
info:
	@echo "════════════════════════════════"
	@echo " 平台偵測結果"
	@echo "════════════════════════════════"
	@echo "  OS        : $(PLATFORM)"
	@echo "  CPU arch  : $(UNAME_M)"
	@echo "  編譯器    : $(CC)"
	@echo "  描述      : $(CC_DESC)"
	@echo "  OpenMP    : $(HAS_OMP)"
	@echo "  執行檔    : $(TARGET)"
	@echo "  ARCH_FLAGS: $(ARCH_FLAGS)"
	@echo "  CFLAGS    : $(CFLAGS)"
	@echo "  LDFLAGS   : $(LDFLAGS)"
	@echo "════════════════════════════════"
ifeq ($(PLATFORM),macOS)
ifeq ($(HAS_OMP),no)
	@echo ""
	@echo "  提示：目前無 OpenMP，以下任一方式可啟用並行加速："
	@echo "    brew install gcc      → make clean && make"
	@echo "    brew install libomp   → make clean && make"
endif
endif

# ════════════════════════════════════════════════════════════
# auto — 批次自動執行（掃 data/net/ 全部 .net 檔）
#
# 與 heta 差異：
#   - heta 是單檔模式：./heta <edge_file> [n_random]
#   - auto 是批次模式：./auto（無參數，掃整個 data/net/）
#   - auto 內建 Pajek 讀取，可直接吃 .net 檔
#
# 使用方式：
#   make auto       - 只編譯
#   make auto_run   - 編譯後自動執行
# ════════════════════════════════════════════════════════════

AUTO_TARGET := auto$(EXE)

AUTO_SRCS := $(SRCDIR)/graph.c      \
             $(SRCDIR)/bitset.c     \
             $(SRCDIR)/cnr.c        \
             $(SRCDIR)/randomize.c  \
             $(SRCDIR)/threshold.c  \
             $(SRCDIR)/heta.c       \
             $(SRCDIR)/auto.c

AUTO_OBJS := $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(AUTO_SRCS))

.PHONY: auto
auto: _show_platform $(OBJDIR) $(AUTO_TARGET)
	@echo ""
	@echo "  auto 建置完成：$(AUTO_TARGET)"
	@echo "  使用方式：./$(AUTO_TARGET)（請在專案根目錄執行）"

# 連結 auto
$(AUTO_TARGET): $(AUTO_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: auto_run
auto_run: auto
	@echo ""
	@echo "── 執行 auto（n_random=1000，掃 data/net/ 全部 .net 檔）──"
	./$(AUTO_TARGET)
