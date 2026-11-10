# 算法竞赛数学教程（C++ 版）

> 面向 OI / ICPC / CCPC 选手的系统性数学讲义。
> 全部代码为 **C++17**，可直接编译（`g++ -std=c++17 -O2`）。
> 目标：不仅能背板子，更理解"为什么这样做"，以及"什么时候该用哪块板子"。

---

## 这份教程适合谁

- 已掌握 C++ 基础语法、递归、STL 容器与常用算法（排序、二分、DP、图论基础）；
- 想系统补齐算法竞赛中**数学板块**的选手；
- 准备 NOI / ICPC 区域赛 / Codeforces Div.1 难度的读者。

数学板块在竞赛中的定位：**它是"工具"而不是"知识点"**。
竞赛里几乎没有纯粹考数学题，数学是拿来给其他问题降复杂度的：
把 $O(n)$ 的暴力变成 $O(\sqrt n)$ 的整除分块，把 $O(n)$ 的 DP 变成 $O(\log n)$ 的矩阵快速幂，
把 $O(n^2)$ 的卷积变成 $O(n\log n)$ 的 NTT……

所以本教程每一章都遵循同一结构：

1. **直觉与动机**：为什么需要这个东西；
2. **定义与定理**：严谨表述，带必要的证明或证明思路；
3. **算法与复杂度**：做法、正确性论证、复杂度分析；
4. **代码模板**：工业级实现，带注释；
5. **例题与套路**：典型题目怎么把问题化为本章工具；
6. **易错点清单**：竞赛中真实踩过的坑。

---

## 目录

| 章节 | 内容 | 关键工具 | 建议难度 |
|---|---|---|---|
| [第 1 章 数论基础](01-数论基础.md) | 整除、gcd、exgcd、快速幂、素数、筛法、质因数分解 | exgcd / 线性筛 / Pollard-Rho | ★★ |
| [第 2 章 同余理论](02-同余理论.md) | 逆元、线性同余方程、CRT/exCRT、欧拉函数、原根、BSGS、二次剩余 | exCRT / BSGS / Cipolla | ★★★ |
| [第 3 章 积性函数与狄利克雷卷积](03-积性函数与狄利克雷卷积.md) | 卷积、莫比乌斯反演、整除分块、杜教筛、min_25 筛 | 杜教筛 / 莫反套路 | ★★★★ |
| [第 4 章 组合数学](04-组合数学.md) | 二项式、容斥、卡特兰/斯特林/贝尔、Lucas/exLucas、生成函数 | Lucas / OGF / EGF | ★★★ |
| [第 5 章 概率与期望](05-概率与期望.md) | 条件概率、期望线性性、期望 DP、高斯消元解期望、随机游走 | 期望线性性 / 期望 DP | ★★★ |
| [第 6 章 线性代数](06-线性代数.md) | 矩阵快速幂、高斯消元、线性基、行列式、矩阵树定理 | 矩阵快速幂 / 线性基 | ★★★ |
| [第 7 章 多项式与快速变换](07-多项式与快速变换.md) | FFT、NTT、多项式求逆/ln/exp/sqrt、分治 FFT、拉格朗日插值、线性递推 | NTT / 多项式全家桶 | ★★★★★ |
| [第 8 章 博弈论](08-博弈论.md) | Nim、SG 函数、反 Nim、阶梯博弈、翻硬币、威佐夫 | SG / Nim 变形 | ★★★ |
| [第 9 章 计算几何](09-计算几何.md) | 向量、凸包、旋转卡壳、半平面交、最小圆覆盖、Pick 定理 | 凸包 / 半平面交 | ★★★★ |
| [第 10 章 杂项与进阶技巧](10-杂项与进阶技巧.md) | Burnside/Pólya、单纯形、三分/模拟退火、递推与特征方程 | Pólya / 单纯形 | ★★★★ |
| [附录 A 常用数学公式速查](附录A-公式速查.md) | 各种求和式、恒等式、常数、量级估算 | —— | —— |

---

## 学习路线建议

**第一阶段（打地基，约 2 周）**

第 1 章全部 + 第 2 章的逆元 / CRT / 欧拉函数 + 第 4 章的排列组合与容斥 + 第 6 章的矩阵快速幂。

这一阶段的目标是形成"条件反射"：看到 `1e9+7` 就想到逆元，看到 $n \le 10^{18}$ 就想到 $\sqrt n$ 或 $\log n$，看到重复线性转移就想到矩阵快速幂。

**第二阶段（补主流，约 3 周）**

第 3 章莫比乌斯反演与整除分块 + 第 4 章生成函数入门 + 第 5 章期望 DP + 第 6 章高斯消元与线性基 + 第 8 章博弈论。

莫反 + 整除分块是**出现频率最高**的进阶数学套路，务必掌握。

**第三阶段（冲高分，约 4 周）**

第 7 章多项式全家桶 + 第 9 章计算几何 + 第 3 章杜教筛 + 第 10 章。

多项式题在 ICPC 区域赛和 CF Div.1 中属于"拉开档次"的题；计算几何则是最考验代码细节的板块。

---

## 代码验证

本教程中**每一个 C++ 代码块都经过编译检查**，核心算法还**与暴力实现做了对拍验证**（共 71 项）。

| 工具 | 作用 | 结果 |
|---|---|---|
| `tools/check.ps1` | 抽取全部 Markdown 中的 C++ 代码块并逐个编译检查 | 103 个代码块，0 失败 |
| `tools/verify_nt.cpp` | 数论模板对拍（exgcd / 筛法 / Miller-Rabin / Pollard-Rho / exCRT / BSGS / Cipolla / 杜教筛 …） | 27 项全部通过 |
| `tools/verify_more.cpp` | 多项式与线性代数对拍（NTT / 多项式全家桶 / 插值 / BM / 高斯消元 / 行列式 / 线性基 …） | 23 项全部通过 |
| `tools/verify_misc.cpp` | 组合、几何、博弈、单纯形、Burnside 对拍 | 21 项全部通过 |

运行方法见 [`tools/README.md`](tools/README.md)。

> **对拍真的有用**：验证过程中发现并修正了教程里单纯形法实现的一个真实缺陷（漏存松弛变量列导致提前终止）。详见 `tools/README.md` 第 3 节。

---

## 代码约定

### 常用宏与头文件

```cpp
#include <bits/stdc++.h>
using namespace std;
using ll  = long long;      // 64 位整数
using ull = unsigned long long;
using i128 = __int128;      // 处理 1e18 * 1e18 这类乘法溢出
```

**约定**：本教程中出现的所有模板都假设以上 `using` 已经存在。

### 关于溢出

竞赛数学题 90% 的 RE / WA 都是溢出。记住三条铁律：

1. `int * int` 可能溢出 `int`，遇到 $a \le 10^9$ 的乘法一律先转 `long long`；
2. `(a * b) % mod` 中若 $a, b < 2^{63}$ 量级（比如模数 $10^{18}$），必须用 `__int128` 或龟速乘；
3. `long long` 的乘法溢出是 **UB（未定义行为）**，不要指望它会绕回成负数然后再取模。

```cpp
// 稳妥的模乘：适用于 mod 达到 1e18 量级
ll mul_mod(ll a, ll b, ll mod) {
    return (ll)((__int128)a * b % mod);
}
```

### 快速取模的常数优化

```cpp
// 法一：少用 %，用条件减法（要求 a+b 不溢出）
inline int add(int a, int b, int mod) { a += b; return a >= mod ? a - mod : a; }

// 法二：无符号自动取模的 Barrett / 技巧，适合 mod < 2^31 且追求极限常数
inline int mul(int a, int b, int mod) { return (ll)a * b % mod; }
```

### 模数封装模板（后续章节反复使用）

```cpp
template <int MOD>
struct ModInt {
    int v;
    ModInt(ll x = 0) { v = (int)(x % MOD); if (v < 0) v += MOD; }
    int val() const { return v; }
    ModInt operator-() const { return ModInt(v ? MOD - v : 0); }
    ModInt &operator+=(const ModInt &o) { v += o.v; if (v >= MOD) v -= MOD; return *this; }
    ModInt &operator-=(const ModInt &o) { v -= o.v; if (v < 0) v += MOD; return *this; }
    ModInt &operator*=(const ModInt &o) { v = (int)((ll)v * o.v % MOD); return *this; }
    ModInt pow(ll e) const {
        ModInt r(1), a(*this);
        for (; e; e >>= 1, a *= a) if (e & 1) r *= a;
        return r;
    }
    ModInt inv() const { return pow(MOD - 2); }   // 仅当 MOD 为质数
    friend ModInt operator+(ModInt a, const ModInt &b) { return a += b; }
    friend ModInt operator-(ModInt a, const ModInt &b) { return a -= b; }
    friend ModInt operator*(ModInt a, const ModInt &b) { return a *= b; }
    friend bool operator==(const ModInt &a, const ModInt &b) { return a.v == b.v; }
};
using mint = ModInt<998244353>;
```

> **注意**：`ModInt` 的 `inv()` 用费马小定理，只对质数模数成立。
> 模数非质数（比如 $10^9$、$10^9+6$）时必须用扩展欧几里得求逆元，见第 2 章。

### 复杂度标记

全文统一：$n$ 为数据规模，$V$ 为值域，$m$ 为模数，$p$ 为质数，$\omega(n)$ 为 $n$ 的不同质因子个数（通常 $\le 9$）。

---

## 如何使用这份教程

1. **不要只读代码**。每个定理自己推一遍，尤其是莫比乌斯反演和生成函数，光看不算等于没学。
2. **代码要自己敲一遍**。抄一遍理解不到边界条件；本教程代码都标了关键边界。
3. **配题服用**。每章末尾的"例题套路"给出的是题目特征 → 方法的映射，试着在 OJ 上按标签搜题练习。
4. **建立自己的板子库**。把每章的模板整理成自己的 `math.cpp`，比赛时直接复制。
