# 第 13 章 bitset、Shift-And 与 FFT/NTT 字符串匹配

> **一句话**：当普通算法的复杂度卡在"×n"时，两种武器可再降一档：
> - **位并行**（bitset / Shift-And）：把一次字符比较变成一次 **64 位字的位运算**，把 n 变成 n/64；
> - **卷积匹配**（FFT/NTT）：把"枚举所有对齐位置"的 n² 次比较变成一次 O(n log n) 的多项式乘法，
>   尤其擅长**带通配符**或**容错**匹配。
> 前置：会用 std::bitset；懂"多项式乘法=卷积"的概念（本册不展开多项式全家桶，给最小 NTT 模板）。

---

## 13.1 位并行思想

普通循环一次只处理一个字符比较；现代 CPU 一条指令可对 64 个 bit 同时运算。很多字符串 DP/匹配算法如果能把"比较结果"编码成 bit，就能一次推进 64 个字符的进度——这就是**位并行（bit-parallel）**。

典型应用：
- Shift-And（单模式匹配，支持通配符）；
- 若干"每次询问 ≤ 1e4 次、文本短"的匹配题用 bitset 暴力反而简单；
- LCS / 编辑距离的 bitset 优化（O(nm/64)）。

---

## 13.2 Shift-And 算法

**思想**：维护位掩码 `D`，第 k 位 = 1 表示"当前文本的某个后缀等于 p[0..k]"（即已匹配长度为 k+1）。读入文本字符 c：
- `D = ((D << 1) | 1) & B[c]`，其中 `B[c]` 的第 k 位 = 1 ⇔ `p[k] == c`（模式里能匹配 c 的所有位置）；
- 若 `D` 的最高位（第 m−1 位）为 1 → 模式 p 完整匹配，结束于当前位置。

### 13.2.1 ull 版（m ≤ 64）

```cpp
// Shift-And：返回 p 在 s 中所有出现位置（0 基起点）；要求 m ≤ 64
vector<int> shiftAndULL(const string& s, const string& p) {
    int n = (int)s.size(), m = (int)p.size();
    if (m == 0 || m > 64) return {};
    unsigned long long B[256] = {};
    for (int i = 0; i < m; ++i)
        B[(unsigned char)p[i]] |= (1ULL << i);
    vector<int> pos;
    unsigned long long D = 0, hit = 1ULL << (m - 1);
    for (int i = 0; i < n; ++i) {
        D = ((D << 1) | 1) & B[(unsigned char)s[i]];
        if (D & hit) pos.push_back(i - m + 1);   // 结束位置 i，起点 i-m+1
    }
    return pos;
}
```

### 13.2.2 通配符与 bitset 版（m 任意）

若模式里有通配符 `?`（匹配任意一个字符），只需让 `?` 所在位对所有字符都可匹配：

```cpp
// bitset 版：m 可到 1e5；复杂度 O(n · m/word)
// 字符集假定小写字母 + '?'
vector<int> shiftAndWildcard(const string& s, const string& p) {
    int n = (int)s.size(), m = (int)p.size();
    static bitset<100005> B[26], D;
    D.reset();
    for (int i = 0; i < m; ++i)
        if (p[i] == '?') {
            for (int c = 0; c < 26; ++c) B[c].set(i);
        } else B[p[i] - 'a'].set(i);

    vector<int> pos;
    for (int i = 0; i < n; ++i) {
        D = ((D << 1) | 1) & B[s[i] - 'a'];
        if (D[m - 1]) pos.push_back(i - m + 1);
    }
    return pos;
}
```

**说明**：bitset 的 `<<1` 整体左移一次是 O(m/word)，所以总复杂度 O(n·m/word)，对"短文本 × 长模式"或"多询问但总移动量可控"的题很实用（如 CF 914F 一类的动态字符集询问）。注意：D 里旧的高位会被 `& B[c]` 清掉（因为模式只有 m 位），不会污染。

> 与 KMP 对比：Shift-And 常数小、支持通配符、写起来短；缺点是复杂度和 m 相关，且不适合"模式串很长且文本也很长"（那是 O(nm/64)，仍可能超时，用 KMP 更稳）。

---

## 13.3 bitset 加速其它 DP（简介）

- **LCS 的 bitset 算法**：利用"行间差分只有 0/±1"的性质把整行转移写成位运算，复杂度 O(nm/64)。实现见 CP-Algorithms 或 Hyyrö 的论文；需要时再抄板子，日常 DP 题 n,m ≤ 1e4 用滚动数组即可。
- **编辑距离的 band / bitset 变体**：允许 k 个错误的近似匹配，可结合位并行（BLAST 思想）。

> 记忆点：*bitset 适合"把一次内层循环压成位运算"的场景；写题时先看能否把条件变成"按位与/或/移位"。*

---

## 13.4 FFT/NTT 匹配（带通配符 / 容错匹配）

### 13.4.1 通配符匹配的卷积化

**问题**（[洛谷 P4173 残缺的字符串](https://www.luogu.com.cn/problem/P4173)）：模式 P（长 m）与文本 T（长 n）都**可以含通配符 `*`**，求 P 在 T 的所有匹配位置。

**关键观察**：把字符映射成整数，通配符映射成 **0**。若两个字符相等则 `(a−b)² = 0`。于是对某个对齐位置（P 对齐到 T[j..j+m-1]），定义

```
S(j) = Σ_{i=0..m-1}  P[i] · T[j+i] · (P[i] − T[j+i])²
```

- 通配符（值 0）那一项贡献 0，自动"什么都匹配"；
- 其余位置字符全相等时每项都是 0 → **S(j) = 0 ⇔ 匹配**（双方任意一方是通配符也算匹配）。

展开：`S(j) = Σ P³T − 2 Σ P²T² + Σ PT³`。每一项都是**卷积**：把 P 反转，则
`Σ_i P[i]·X[j+i]` 这类交叉和 = 卷积数组的第 `j+m−1` 项。于是三次多项式乘法即可一次算出所有位置的 S。

```cpp
// 约定 multiply(a, b) 返回 a 与 b 的卷积（下面给出最小 NTT 模板）
// P 长 m，T 长 n；值 1..26，通配符 '*' → 0
vector<int> wildcardMatch(const string& P, const string& T) {
    int m = P.size(), n = T.size();
    vector<long long> pa(m), pb2(m), pc3(m), ta(n), tb2(n), tc3(n);
    for (int i = 0; i < m; ++i) {
        long long v = (P[i] == '*') ? 0 : (P[i] - 'a' + 1);
        pa[m - 1 - i] = v;                  // 反转模式
        pb2[m - 1 - i] = v * v % MOD;
        pc3[m - 1 - i] = v * v % MOD * v % MOD;
    }
    for (int j = 0; j < n; ++j) {
        long long v = (T[j] == '*') ? 0 : (T[j] - 'a' + 1);
        ta[j] = v;
        tb2[j] = v * v % MOD;
        tc3[j] = v * v % MOD * v % MOD;
    }
    auto c1 = multiply(pc3, ta);            // Σ P³·T
    auto c2 = multiply(pb2, tb2);           // Σ P²·T²
    auto c3 = multiply(pa, tc3);            // Σ P·T³

    vector<int> res;
    for (int j = 0; j + m <= n; ++j) {      // 对齐起点 j
        int k = j + m - 1;                  // 卷积下标
        long long s = (c1[k] - 2 * c2[k] + c3[k]) % MOD;
        if (s == 0) res.push_back(j);       // mod 意义下为 0（详见 13.4.3 防碰撞）
    }
    return res;
}
```

**复杂度**：3 次卷积 O(n log n)，对 n,m ≤ 1e5 轻松通过；这是 P4173 的正解。

### 13.4.2 最小 NTT 模板（mod = 998244353，原根 3）

```cpp
using ll = long long;
const ll MOD = 998244353, G = 3;

ll modPow(ll a, ll e) {
    ll r = 1;
    while (e) { if (e & 1) r = r * a % MOD; a = a * a % MOD; e >>= 1; }
    return r;
}
void ntt(vector<ll>& a, bool inv) {
    int n = (int)a.size();
    for (int i = 1, j = 0; i < n; ++i) {          // 位反转
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        ll wlen = modPow(G, (MOD - 1) / len);
        if (inv) wlen = modPow(wlen, MOD - 2);
        for (int i = 0; i < n; i += len) {
            ll w = 1;
            for (int j = 0; j < len / 2; ++j) {
                ll u = a[i + j], v = a[i + j + len / 2] * w % MOD;
                a[i + j] = (u + v) % MOD;
                a[i + j + len / 2] = (u - v + MOD) % MOD;
                w = w * wlen % MOD;
            }
        }
    }
    if (inv) {
        ll invN = modPow(n, MOD - 2);
        for (auto& x : a) x = x * invN % MOD;
    }
}
vector<ll> multiply(vector<ll> a, vector<ll> b) {
    int need = (int)a.size() + (int)b.size() - 1, n = 1;
    while (n < need) n <<= 1;
    a.resize(n); b.resize(n);
    ntt(a, false); ntt(b, false);
    for (int i = 0; i < n; ++i) a[i] = a[i] * b[i] % MOD;
    ntt(a, true);
    a.resize(need);
    return a;
}
```

（值域大或需要精确整数时改用 FFT+long double 或 Python 大整数；数值都非负时也可把三次卷积换成"每字符一次卷积"的字符分解法：对小字符集，逐字符做指示数组卷积并累加，能避免大值取模碰撞。）

### 13.4.3 模碰撞的小心机

`S(j)` 的真实值是**非负整数**，可能超过 MOD（如 m=1e5、字符值 26 时，ΣPT(P−T)² 最大 ≈ 4×10¹⁰）。只在模意义下判 0 有极小概率误判。稳妥做法：**用两组不同的模数/原根各算一遍**，或先按 13.4.2 注释里的"逐字符分解"做法（每个字符的匹配互相独立、和值 ≤ m×1 ≤ 1e5 < MOD，无碰撞）。逐字符做法如下：

```
对每个字符 c ∈ Σ：
  令 Pa[i] = (P[i]=='*' || P[i]==c) ? 1 : 0（反转后）
     Ta[j] = (T[j]=='*' || T[j]==c) ? 1 : 0
  卷积 Ca = conv(Pa, Ta)
  S[j] += Ca[j+m-1]     // 该字符在所有对齐位置的对齐个数
最后 S[j] == m 的位置即完全匹配（m 次 ≤ 1e5 不会溢出 mod）
```

这样每个卷积的系数都是"对齐计数 ≤ m"，安全无碰撞，代价是 |Σ| 次卷积（26 次 NTT 仍为 O(|Σ| n log n)，1e5 内可接受）。两种写法任选其一。

### 13.4.4 其它卷积匹配变体

- **容错（模糊）匹配**：模式与文本每个位置给分，Σ(P−T)² 最小化 → 允许 k 个错配的近似匹配、生物序列比对；
- **二进制串近似匹配**：把 0/1 串直接卷 Σ(P xor T)；
- 把 13.4.1 的 `(a−b)²` 换成其它代价函数可得到不同语义的匹配。

---

## 13.5 经典例题

| 题目 | 考点 | 要点 |
|------|------|------|
| [洛谷 P4173 残缺的字符串](https://www.luogu.com.cn/problem/P4173) | FFT 通配符匹配 | 三次卷积或逐字符分解 |
| CF 528D *Fuzzy Search* | FFT 容错匹配 | 每字符卷积，窗口错配数 ≤ k |
| CF 914F *Substrings in a String* | bitset 匹配 | 动态字符、短串询问：bitset 位移法 |
| SPOJ NHAY / 洛谷 KMP 模板（复习） | 普通匹配 | Shift-And 可作对拍 |
| UVA 1056?（近似匹配类，选做） | 位并行 | 视兴趣 |

---

## 13.6 小结

```
匹配问题的武器谱：
  精确匹配（单模式）        KMP / Z / Shift-And
  精确匹配（多模式）        AC 自动机
  前缀/后缀查询             SA / SAM
  带通配符匹配              FFT 卷积（O(n log n)）或 bitset（O(nm/64)）
  容错匹配                  FFT 距离卷积
  小范围多次询问            bitset 位移 / 哈希
```

> 下一章（14）把本册各种工具与**数据结构**组装起来：后缀树、后缀平衡树、"线段树/平衡树维护哈希"、树上字符串问题——字符串进阶题的常客。
