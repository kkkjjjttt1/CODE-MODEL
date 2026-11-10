# 第 02 章 字符串哈希与 Rabin–Karp

> **一句话**：把字符串映射成"大整数"，从而 **O(1) 判断任意两个子串是否相等**。
> 这是性价比最高的字符串工具——简单、强大、易错点少，但要注意**碰撞**与**被卡**的问题。

---

## 2.1 核心思想：多项式滚动哈希

把串看成 **base 进制的大整数**（字符值即每一位的数字），再对一个**模数 M** 取模：

```
s[0..n-1]  →  H(s) = ( s[0]*base^{n-1} + s[1]*base^{n-2} + ... + s[n-1] ) mod M
```

这样相等的串哈希必相等；不相等的串哈希有极小概率相等（**碰撞**），所以它是"概率算法"。

### 前缀哈希（1 基下标，最常用）

设 `h[i]` 表示 `s` 的前缀 `s[0..i-1]`（长度 `i`）的哈希，递推：

```
h[0] = 0
h[i] = h[i-1] * base + val(s[i-1])        (mod M)
```

其中 `val(c)` 可取 `c` 的 ASCII 值（不取 0，避免 "a" 与 "aa" 前导零问题；`'a'` 用 ASCII 97 本身也安全）。

### 区间哈希（关键公式）

要取 `s[l..r]`（**1 基、闭区间**）的哈希，只需：

```
hash(s[l..r]) = h[r] - h[l-1] * base^{r-l+1}     (mod M)
```

直观理解：`h[r]` 是把 `s[0..r-1]` 当 base 进制数，`h[l-1]*base^{len}` 就是把高位的 `l-1` 个字符"左移"出区间长度后再减掉，剩下的就是区间内的数字。

> **为什么用前缀而不是直接枚举子串**：预处理 `h[1..n]` 只需 O(n)，之后每个区间查询 O(1)。若把每个子串单独算哈希，总复杂度会是 O(n²)。

### 两种取模实现

| 方式 | 写法 | 优点 | 缺点 |
|------|------|------|------|
| **自然溢出** | `unsigned long long`，乘法自然溢出即模 `2^64` | 快、不用写 mod | 可被**构造数据卡**（见 2.4） |
| **大质数取模** | 模 `1e9+7` / `1e9+9` | 安全得多 | 稍慢；乘法注意 `long long` |

---

## 2.2 代码模板

### 模板一：自然溢出（快，模板题首选）

```cpp
using ull = unsigned long long;
struct StrHash {
    vector<ull> h, pw;
    // base 建议选一个比字符集大的"随机数"，如 131 / 13331 / 随机
    StrHash(const string& s, ull base = 13331) { build(s, base); }

    void build(const string& s, ull base) {
        int n = (int)s.size();
        h.assign(n + 1, 0);
        pw.assign(n + 1, 1);
        for (int i = 1; i <= n; ++i) {
            h[i] = h[i - 1] * base + (ull)s[i - 1];  // 自然溢出 = mod 2^64
            pw[i] = pw[i - 1] * base;
        }
    }
    // s[l..r] 的哈希，l、r 均为 1 基闭区间
    ull get(int l, int r) const {
        return h[r] - h[l - 1] * pw[r - l + 1];
    }
};
```

### 模板二：双哈希（抗卡首选）

用两个不同的大质数分别算，碰撞概率降到约 `1/(M1*M2)`，几乎不可能被卡。

```cpp
struct DoubleHash {
    using ll = long long;
    static constexpr ll M1 = 1000000007LL;   // 1e9+7
    static constexpr ll M2 = 1000000009LL;   // 1e9+9
    ll base;
    vector<ll> h1, h2, pw1, pw2;

    DoubleHash(const string& s, ll base_ = 911382323LL) : base(base_) { build(s); }

    void build(const string& s) {
        int n = (int)s.size();
        h1.assign(n + 1, 0); h2.assign(n + 1, 0);
        pw1.assign(n + 1, 1); pw2.assign(n + 1, 1);
        for (int i = 1; i <= n; ++i) {
            h1[i] = (h1[i - 1] * base + s[i - 1]) % M1;
            h2[i] = (h2[i - 1] * base + s[i - 1]) % M2;
            pw1[i] = pw1[i - 1] * base % M1;
            pw2[i] = pw2[i - 1] * base % M2;
        }
    }
    pair<ll, ll> get(int l, int r) const {   // 1 基闭区间
        ll a = (h1[r] - h1[l - 1] * pw1[r - l + 1] % M1 + M1) % M1;
        ll b = (h2[r] - h2[l - 1] * pw2[r - l + 1] % M2 + M2) % M2;
        return {a, b};
    }
};
```

> 说明：两个结构都假设字符值不为 0 的倍数问题不大；若用 `mod` 版本且 `base` 与 `M` 接近（`911382323 < 1e9+7`），乘法 `h*base` 均 < 1e18，在 `long long` 内安全。

---

## 2.3 哈希的五大基本应用

### 2.3.1 O(1) 判子串相等 → 字符串匹配（Rabin–Karp）

把模式串 `p`（长 `m`）的哈希算好，在文本 `s` 上滑动长度 `m` 的窗口，逐一比较哈希：

```cpp
// 返回 p 在 s 中出现的所有起始下标（0 基）
vector<int> rabinKarp(const string& s, const string& p) {
    int n = (int)s.size(), m = (int)p.size();
    vector<int> res;
    if (m > n) return res;
    ull base = 911382323;                       // 或换 DoubleHash
    StrHash hs(s, base), hp(p, base);
    for (int l = 1; l + m - 1 <= n; ++l)        // 窗口 [l, l+m-1]，1 基
        if (hs.get(l, l + m - 1) == hp.get(1, m))
            res.push_back(l - 1);               // 转回 0 基
    return res;
}
```

复杂度 O(n + m)。它比 KMP 好写，且天然支持**多模式**（把每个模式串的哈希放进哈希表）。最坏情况理论上是 O(nm)（每次窗口哈希都碰撞才需验证），实际中不验证也极少出错。

### 2.3.2 判断回文：正反两个哈希

对原串 `s` 与反串 `rev(s)` 各建一套哈希，则

```
s[l..r] 是回文  ⇔  hash(s[l..r]) == hash(反串中对应段)
```

```cpp
// 判断 s[l..r] 是否回文（1 基闭区间）
bool isPal(const StrHash& hf, const StrHash& hb, int l, int r) {
    // hf: 正串哈希；hb: 反串哈希（hb 的坐标按反串计算）
    // 原串下标 [l-1, r-1]（0 基）在反串中的下标为 [n-r, n-l]（0 基）
    int L = (int)hf.h.size() - 1 - r + 1;   // = n - r + 1（1 基）
    int R = (int)hf.h.size() - 1 - l + 1;   // = n - l + 1（1 基）
    return hf.get(l, r) == hb.get(L, R);
}
// 建反串哈希：
// string rs = s; reverse(rs.begin(), rs.end()); StrHash hb(rs);
```

### 2.3.3 最长公共前缀 LCP（二分答案）

想求同一串内两个位置 `i, j`（0 基）出发的**最长公共前缀长度**：若 `s[i..i+k-1] == s[j..j+k-1]` 则 `k` 可行——而相等恰好可以用哈希 O(1) 判断，`k` 有单调性，因此二分：

```cpp
// s 内从 0 基位置 a、b 开始的 LCP（配合 StrHash 的 h/pw）
int lcp(const StrHash& H, int a, int b) {
    int n = (int)H.h.size() - 1;             // 串长
    int lo = 0, hi = min(n - a, n - b);      // [0, hi]
    while (lo < hi) {
        int mid = (lo + hi + 1) >> 1;
        if (H.get(a + 1, a + mid) == H.get(b + 1, b + mid)) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}
```

LCP 用处极大：**比较两个子串的字典序**（先比 LCP，再看下一个字符）→ 可以 O(n log n) 对"若干子串"排序；也是后缀数组（09 章）里 RMQ 做法的哈希平替。

### 2.3.4 找最小周期 / 判断循环节

`p` 是 `s`（长 `n`）的周期 ⇔ 首尾两段长 `n-p` 的子串相等：

```cpp
// p 是否是 s 的周期（p < n）
bool isPeriod(const StrHash& H, int n, int p) {
    return H.get(1, n - p) == H.get(p + 1, n);
}
```

想找**最小周期**：枚举能整除/不整除的候选或对长度因子做检查；`p` 整除 `n` 时 `s` 完全由 `s[0..p-1]` 重复得到（"整周期"，最小整循环节）。例：`"abcabcabc"` 周期 3；`"abababa"` 的**最小周期是 2**（虽然 7 不被 2 整除，仍有周期 2——别被整除条件误导）。线性求法见 03 章 KMP 的 `n - border` 结论。

### 2.3.5 两个串的最长公共子串（哈希 + 二分）

二分长度 `L`：把串 A 中所有长度 `L` 的子串哈希插入 `unordered_set`，再检查串 B 中是否存在长度 `L` 的子串哈希命中。

```cpp
int longestCommonSubstring(const string& A, const string& B) {
    // 用 DoubleHash 或 (自然溢出 + unordered_set<ull>)
    auto ok = [&](int L) -> bool {
        if (L == 0) return true;
        if (L > (int)A.size() || L > (int)B.size()) return false;
        StrHash ha(A, 911382323), hb(B, 911382323);
        unordered_set<ull> st;
        for (int l = 1; l + L - 1 <= (int)A.size(); ++l)
            st.insert(ha.get(l, l + L - 1));
        for (int l = 1; l + L - 1 <= (int)B.size(); ++l)
            if (st.count(hb.get(l, l + L - 1))) return true;
        return false;
    };
    int lo = 0, hi = min((int)A.size(), (int)B.size());
    while (lo < hi) { int mid = (lo + hi + 1) >> 1; if (ok(mid)) lo = mid; else hi = mid - 1; }
    return lo;
}
```

复杂度 O((n+m) log n)。注意每次 check 重建哈希略浪费，可提前建好后只枚举窗口（复杂度不变）。O(n) 线性做法见 09 章后缀数组 / 10 章后缀自动机。

---

## 2.4 碰撞与"被卡"：必须知道的坑

1. **自然溢出可以被构造数据击穿**：CF 上有著名讨论（anti-hash test），用 Thue–Morse 序列等可构造大量碰撞，UOJ 甚至出过专门卡 `base=131` + 自然溢出的数据。**省选/ICPC 级别建议双哈希**；单哈希选大质数模也行。
2. **固定 base + 固定单模（如只有 1e9+7）** 也可能被针对性构造碰撞，只是成本高；稳妥是 `双模` 或 `随机 base`（运行时 `mt19937` 生成 base）。
3. 用 `int` 存乘积会溢出：模版乘法务必用 `long long`（`(h*base)%M` 中 `h*base` 约为 1e18，必须 `long long`）。
4. `unordered_set` / `unordered_map` 存哈希做统计时，**碰撞会直接导致答案错**，本地小数据测不出来——所以统计类题目（如数本质不同子串）能用双哈希就用双哈希，或用 SA/SAM。
5. 模减法出现负数要 `+M` 再 `%M`（模板里已处理）。
6. 字符 `val(c)` 建议用 `(ull)c` 或 `c - 'a' + 1`：若映射到 0，串 `"a"` 与 `"aa"` 会得到相同前缀值导致某些公式出错，映射到 1..26 更安全。

---

## 2.5 经典例题

| 题目 | 考点 | 要点 |
|------|------|------|
| [洛谷 P3370 【模板】字符串哈希](https://www.luogu.com.cn/problem/P3370) | 去重计数 | 每串求 (h1,h2) 插入 `set<pair<ll,ll>>`，输出 size |
| [洛谷 P2957 [USACO09OCT] Barn Echoes G](https://www.luogu.com.cn/problem/P2957) | 哈希判相等 | 枚举公共前后缀长度：A 前缀 == B 后缀，或反之；哈希 O(1) 判 |
| [洛谷 P3805 【模板】Manacher](https://www.luogu.com.cn/problem/P3805)（对照用） | 哈希二分判回文 | 可先用"正反哈希+二分"过小数据理解题意，再学 07 章线性做法 |
| UVA 11475 *Extend to Palindrome* | 最长回文后缀 | 哈希 O(1) 判回文：找最长的 `l` 使 `s[l..n-1]` 是回文，把前半补到末尾 |
| [洛谷 P3375 【模板】KMP](https://www.luogu.com.cn/problem/P3375)（对照用） | 匹配 | Rabin–Karp 也能 AC 此题（求所有出现位置） |
| CF 1320D *Reachable Strings*（进阶） | 奇偶位哈希 | 把 0 的奇偶位置分别哈希；难度较高，入门后可挑战 |

**双哈希去重模板写法**（P3370 常用）：

```cpp
int n; cin >> n;
set<pair<long long, long long>> st;
for (int i = 0; i < n; ++i) {
    string s; cin >> s;
    DoubleHash H(s);
    st.insert(H.get(1, (int)s.size()));
}
cout << st.size() << '\n';
```

---

## 2.6 延伸预告

- 哈希 + **线段树/平衡树**：支持单点修改的动态字符串比较（见 14 章）；
- 哈希 + **树上路径**：树上任意路径串的哈希与 LCP（见 14 章）；
- **二维哈希**：矩阵子矩阵匹配；
- 哈希判回文 O(n log n)，但线性方法见 07 章 Manacher。

> 本章最该记住的两句话：
> 1. **区间哈希公式 `h[r] - h[l-1] * base^{len}`（1 基）**；
> 2. **大比赛别用裸自然溢出，双哈希或随机 base**。
