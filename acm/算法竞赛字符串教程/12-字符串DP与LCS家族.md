# 第 12 章 字符串动态规划与 LCS 家族

> **一句话**：字符串 DP 分两大类——
> ① **经典序列 DP**：LCS、LCIS、编辑距离、回文子序列，套路固定、模板好抄；
> ② **状态机 DP**：把"当前匹配到哪"（KMP/AC/SAM/子序列自动机的状态）写进 DP 下标，处理"构造/修改串使其满足或不满足某些子串约束"的问题——这是字符串大题的真正门槛。
> 本章给出两类问题的模板、公式与例题导航。

---

## 12.1 LCS：最长公共子序列

### 12.1.1 经典 O(nm) DP

```
dp[i][j] = LCS(A[1..i], B[1..j])
A[i]==B[j] : dp[i][j] = dp[i-1][j-1] + 1
否则        : dp[i][j] = max(dp[i-1][j], dp[i][j-1])
```

```cpp
int LCS(const string& A, const string& B) {
    int n = A.size(), m = B.size();
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= m; ++j) {
            if (A[i - 1] == B[j - 1]) dp[i][j] = dp[i - 1][j - 1] + 1;
            else dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
        }
    return dp[n][m];
}
// 空间优化：只保留上一行（滚动数组）；回溯方案：从 (n,m) 按转移反向走并记录字符
```

### 12.1.2 特例：两序列都是排列 → 转 LIS O(n log n)

若 A、B 各自是 0..n-1 的排列：把 A 中每个值的"位置"记下来，求 B 中这些位置构成的 **LIS**（最长上升子序列），长度即 LCS。LIS 用经典 O(n log n) 贪心+二分。

> [洛谷 P1439 最长公共子序列](https://www.luogu.com.cn/problem/P1439)（排列版）就是这道模板题。

### 12.1.3 bitset 优化 LCS（O(nm/w)）

把行 DP 写成位运算形式（"seaweed 算法"），n,m ≤ 5e4 级别可行；实现见第 13 章的 bitset 小节。

---

## 12.2 最长公共子串（DP vs 后缀结构）

| 方法 | 复杂度 | 说明 |
|---|---|---|
| 二维 DP `dp[i][j]=A[i]==B[j]? dp[i-1][j-1]+1 : 0` | O(nm) | 只适合 n,m ≤ 5e3 |
| 哈希 + 二分 | O((n+m) log n) | 见 02 章 |
| **后缀数组/后缀自动机** | **O(n+m)** | 见 09 / 10 章，大数据的正解 |

（提示：问"最长公共子串"时先看数据范围，>1e4 基本就是 SA/SAM 题。）

---

## 12.3 LCIS：最长公共上升子序列（O(nm)）

```cpp
// A[1..n], B[1..m]；dp[j] = 以 B[j] 结尾的 LCIS 长度（滚动一维）
int LCIS(const vector<int>& A, const vector<int>& B) {
    int n = A.size(), m = B.size();          // 若从 1 基读入请相应调整
    vector<int> dp(m + 1, 0);
    for (int i = 0; i < n; ++i) {
        int best = 0;
        for (int j = 0; j < m; ++j) {
            if (B[j] < A[i]) best = max(best, dp[j]);      // 可接在后面且递增
            else if (B[j] == A[i]) dp[j] = max(dp[j], best + 1);
        }
    }
    return *max_element(dp.begin(), dp.end());
}
```

---

## 12.4 编辑距离（Levenshtein Distance）

允许三种操作：插入、删除、替换（代价 1）：

```
dp[i][j] = 把 A[1..i] 变成 B[1..j] 的最小代价
A[i]==B[j] : dp[i][j] = dp[i-1][j-1]
否则        : dp[i][j] = 1 + min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1])
                                删除A[i]   插入B[j]     替换A[i]
```

```cpp
int editDistance(const string& A, const string& B) {
    int n = A.size(), m = B.size();
    vector<vector<int>> dp(n + 1, vector<int>(m + 1));
    for (int i = 0; i <= n; ++i) dp[i][0] = i;
    for (int j = 0; j <= m; ++j) dp[0][j] = j;
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= m; ++j)
            dp[i][j] = (A[i - 1] == B[j - 1] ? dp[i - 1][j - 1]
                        : 1 + min({dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]}));
    return dp[n][m];
}
```

**变体**：只允许插入/删除（不允许替换）→ `min(dp[i-1][j], dp[i][j-1]) + 1` 且不相等时不可用对角线；带权编辑距离把 1 换成权；**LCS 与编辑距离关系**：只允许删除的"编辑距离"= n + m − 2·LCS。

**与自动机结合**：把 A 换成"AC 自动机状态"（DNA repair 类题），见 12.6。

---

## 12.5 回文类 DP

### 12.5.1 最长回文子序列 LPS（O(n²)）

```cpp
// dp[l][r]：s[l..r] 的最长回文子序列长度
for (int len = 1; len <= n; ++len)
    for (int l = 0, r = l + len - 1; r < n; ++l, ++r) {
        if (s[l] == s[r]) dp[l][r] = (len == 1 ? 1 : dp[l + 1][r - 1] + 2);
        else dp[l][r] = max(dp[l + 1][r], dp[l][r - 1]);
    }
// 注：LPS(s) 也等于 LCS(s, reverse(s))
```

### 12.5.2 最少回文分割 / 回文划分

- 朴素 O(n²)：`cut[i] = min{ cut[j-1] + 1 | s[j..i] 是回文 }`（[UVA 11584 Partitioning](https://onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&page=show_problem&problem=2630)）；
- Manacher / PAM 加速：回文后缀按长度分段（见 07、08 章），O(n log n)（CF 932G 一类的 slink 优化）。

---

## 12.6 状态机 DP（本册进阶核心）

通用套路：**构造/修改一个目标串 T，使其满足（或不满足）某些"包含/不包含模式"的约束，求方案数/最值**。
记 `dp[i][状态]`：处理完 T 的前 i 个字符后，T 与模式集合的"匹配状态"为 `状态`（KMP 的已匹配前缀长度 / AC 自动机节点 / SAM 节点），转移就是自动机读入一个字符。**自动机的确定性转移让 DP 天然不重不漏**。

### 12.6.1 KMP 自动机 + DP（单个禁止模式）

先补全 KMP 转移表 `go[state][c]`（state = 已匹配模式前缀长度 0..m）：

```cpp
// 对模式 p 建 go[state][c]：state 表示已匹配 p[0..state-1]
// go 满足：读入字符 c 后新的已匹配长度
vector<array<int, 26>> buildGo(const string& p) {
    int m = p.size();
    auto pi = prefixFunction(p);                 // 03 章
    vector<array<int, 26>> go(m + 1);
    for (int state = 0; state <= m; ++state)
        for (int c = 0; c < 26; ++c) {
            if (state < m && c == p[state] - 'a') go[state][c] = state + 1;
            else if (state == 0) go[state][c] = 0;
            else go[state][c] = go[pi[state - 1]][c];   // 已算过（编号更小）
        }
    return go;
}
```

**应用**：数"长度 n、不含模式 p 的串"个数：

```cpp
// dp[i][state]：长度 i、匹配状态 state（<m）的串数；禁止进入状态 m
vector<vector<int>> dp(n + 1, vector<int>(m, 0));
dp[0][0] = 1;
for (int i = 0; i < n; ++i)
    for (int s = 0; s < m; ++s)
        for (int c = 0; c < 26; ++c) {
            int ns = go[s][c];
            if (ns < m) dp[i + 1][ns] = (dp[i + 1][ns] + dp[i][s]) % MOD;
        }
ans = Σ_s dp[n][s];
// n 很大（1e9）时把 go 的 0/1 矩阵做矩阵快速幂，复杂度 O(m^3 log n)
```

### 12.6.2 AC 自动机 + DP（多禁止模式）

见 06 章 6.5：先标记 "bad"（含任一模式的状态），再在上面做同样 DP。
经典题：[P4052 [JSOI2007] 文本生成器](https://www.luogu.com.cn/problem/P4052)（含至少一个模式：总数 − 不含）、POJ 3691 / HDU 2457 DNA repair（**最少修改**，用 min 转移）。

### 12.6.3 SAM / SA 上的 DP

- "第 k 小子串"：SAM 路径计数 + 贪心（10 章⑤）；
- "统计字典序第 k 大路径 / 至少出现 k 次子串的最长长度"等，本质都是 **DAG 上 DP**（SAM 的转移图是 DAG）；
- SAM link 树 + 树上 DP / 线段树合并处理区间约束（进阶）。

### 12.6.4 子序列自动机 + DP

**本质不同子序列计数**（含空串，O(n)）：

```cpp
// s 长 n，字符集小写；f[i] = s[0..i-1] 的本质不同子序列数（含空串）
long long distinctSubsequences(const string& s, long long MOD) {
    vector<long long> f(s.size() + 1, 1);      // f[0] = 1（空）
    vector<int> last(26, -1);                  // 该字符上一次出现的位置
    for (int i = 0; i < (int)s.size(); ++i) {
        int c = s[i] - 'a';
        f[i + 1] = (2 * f[i]) % MOD;
        if (last[c] != -1)                     // 去掉重复：上次以该字符结尾的集合
            f[i + 1] = (f[i + 1] - f[last[c]] + MOD) % MOD;
        last[c] = i;                           // 记录本次位置
    }
    return f[s.size()];
}
// 验证：s="aaa" → f = [1,2,3,4]（空,a,aa,aaa）✓；s="aba" → 7 ✓
```

> 子序列自动机（05 章 5.6 的 nxt 表）与这类 DP 结合，还能做"长为 k 的本质不同子序列计数"等变体。

### 12.6.5 字符串上的"区间覆盖/最优化"DP

- 回文覆盖、区间分段、词拆分（word break）：`dp[i] = f(dp[j] | s[j..i] 满足条件)`，配合哈希/Trie/自动机 O(1) 判定条件；
- [P2292 [HNOI2004] L语言](https://www.luogu.com.cn/problem/P2292)：词典 Trie + 对文本前缀做可行性 DP/记忆化 BFS。

---

## 12.7 经典例题

| 题目 | 类型 | 要点 |
|------|------|------|
| [洛谷 P1439 最长公共子序列](https://www.luogu.com.cn/problem/P1439) | 排列 LCS→LIS | O(n log n) |
| [洛谷 P3902 递增](https://www.luogu.com.cn/problem/P3902)（顺带） | LIS | 模板 |
| [洛谷 P2758 编辑距离](https://www.luogu.com.cn/problem/P2758) | 编辑距离 | 模板直出 |
| UVA 11584 | 最少回文分割 | 12.5.2 |
| [洛谷 P4052 [JSOI2007] 文本生成器](https://www.luogu.com.cn/problem/P4052) | AC+计数 DP | 补集思想 |
| POJ 3691 / HDU 2457 | AC+最小修改 | DNA repair |
| [洛谷 P3975 [TJOI2015] 弦论](https://www.luogu.com.cn/problem/P3975) | SAM 路径 DP | 第 k 小 |
| SPOJ DQUERY（变体）/ CF 各种 | 子序列计数 | 12.6.4 公式 |
| [洛谷 P2292 [HNOI2004] L语言](https://www.luogu.com.cn/problem/P2292) | Trie + DP | 词典匹配可行性 |
| [洛谷 P2463 [SDOI2008] Sandy的卡片](https://www.luogu.com.cn/problem/P2463)（进阶） | 差分 + 公共子串 | 数串处理技巧 |

---

## 12.8 小结：字符串 DP 检查清单

1. 先问自己：**约束是"子串"还是"子序列"**？（子串→哈希/KMP/AC/SAM/Manacher 辅助判定；子序列→LCS/子序列自动机）
2. 是"匹配类约束"吗？→ 把**自动机状态**（KMP 前缀长度 / AC 节点 / SAM 节点）放进 DP 下标；
3. 是"区间划分类"吗？→ 想清楚 `s[j..i]` 满足的条件能否 O(1) 判定（预处理 01 表 / 哈希 / 回文半径）；
4. 关注**数据范围**决定 DP 维度与是否要优化（滚动数组、单调性、bitset、矩阵快速幂、slink）。

> 下一章（13）介绍字符串问题的另一类"非主流但高价值"工具：**bitset / Shift-And / FFT（NTT）**——位并行与卷积把匹配做成 O(n²/w) 或 O(n log n)。
