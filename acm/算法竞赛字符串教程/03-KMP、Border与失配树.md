# 第 03 章 KMP、Border 与失配树

> **一句话**：KMP 用"已匹配部分的信息"避免文本指针回退，把单模式匹配从 O(nm) 降到 **O(n + m)**。
> 本章同时讲透 **border / 周期** 理论与 **失配树（fail 树）**——它们是把 KMP 从"模板"升级为"武器"的关键。
> 本章代码**默认 0 基下标**，字符串是 `std::string`。

---

## 3.1 朴素匹配为什么慢

文本 `s`（长 n）、模式 `p`（长 m）逐位对齐比较，失败时文本指针 i 和模式指针 j **都回退**，最坏 O(nm)（如 `s="aaaa...b", p="aaaa...a"`）。

KMP 的洞察：失配时，**`p[0..j-1]` 已经与 `s[i-j..i-1]` 匹配**，这段已匹配文本可能已包含模式串的某个"真前缀"，因此 j 只需回退到该前缀长度，i **永不回退**。

---

## 3.2 前缀函数 `pi`（Prefix Function）

### 定义

对模式串（或任意串）`p`，`pi[i]` = 前缀 `p[0..i]` 的**最长真 border 长度**，即最大的 `k < i+1`，使 `p[0..k-1] == p[i-k+1..i]`。

例：`p = "aabaaab"`：
```
p:     a a b a a a b
pi:    0 1 0 1 2 2 3
```
（`pi[6]=3`，因为 `"aab"` 既是前缀也是后缀 `p[4..6]`。）

### O(n) 算法（核心代码，务必背熟）

```cpp
vector<int> prefixFunction(const string& p) {
    int n = (int)p.size();
    vector<int> pi(n, 0);
    for (int i = 1; i < n; ++i) {
        int j = pi[i - 1];                 // 上一位置的 border 长度（候选起点）
        while (j > 0 && p[i] != p[j])      // 失配：退到 p[j-1] 的 border
            j = pi[j - 1];
        if (p[i] == p[j]) ++j;             // 匹配成功则长度 +1
        pi[i] = j;
    }
    return pi;
}
```

**正确性要点**（为什么从 `pi[i-1]` 开始、为什么失配跳 `pi[j-1]`）：
- 若 `p[0..i]` 有长度为 k 的 border，则去掉最后一个字符后，`p[0..i-1]` 有长度为 k−1 的 border；
- 因此新 border 长度 k 只能来自 `p[0..i-1]` 的 border 链（`pi[i-1]` → `pi[pi[i-1]-1]` → …）逐项试探，第一条匹配成功的链就是最长的。
- 均摊 O(n)：j 每次 while 至少减 1，而每条边最多回退一次。

---

## 3.3 用前缀函数做匹配

把模式与文本用**一个分隔符 `#`** 拼起来：`T = p + '#' + s`，跑前缀函数。则 `pi` 值等于 m 的位置，说明该处模式完整出现（`#` 保证不会跨段匹配）：

```cpp
// 返回 p 在 s 中所有出现位置的起始下标（0 基）
vector<int> kmpMatch(const string& s, const string& p) {
    int m = (int)p.size();
    string t = p + "#" + s;            // 分隔符需不在字符集中
    vector<int> pi = prefixFunction(t), pos;
    for (int i = m + 1; i < (int)t.size(); ++i)
        if (pi[i] == m) pos.push_back(i - 2 * m);
    return pos;
}
// 注意：i-2*m 的推导：匹配结束于 t 的下标 i，模式起点 = i - m，再减去前导 "p#"(长 m+1)
```

复杂度 O(n+m)。文本指针从不回退，这一性质在"在线流式匹配"中也有用。

---

## 3.4 经典 KMP 的 `next` 数组写法（中文教材版）

很多题解用 1 基前缀长度的 `nxt`（`nxt[i]` = 前缀 `p[0..i-1]` 的最长 border 长度，即 `pi[i-1]`），配一个 `-1` 哨兵。两者完全等价，**认准一套下标别混**：

```cpp
const int N = 1000005;
int nxt[N];                                   // nxt[i]: 前缀长度 i 的最长 border 长度

void buildNext(const string& p) {
    nxt[0] = -1;                              // 哨兵
    int i = 0, j = -1, m = (int)p.size();
    while (i < m) {
        if (j == -1 || p[i] == p[j]) nxt[++i] = ++j;
        else j = nxt[j];
    }
}

vector<int> kmpSearch(const string& s, const string& p) {
    buildNext(p);
    vector<int> pos;
    int n = (int)s.size(), m = (int)p.size(), i = 0, j = 0;
    while (i < n) {
        if (j == -1 || s[i] == p[j]) { ++i; ++j; }
        else j = nxt[j];
        if (j == m) {                         // 找到一次（位置 i-m）
            pos.push_back(i - m);
            j = nxt[j];                       // 允许重叠继续找
        }
    }
    return pos;
}
```

| 写法 | 下标含义 | 跳转 |
|------|---------|------|
| `pi[i]`（0 基字符下标） | `p[0..i]` 的 border | 失配跳 `pi[j-1]` |
| `nxt[i]`（1 基长度） | `p[0..i-1]` 的 border | 失配跳 `nxt[j]` |

**两数组关系**：`nxt[i] == pi[i-1]`（`i≥1`）。

> 禁止重叠的匹配：找到一次后把 `j` 置 0 而不是 `nxt[j]`；或统计时若 `pos[i] - last < m` 则跳过。

---

## 3.5 Border 与周期：KMP 的灵魂结论

### 3.5.1 最小周期

> **定理**：串 `s`（长 n）的**最小周期**为 `p = n - pi[n-1]`（也即 `n - nxt[n]`）。
> 若 `n % p == 0`，则 `s` 是循环节 `s[0..p-1]` 的完整重复（整周期）；否则周期仍然存在（如 `"abababa"` 最小周期 2，但 7 不是 2 的倍数）。

直觉：border 是"前后重叠相等"，把 border 重叠掉之后剩下的部分就是周期。

```cpp
// 求 s 的最小周期（P4391 无线传输 直接输出 n - pi[n-1]）
int minPeriod(const string& s) {
    auto pi = prefixFunction(s);
    return (int)s.size() - pi.back();
}
```

### 3.5.2 全部 border

`pi[n-1]` 给出最长 border；其 border 又是原串的 border……迭代下去得到**所有 border 长度**：

```cpp
vector<int> allBorders(const vector<int>& pi) {
    vector<int> res;
    int j = pi.back();               // 从最长 border（长度 pi[n-1]）开始
    while (j > 0) {
        res.push_back(j);            // 记录一个 border 长度
        j = pi[j - 1];               // 跳向更短的 border
    }
    return res;                      // 降序返回所有 border 长度
}
```

例：`"aaaa"` → border 链 3,2,1。

**周期与 border 对偶**：长度 k 是 border ⇔ `n-k` 是周期（`k<n`）。所以"找最长 border"= "找最短周期"。

---

## 3.6 失配树（fail 树 / KMP 树）

### 3.6.1 构造

把"前缀长度"看作节点（`0..n`，节点 0 是空串/根），对每个 `i ∈ [1, n]` 连边 `nxt[i] → i`（或 `pi[i-1] → i`）：

```cpp
// 建失配树：fa[i] = nxt[i]（i = 1..n），fa[0] = 0
// 边 fa[i] -> i；fa[i] < i，所以这是一棵以 0 为根的树
```

**重要性质**：
1. 节点 `u` 到根的路径上依次是前缀 `u` 的**全部 border**（长度递减）；
2. 节点 `u` 的**子树** = 所有"以 `u` 为 border"的前缀 → 在失配树上做子树统计，可以得到 **border 出现次数 / 前缀出现次数**。

### 3.6.2 应用一：统计每个模式（=某前缀）出现次数

对整串 `s` 建失配树，每个节点 `i(1..n)` 初始权值 1（代表"文本在位置 i 结束"这一事件）。**按 i 从大到小**把 `cnt[fa[i]] += cnt[i]`（等价于拓扑/按深度），最后：

```
cnt[u] = 前缀 u（即子串 s[0..u-1]）在 s 中出现的次数
```

```cpp
vector<int> cnt(n + 1, 1);              // 节点 0 也初始为 1？不，0 表示空串不参与
for (int i = 1; i <= n; ++i) cnt[i] = 1;
for (int i = n; i >= 1; --i) cnt[nxt[i]] += cnt[i];   // 前提 nxt[i] < i 且按 i 降序即拓扑序
// 之后 cnt[L] = 子串 s[0..L-1] 的出现次数
```

例：`s = "ababa"`，前缀 `"aba"`（u=3）出现 2 次（位置 0、2）→ `cnt[3] = 2`。✓

### 3.6.3 应用二：失配树上的 LCA

两前缀的**最长公共 border** = 它们在失配树上的 **LCA**（若 LCA 恰为某个端点，则答案取该端点的父亲——因为 border 必须是真前缀）。配合倍增 LCA，可 O(n log n) 预处理、O(log n) 回答。

> 例题 [洛谷 P5829 【模板】失配树](https://www.luogu.com.cn/problem/P5829)：查询前缀 a 与 b 的最长公共 border。

---

## 3.7 KMP 自动机（用于 DP 的视角）

匹配过程可以看作一个**自动机**：状态 = "当前已匹配的模式前缀长度 j"（0..m），读入字符 c 后转移到新状态：

```
while (j > 0 && p[j] != c) j = nxt[j];   // 沿 border 链退
if (p[j] == c) ++j;
if (j == m) { 命中一次; j = nxt[j]; }     // m 状态是"终止/命中"
```

若把转移表预处理成 `go[j][c]`（O(m·|Σ|)），就得到完整的 DFA。之后"生成长度为 n 的串、不能包含模式 p"这类计数/最优化问题，只需在这个自动机上做 DP（详见 12 章）。例如"长度为 n、不含子串 p 的 01 串个数"，可以进一步用**矩阵快速幂**加速。

---

## 3.8 本章模板汇总（可直接抄）

```cpp
#include <bits/stdc++.h>
using namespace std;

vector<int> prefixFunction(const string& p);  // 见 3.2

int main() {  // 例：读入文本与模式，输出所有匹配位置
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string s, p;
    cin >> s >> p;
    auto pos = kmpMatch(s, p);                // 见 3.3
    for (int x : pos) cout << x << '\n';
    auto pi = prefixFunction(p);
    for (int i = 0; i < (int)p.size(); ++i) cout << pi[i] << " \n"[i + 1 == (int)p.size()];
    return 0;
}
```

---

## 3.9 经典例题

| 题目 | 考点 | 要点 |
|------|------|------|
| [洛谷 P3375 【模板】KMP](https://www.luogu.com.cn/problem/P3375) | 匹配 + border | 输出所有出现位置与模式串的 `nxt/pi` |
| [洛谷 P4391 [BOI2009] Radio Transmission](https://www.luogu.com.cn/problem/P4391) | 最小周期 | 答案 = n − pi[n−1]（题目保证 s 由某循环节重复/截断而来） |
| [CF 432D Prefixes and Suffixes](https://codeforces.com/problemset/problem/432/D) | border + 出现次数 | 枚举 border 链上的长度，失配树子树和求出现次数 |
| [洛谷 P5829 【模板】失配树](https://www.luogu.com.cn/problem/P5829) | 失配树 + LCA | 公共 border = fail 树 LCA，注意端点取父 |
| [洛谷 P2375 [NOI2014] 动物园](https://www.luogu.com.cn/problem/P2375) | 不重叠 border 计数 | 双指针：主循环维护 j，另用 `j2` 限制 `2*j2 ≤ i`，答案与 border 链深度有关 |
| [洛谷 P3435 [POI2006] OKR-Periods of Words](https://www.luogu.com.cn/problem/P3435) | 周期 + 路径压缩 | 每个前缀贡献 = 长度 − 其"最小 border 链上的最短 border"，记忆化跳 next 到最靠近根的 border |
| [洛谷 P3426 [POI2005] SZA-Template](https://www.luogu.com.cn/problem/P3426)（进阶） | border 链 + DP | 印章问题，需维护 border 集合的若干性质 |

---

## 3.10 总结：KMP 全家桶关系图

```
prefixFunction(pi)  --等价-->  nxt(经典版)
        │
        ├── KMP 匹配（O(n+m)，文本指针不回退）
        ├── 最小周期 = n - pi[n-1]
        ├── border 链 = 反复 pi[j-1]
        └── 失配树 = 树边 nxt[i]→i
              ├── 子树和 → 出现次数
              └── LCA → 最长公共 border
```

> 下一章（04）的 **Z 函数**与 KMP 殊途同归：它求"每个后缀与整个串的 LCP"，在求 LCP、压缩、最小周期等场景有时更好写。
