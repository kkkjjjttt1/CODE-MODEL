# 第 07 章 Manacher 与回文专题

> **一句话**：Manacher 在 **O(n)** 内求出以每个位置为中心的最长回文半径，是"最长回文子串 / 回文计数 / 双回文"一类题的标配。
> 本章 0 基下标；先讲 O(n²) 中心扩展建立直觉，再讲线性算法与它的经典衍生应用。

---

## 7.1 回文与"中心扩展"直觉

回文（如 `"ababa"`、`"abba"`）有一个独有特征：**从中心向两边对称相等**。
- 奇数长回文中心是一个字符：`"aba"` 中心 `b`；
- 偶数长回文中心是两个字符之间的空隙：`"abba"` 中心在 `b|b` 之间。

朴素方法：枚举每个"中心"（n 个字符中心 + n−1 个空隙中心），向两边扩展比较，O(n²)。

例：`s = "abababa"`，以第 3 个字符 `b` 为中心向外扩 → 半径 4，得最长回文 `"abababa"`（7 个字符）。

---

## 7.2 Manacher 线性算法

### 7.2.1 思路

维护**最靠右的回文区间** `[l, r]`（当前所有中心中能覆盖到的最右边界）。处理新中心 i 时：
- 若 `i ≤ r`，则 i 关于当前区间中心（设为 c）的**镜像点** `m = l + r - i` 已在区间内；以 m 为中心的回文半径可以**直接复制**给 i（在区间限制内），避免重复比较 → 这就是线性性的来源；
- 然后继续暴力扩展（扩出去的字符每个最多被比较常数次，因为 r 单调右移）。

### 7.2.2 模板：奇数回文半径 `d1`、偶数回文半径 `d2`

**推荐直接记这两段**（无需在串里插 `#`，最不容易错）：

```cpp
// d1[i]：以 i 为中心的最长奇回文半径（含中心），回文为 s[i-d1[i]+1 .. i+d1[i]-1]
vector<int> manacherOdd(const string& s) {
    int n = (int)s.size();
    vector<int> d1(n);
    for (int i = 0, l = 0, r = -1; i < n; ++i) {
        int k = (i > r) ? 1 : min(d1[l + r - i], r - i + 1);
        while (i - k >= 0 && i + k < n && s[i - k] == s[i + k]) ++k;
        d1[i] = k--;
        if (i + k > r) { l = i - k; r = i + k; }
    }
    return d1;
}

// d2[i]：以 i-1 与 i 之间空隙为中心的最长偶回文半径，回文为 s[i-d2[i] .. i+d2[i]-1]
vector<int> manacherEven(const string& s) {
    int n = (int)s.size();
    vector<int> d2(n);
    for (int i = 0, l = 0, r = -1; i < n; ++i) {
        int k = (i > r) ? 0 : min(d2[l + r - i + 1], r - i + 1);
        while (i - k - 1 >= 0 && i + k < n && s[i - k - 1] == s[i + k]) ++k;
        d2[i] = k--;
        if (i + k > r) { l = i - k - 1; r = i + k; }
    }
    return d2;
}
```

**另一种常见写法（插 `#` 统一奇偶）**：构造 `t = "#a#b#b#a#"`，对 t 跑奇数版 Manacher 得半径 `rad[i]`（含中心）。则：
- 原串的回文长度 = `rad[i] - 1`（中心在 `#` 上时是偶数长回文，在字母上时是奇数长回文）；
- 原串中心位置 = `(i - 1) / 2`（近似）。

两种写法等价，推荐掌握 d1/d2 版。

> **正确性验证小数据**：`s="abba"` → `d1=[1,1,1,1]`（无奇回文长于1）；`d2=[0,0,2,0]`？算一下：i=2（空隙在 s[1],s[2] 之间）扩展得 "abba" 半径 2 ✓，故最长回文 = 4。
> `s="aba"` → `d1=[1,2,1]`，最长 = 2×2−1 = 3 ✓。

---

## 7.3 基本应用

### 7.3.1 最长回文子串长度

```cpp
int longestPalindromeLen(const string& s) {
    auto d1 = manacherOdd(s), d2 = manacherEven(s);
    int n = (int)s.size(), ans = 1;
    for (int i = 0; i < n; ++i) {
        ans = max(ans, 2 * d1[i] - 1);   // 奇
        ans = max(ans, 2 * d2[i]);       // 偶
    }
    return ans;
}
```

### 7.3.2 回文子串总数（按中心计数）

每个中心 i 的每个半径 k（1..d1[i]）都对应一个**不同的**奇回文（中心不同或半径不同 → 子串不同）。同理偶中心。所以：

```cpp
long long countPalindromicSubstrings(const string& s) {
    auto d1 = manacherOdd(s), d2 = manacherEven(s);
    long long ans = 0;
    for (int i = 0; i < (int)s.size(); ++i) { ans += d1[i] + d2[i]; }
    return ans;   // 例："aaa" → d1=[1,2,1], d2=[0,1,0] → 6 个回文子串（"a","a","a","aa","aa","aaa"）
}
```

### 7.3.3 O(1) 判断任意区间是否回文

```cpp
// s[l..r] 是否回文（0 基闭区间）
bool isPalByRadius(const vector<int>& d1, const vector<int>& d2, int l, int r) {
    int len = r - l + 1;
    if (len & 1) {
        int c = (l + r) >> 1;                 // 奇回文中心
        return d1[c] >= (len + 1) / 2;
    } else {
        int c = (l + r + 1) >> 1;             // 偶回文中心（空隙在 c-1 与 c 之间）
        return d2[c] >= len / 2;
    }
}
```

### 7.3.4 按长度统计回文个数（差分）

例如"长度为奇数 L 的回文子串有多少个"（[P1659 拉拉队排练](https://www.luogu.com.cn/problem/P1659) 等）：中心 i 的半径 d1[i] 说明它贡献了**每条奇数长度** `1,3,...,2*d1[i]-1` 一个回文。用差分数组在"半径桶"上打点，再前缀和得到每个半径的出现次数：

```cpp
// 需求：统计"长度为 2t-1 的奇数回文子串"各有多少个
vector<int> diff(n + 2, 0);
for (int i = 0; i < n; ++i) { diff[1]++; diff[d1[i] + 1]--; }
for (int t = 1; t <= n; ++t) diff[t] += diff[t - 1];
// 前缀和之后 diff[t] = 半径 ≥ t 的中心个数。
// 半径 d1[i] 的中心恰好贡献长度 1,3,…,2*d1[i]-1 各一个回文，
// 所以"长度恰为 2t-1 的回文个数" = diff[t]。
```

（更常用的写法：`cnt[1]++, cnt[d1[i]+1]--` 后求前缀和，得到"半径 ≥ t 的中心总数"；而"长度恰为 2t-1 的回文数"= 半径≥t 的中心数总和变化量的累计，按题目口径处理。）

---

## 7.4 进阶应用：最长双回文串（[洛谷 P4555](https://www.luogu.com.cn/problem/P4555)）

**问题**：求最长子串 `s[i..j]`，可被切成两段**均为回文**（两段都不为空）。

**做法**：
1. 预处理 `lenEnd[x]` = 以 x 结尾的最长回文长度；`lenSt[x]` = 以 x 开头的最长回文长度；
2. 答案 = `max_x ( lenEnd[x] + lenSt[x+1] )`。

**如何 O(n) 求 lenEnd / lenSt**：每个回文中心 i 贡献一串"端点值线性增长"的回文。关键观察：右端点序列随 x 增大，**活跃中心的最小下标单调不减**，于是可以用队列维护（先进先出 + 过期弹出），总 O(n)：

```cpp
#include <bits/stdc++.h>
using namespace std;

vector<int> manacherOdd(const string& s);   // 同 7.2.2
vector<int> manacherEven(const string& s);

// lenEnd[i]：以 i 结尾的最长回文长度；lenSt[i]：以 i 开头的最长回文长度
void longestPalEndStart(const string& s, vector<int>& lenEnd, vector<int>& lenSt) {
    int n = (int)s.size();
    auto d1 = manacherOdd(s);
    auto d2 = manacherEven(s);
    lenEnd.assign(n, 0);
    lenSt.assign(n, 0);

    // ---- 以 x 结尾 ----
    queue<int> q;
    for (int x = 0; x < n; ++x) {                 // 奇中心：覆盖右端点 [i, i+d1[i]-1]
        q.push(x);
        while (!q.empty() && q.front() + d1[q.front()] - 1 < x) q.pop();
        if (!q.empty()) lenEnd[x] = max(lenEnd[x], 2 * (x - q.front()) + 1);
    }
    while (!q.empty()) q.pop();
    for (int x = 0; x < n; ++x) {                 // 偶中心（空隙 i-1|i）：覆盖右端点 [i, i+d2[i]-1]
        q.push(x);
        while (!q.empty() && q.front() + d2[q.front()] - 1 < x) q.pop();
        if (!q.empty()) lenEnd[x] = max(lenEnd[x], 2 * (x - q.front() + 1));
    }

    // ---- 以 x 开头（对称地从右往左扫）----
    while (!q.empty()) q.pop();
    for (int x = n - 1; x >= 0; --x) {            // 奇中心：覆盖左端点 [i-d1[i]+1, i]
        q.push(x);
        while (!q.empty() && q.front() - d1[q.front()] + 1 > x) q.pop();
        if (!q.empty()) lenSt[x] = max(lenSt[x], 2 * (q.front() - x) + 1);
    }
    while (!q.empty()) q.pop();
    for (int x = n - 1; x >= 0; --x) {            // 偶中心 i=x+1 的覆盖起点区间 [i-d2[i], i-1]
        int i = x + 1;
        if (i < n) q.push(i);
        while (!q.empty() && q.front() - d2[q.front()] > x) q.pop();
        if (!q.empty()) lenSt[x] = max(lenSt[x], 2 * (q.front() - x));
    }
}

int main() {   // P4555：最长双回文串
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string s; cin >> s;
    vector<int> lenEnd, lenSt;
    longestPalEndStart(s, lenEnd, lenSt);
    int n = (int)s.size(), ans = 0;
    for (int x = 0; x + 1 < n; ++x)
        ans = max(ans, lenEnd[x] + lenSt[x + 1]);
    cout << ans << '\n';
    return 0;
}
```

> 队列正确性直观说明：当 x 向右推进，中心 i 一旦"够不到"x 就再也不会回来（右端点单调），所以先进先出的中心队列 + 队首过期弹出正好维护"仍然覆盖 x 的最小/最大下标"。偶/奇、左/右四个方向只是公式里的常数不同。

---

## 7.5 其他衍生题型

- **相交回文对计数**（CF 17E *Palisection*）：总数 C(回文数, 2) − 不相交对数；不相交 = 枚举分界，左边"以 ≤i 结尾的回文数"× 右边"以 >i 开头的回文数"（用 7.3.2 的差分统计端点）。
- **回文子串区间覆盖 / 最少回文分割**：Manacher 得到若干"回文区间"，问题转化为区间覆盖/最短路 DP（另见 12 章 DP 部分）。
- **本质不同回文子串 / 回文出现次数**：交给**回文自动机 PAM**（下一章）。
- **哈希二分判回文**（回顾 02 章）：O(n log n)，代码简单但慢于 Manacher，适用于"需要在线/带修改"的变体。

---

## 7.6 经典例题

| 题目 | 考点 | 要点 |
|------|------|------|
| [洛谷 P3805 【模板】Manacher](https://www.luogu.com.cn/problem/P3805) | 最长回文 | d1/d2 模板直出 |
| [洛谷 P4555 最长双回文串](https://www.luogu.com.cn/problem/P4555) | lenEnd + lenSt | 7.4 的队列法 |
| [洛谷 P1659 [国家集训队] 拉拉队排练](https://www.luogu.com.cn/problem/P1659) | 按长度统计回文 + 快速幂 | 差分统计每种奇数长度的回文个数，从大到小取前 k 个乘起来 |
| CF 17E *Palisection* | 回文计数 + 容斥 | 差分算端点数量，总 − 不相交 |
| [洛谷 P4287 [SHOI2011] 双倍回文](https://www.luogu.com.cn/problem/P4287)（进阶） | 回文+border 思想 | 需沿回文半径做额外判定，可用 PAM/Manacher+枚举 |

---

## 7.7 小结

```
中心扩展 O(n²) ──重复比较太多
Manacher：复用"镜像点半径"，维护最右回文区间 [l, r]
  d1[i]：奇回文半径（含中心）
  d2[i]：偶回文半径（空隙 i-1|i）
复杂度 O(n)，线性求出：
  · 最长回文子串        max(2*d1[i]-1, 2*d2[i])
  · 回文子串总数        Σ d1[i] + Σ d2[i]
  · 区间是否回文        O(1)（见 7.3.3）
  · lenEnd/lenSt        （队列法，双回文等题）
  · 各长度回文个数      差分
本质不同的回文 → 第 08 章 回文自动机 PAM
```

> 下一章（08）用**回文自动机（Eertree / PAM）**回答 Manacher 答不了的问题：本质不同回文子串的种类、每种的出现次数。
