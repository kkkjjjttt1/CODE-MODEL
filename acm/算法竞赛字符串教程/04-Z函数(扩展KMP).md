# 第 04 章 Z 函数（Z-algorithm）与扩展 KMP

> **一句话**：`z[i]` = 串 `s` 的**每个后缀** `s[i..]` 与**整个串** `s` 的最长公共前缀（LCP）。
> Z 函数与 KMP 是"同一枚硬币的两面"：KMP 关注"已匹配前缀能退到哪"，Z 关注"每个后缀与全串有多像"。
> 本章 0 基下标；`z[0]` 习惯上记为 0（也可记 n，各题模板不同，抄代码时留意）。

---

## 4.1 定义

对 `s`（长 n），定义 `z[i] = LCP(s, s[i..])`，即最大 `k` 使得 `s[0..k-1] == s[i..i+k-1]`。

例：`s = "aabcaabxaaaz"`：
```
i:   0  1  2  3  4  5  6  7  8  9 10 11
s:   a  a  b  c  a  a  b  x  a  a  a  z
z:   -  1  0  0  3  1  0  0  2  1  0  0
```
（`z[4]=3`：`s[4..6]="aab"` 与全串前缀 `"aab"` 相同。）

---

## 4.2 线性算法：维护"匹配盒" `[l, r)`

核心思想：维护当前**最靠右的已匹配段** `[l, r)`（`s[l..r-1] == s[0..r-l-1]`）。当计算 `z[i]` 时，若 `i ∈ [l, r)`，则 `s[i..]` 与 `s[i-l..]` 开头相同，`z[i]` **至少**为 `min(r-i, z[i-l])`，直接以此为初值再暴力扩展，保证总复杂度 O(n)。

```cpp
vector<int> zFunction(const string& s) {
    int n = (int)s.size();
    vector<int> z(n, 0);
    int l = 0, r = 0;                       // 匹配盒 [l, r)，初始为空
    for (int i = 1; i < n; ++i) {
        if (i < r) z[i] = min(r - i, z[i - l]);   // 借用盒内信息
        while (i + z[i] < n && s[z[i]] == s[i + z[i]]) ++z[i];  // 暴力扩展
        if (i + z[i] > r) { l = i; r = i + z[i]; }               // 更新最右盒
    }
    return z;
}
```

**复杂度**：每次 while 扩展都会让 `r` 增大，`r` 至多到 n，所以总 O(n)。

> 与 KMP 的对比：KMP 是"失配时用 border 回退模式指针"，适合**在线**（边读文本边匹配）；Z 需要整串在手（**离线**），但写起来更直观，尤其适合"每个后缀与全串的 LCP"这类批量问题。

---

## 4.3 Z 函数的经典应用

### 4.3.1 字符串匹配（找所有出现位置）

`T = p + '#' + s`，求 T 的 z 数组，则 `z[i] == m` 的位置 i 说明 p 完整出现（同理 KMP 一章的做法）：

```cpp
vector<int> matchPositions(const string& s, const string& p) {
    int m = (int)p.size();
    string t = p + '\1' + s;               // 分隔符：任意不在字符集中的字符（'\1' 等）
    vector<int> z = zFunction(t), pos;
    for (int i = m + 1; i < (int)t.size(); ++i)
        if (z[i] >= m) pos.push_back(i - m - 1);   // 转回 s 的下标
    return pos;
}
```

> `z[i] >= m` 用 `>=` 而非 `==`：因为后面紧跟的字符若相同，z 会更大，但匹配盒已经能确认"长度 ≥ m"。

### 4.3.2 求最小周期 / 判整循环节

若 `z[p] >= n - p`，说明 `s[p..n-1] == s[0..n-p-1]`，即 **p 是一个周期**。因此：

```cpp
// 最小周期 p（p<n）：从小到大第一个满足 z[p] >= n-p 的 p
int minPeriodByZ(const string& s) {
    int n = (int)s.size();
    vector<int> z = zFunction(s);
    for (int p = 1; p < n; ++p)
        if (z[p] >= n - p) return p;
    return n;   // 无真周期
}
// 若题目要求 s 由 block 完整重复（整周期），再加条件 n % p == 0
```

### 4.3.3 压缩 / 去重："该串是否完全由某个前缀重复而成"

`n % p == 0 && z[p] >= n - p` ⇔ `s` 是 `s[0..p-1]` 的重复。数据压缩题常用（如"字符串压缩到最短" → 取最小的整周期 p）。

### 4.3.4 每个前缀的出现次数

用 z 数组可以统计**模式串每个前缀**在 s 中作为子串出现多少次？不能直接做——这类问题请用失配树（03 章 3.6），或后缀自动机（10 章）。

---

## 4.4 扩展 KMP（exKMP）：两串版本的 Z

### 4.4.1 问题

给定文本 `a`（长 n）与模式 `b`（长 m），对 a 的每个后缀求与 b 的 LCP：

```
ex[i] = LCP(b, a[i..])    (0 ≤ i < n)
```

例：`a = "ababc", b = "ab"` → `ex = [2, 0, 2, 0, 0]`。

### 4.4.2 模板（无需分隔符，任何字符集通用）

思路与单串 Z 完全一样，只是"参照串"换成 b 自己的 z 数组：

```cpp
// z: 模式 b 自身的 z 函数；ex[i]: a 的后缀 i 与 b 的 LCP
vector<int> exKMP(const string& a, const string& b) {
    int n = (int)a.size(), m = (int)b.size();
    vector<int> z(m, 0), ex(n, 0);

    // (1) 先求 b 自己的 z（参照 4.2 的单串代码，可复用 zFunction(b)）
    z = zFunction(b);

    // (2) 求 a 的每个后缀与 b 的 LCP
    int l = 0, r = 0;                    // 匹配盒 [l, r)：a[l..r-1] == b[0..r-l-1]
    for (int i = 0; i < n; ++i) {
        if (i < r) ex[i] = min(r - i, z[i - l]);
        while (i + ex[i] < n && ex[i] < m && b[ex[i]] == a[i + ex[i]]) ++ex[i];
        if (i + ex[i] > r) { l = i; r = i + ex[i]; }
    }
    return ex;
}
```

**要点**：
- `while` 里多一个 `ex[i] < m` 的限制：`a[i..]` 与 `b` 的 LCP 不可能超过 m；
- 与单串 z 唯一的区别是：盒内借用的是 `z[i-l]`（b 自己的 z），而不是"ex 自己"。

**简易平替**：如果不介意分隔符，`t = b + '\1' + a` 再跑单串 z 也能得到 ex（`z[i+m+1]`），代码少一半；字符集"所有可打印字符"的题用上面无分隔符版更稳。

### 4.4.3 应用

- **两串所有后缀的 LCP**：批量比较两串的后缀字典序、找"文本中每个位置能匹配模式多长"（如病毒/生物序列的拼接比对题）；
- **最长公共子串**：把两串拼起来跑 SA（09 章）也行；用 exKMP 需配合二分/分治；
- 计算 **B 的每个前缀在 A 中出现多少次**：出现位置 i 需 `ex[i] >= len`，统计可用差分（`diff[1]++, diff[ex[i]+1]--` 之类按长度贡献），见例题 P5410 的思路。

---

## 4.5 经典例题

| 题目 | 考点 | 要点 |
|------|------|------|
| [洛谷 P5410 【模板】扩展 KMP（Z 函数）](https://www.luogu.com.cn/problem/P5410) | z + ex | 输出 `z[1..m-1]`（b 的 z）与 `ex[0..n-1]`（a 各后缀与 b 的 LCP） |
| [CF 126B Password](https://codeforces.com/problemset/problem/126/B) | z 求 border + 中部出现 | 候选 = 满足 `z[i] == n-i` 的 border；再检查中间是否存在 ≥ 该长度的 z 值 |
| [洛谷 P4391 Radio Transmission](https://www.luogu.com.cn/problem/P4391)（复习） | 最小周期 | KMP / Z 均可：最小 p 使 `z[p] ≥ n-p` |
| [洛谷 P3503 [POI2010] Blocks](https://www.luogu.com.cn/problem/P3503)（进阶） | 单调栈 | 与 z 无关，属"字符串+数据结构"综合，适合后期练习 |

> CF 126B 提示：先用 z 找出所有"既是前缀又是后缀"的长度（`z[i] = n - i`），记为 border 集合；若某个 border 长度 L 在**中间的某个 z[j] ≥ L**（`L ≤ j ≤ n-L-1` 之外的位置也行，需 j 不在两端）中存在，即答案为该 L，取最大者；否则输出 *Just a legend*。

---

## 4.6 Z 与 KMP 的关系小结

| | KMP（前缀函数 pi） | Z 函数 |
|---|---|---|
| 求什么 | 每个前缀的**最长 border** | 每个后缀与全串的 **LCP** |
| 语义侧重 | "前缀内部" | "后缀 vs 全串" |
| 在线？ | 可在线匹配 | 需整串 |
| 求周期 | `n - pi[n-1]` | 最小 p：`z[p] ≥ n-p` |
| 匹配模式 | 需拼接分隔符或经典 next | 需拼接分隔符 |
| 两串 LCP | 不方便 | **exKMP 直接支持** |

> 还有一条深层联系：`pi` 与 `z` 可通过"反串"互相转化（`z` 数组倒序对应反串的 pi 链），竞赛中几乎不需要手动互转，两者都会写即可。

---

## 4.7 练习题建议（由易到难）

1. P5410（模板，务必 AC）；
2. CF 126B（border + 出现性判断，训练 z 的灵活运用）；
3. 用 z 重做 P3375（对照 KMP 写法）；
4. 尝试："给定 s，判断它能否由某循环节重复得到，并输出最小循环节长度"（z 版与 KMP 版各写一遍）。

> 下一章（05）进入**树形结构**的第一站：Trie 字典树——AC 自动机、可持久化结构的地基。
