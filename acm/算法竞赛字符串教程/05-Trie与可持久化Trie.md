# 第 05 章 Trie 字典树（含 01-Trie 与可持久化 Trie）

> **一句话**：Trie 把一组字符串按**公共前缀**组织成树——从根到某个节点的路径就是一个前缀。
> 它是 AC 自动机（06 章）的地基，本身也是"前缀统计 / 最大异或 / 可持久化查询"的利器。
> 本章下标 0 基；节点数组从 0（根）开始。

---

## 5.1 基本结构与插入

每个节点存"通向下一字符的指针 + 附加信息（计数）"。小写字母集用定长数组 `ch[26]` 最快：

```cpp
struct Trie {
    struct Node {
        int ch[26];      // 子节点编号，0 表示不存在
        int cnt;         // 经过该节点的字符串个数（前缀计数）
        int end;         // 以该节点结尾的字符串个数（模式串去重计数用）
        Node() { memset(ch, 0, sizeof ch); cnt = end = 0; }
    };
    vector<Node> tr;
    Trie() { tr.emplace_back(); }        // 根节点 0

    void insert(const string& s) {
        int u = 0;
        for (char c : s) {
            int x = c - 'a';
            if (!tr[u].ch[x]) {
                tr[u].ch[x] = (int)tr.size();
                tr.emplace_back();
            }
            u = tr[u].ch[x];
            tr[u].cnt++;                 // 统计该前缀出现次数
        }
        tr[u].end++;                     // 完整插入一个串
    }

    // 查询 s 作为“前缀”出现多少次（即有多少个插入串以 s 开头）
    int queryPrefix(const string& s) const {
        int u = 0;
        for (char c : s) {
            u = tr[u].ch[c - 'a'];
            if (!u) return 0;
        }
        return tr[u].cnt;
    }

    // 查询 s 作为“完整串”被插入了多少次（去重/点名用）
    int queryWord(const string& s) const {
        int u = 0;
        for (char c : s) {
            u = tr[u].ch[c - 'a'];
            if (!u) return 0;
        }
        return tr[u].end;
    }
};
```

- 空间：节点数 ≤ **所有插入串的总长度 + 1**（有共享前缀时更少），一般开 `总长 + 5`；
- 时间：插入/查询 O(|s|)；
- 删除：竞赛中常用"逻辑删除"（`cnt--`）而不是物理释放，避免麻烦。

---

## 5.2 定长数组的局限与替代

- 字符集是 `a..z`（26）：`ch[26]` 数组，**最快**；
- 字符集是大写+数字（36/62）或更大（如所有可见 ASCII 约 95）：数组会浪费，改用：
  ```cpp
  // 每节点一个 map：内存 O(边数 log)，查询多一个 log
  struct Node2 { map<int, int> ch; int cnt = 0; };
  // 或把所有串排序去重后手动模拟（离线场景常用）
  ```
- 只在"读入全部串"后需要批量处理时，也可以**排序 + 相邻比较**替代一部分 Trie 功能——但 Trie 的"在线插入 + 前缀查询"不可替代。

---

## 5.3 经典应用

### 5.3.1 前缀统计 / 点名

- 插入 n 个名字，随后 m 次点名：`queryWord` 区分"没这个人 / 点到一次 / 重复点名"（[洛谷 P2580](https://www.luogu.com.cn/problem/P2580)）；
- 问"有多少单词以某前缀开头"：`queryPrefix`；
- 判断一个串**是否是一个集合中某个串的前缀 / 是否集合中存在一个串是它的前缀**：插入时若经过 `end>0` 的节点、或插入结束后发现新节点没增长（全被共享）即可判断。

### 5.3.2 字典序相关

Trie 的 DFS 输出即为**字典序升序**的所有串；根到节点 u 的路径 = 一组串的公共前缀，`cnt` 最大的深节点 ≈ 最长公共前缀。

### 5.3.3 自动补全 / 浏览器联想

查询以某前缀开头的所有串：DFS 该前缀节点下的子树，取 `end>0` 的节点。

---

## 5.4 01-Trie：处理"异或最大值"

把每个**整数看成 30/31 位二进制串**插入 Trie，查询时**贪心**：从高位到低位尽量走向与 x 相反的位（`x 的某位是 0 就走 1`），保证异或结果最大。每一位都取到最优 = 全局最优（贪心正确性来自高位权重大）。

### 5.4.1 最大异或对

```cpp
const int MAXB = 30;                      // a_i < 2^31
struct Trie01 {
    vector<array<int, 2>> ch;             // ch[u][0/1]
    Trie01() { ch.emplace_back(); ch[0] = {0, 0}; }

    void insert(int x) {
        int u = 0;
        for (int b = MAXB; b >= 0; --b) {
            int bit = (x >> b) & 1;
            if (!ch[u][bit]) {
                ch[u][bit] = (int)ch.size();
                ch.emplace_back(); ch.back() = {0, 0};
            }
            u = ch[u][bit];
        }
    }
    // 与 x 异或能得到的最大值
    int maxXor(int x) const {
        int u = 0, res = 0;
        for (int b = MAXB; b >= 0; --b) {
            int bit = (x >> b) & 1;
            if (ch[u][bit ^ 1]) {          // 优先走相反位
                res |= 1 << b;
                u = ch[u][bit ^ 1];
            } else u = ch[u][bit];
        }
        return res;
    }
};
// 最大异或对 = 依次插入每个 a[i]，查询 maxXor(a[i]) 取最大
```

复杂度 O(n·B)，B=30/31。

### 5.4.2 树上最大异或路径

经典性质：记 `d[u]` = 从根到 u 路径上所有边权的异或和，则 **u→v 路径异或 = d[u] xor d[v]**（LCA 以上部分异或两次抵消）。于是问题又变成"一堆数（所有 d）中选两个的最大异或" → 直接 01-Trie。

> 例题：[洛谷 P4551 最长异或路径](https://www.luogu.com.cn/problem/P4551)。

---

## 5.5 可持久化 Trie

### 5.5.1 思想

每次插入**新建一条从根到叶的链**（O(B) 个新节点），旧节点全部保留共享，`root[i]` 记下第 i 次插入后的根。两个版本相减即得"区间内有哪些数"——与可持久化线段树同理。

### 5.5.2 可持久化 01-Trie（查询 [l, r] 内与 x 的最大异或）

```cpp
struct PersistentTrie01 {
    struct Node { int ch[2], sz; };        // sz: 该子树中的数字个数
    vector<Node> tr;                       // tr[0] = 空节点
    vector<int> root;                      // root[i]: 前 i 个数的版本
    int MAXB = 30;

    PersistentTrie01() { tr.push_back({0, 0, 0}); root.push_back(0); }

    int insert(int pre, int x) {           // 在版本 pre 基础上插入 x
        int cur = (int)tr.size();
        tr.push_back(tr[pre]);             // 复制根
        tr[cur].sz++;
        int u = cur;
        for (int b = MAXB; b >= 0; --b) {
            int bit = (x >> b) & 1;
            int nxt = (int)tr.size();
            tr.push_back(tr[tr[u].ch[bit]]);   // 复制将要走的子节点
            tr[nxt].sz++;
            tr[u].ch[bit] = nxt;
            u = nxt;
        }
        return cur;
    }

    // 在 a[l..r]（1 基，版本 root[r] 减 root[l-1]）中找与 x 异或的最大值
    int query(int l, int r, int x) const {
        int u = root[r], v = root[l - 1], ans = 0;
        for (int b = MAXB; b >= 0; --b) {
            int want = ((x >> b) & 1) ^ 1;     // 想要 1 的位
            int cnt = tr[tr[u].ch[want]].sz - tr[tr[v].ch[want]].sz;
            if (cnt > 0) {                     // 区间内存在该分支
                ans |= (1 << b);
                u = tr[u].ch[want]; v = tr[v].ch[want];
            } else {
                u = tr[u].ch[want ^ 1]; v = tr[v].ch[want ^ 1];
            }
        }
        return ans;
    }
};
```

**例题 [洛谷 P4735 最大异或和](https://www.luogu.com.cn/problem/P4735)**：维护可追加的序列，操作 ① 末尾追加 x，② 询问 `max_{p∈[l-1, r-1]} (s[p] xor s[n] xor x)`，其中 s 为前缀异或数组。把每个 s[i] 插入可持久化 01-Trie，区间查询即上面 `query(l, r, s[n] xor x)`。（把 `s[0]=0` 也插入。）

- 空间：每次插入新建 MAXB+1 ≈ 31 个节点，n 次插入约 31·n 个，注意开够或 `reserve`；
- 普通字符 Trie 也可持久化（如"区间内某前缀是否存在"），写法相同，只是把 2 分支换成多分支。

---

## 5.6 番外：子序列自动机（判断子序列 / 子序列 DP）

与 Trie 无关但名字相近、且是常见小题：预处理 `nxt[i][c]` = 位置 i 之后（含 i）字符 c **第一次出现**的位置，则可以在 O(|t|·?) 内判断 t 是否是 s 的子序列，还能做"不同子序列个数"DP：

```cpp
// s 长 n，字符集小写
vector<array<int, 26>> buildSeqAuto(const string& s) {
    int n = (int)s.size();
    vector<array<int, 26>> nxt(n + 1);
    nxt[n].fill(-1);                            // nxt[n][c] = -1：之后没有该字符
    for (int i = n - 1; i >= 0; --i) {
        nxt[i] = nxt[i + 1];
        nxt[i][s[i] - 'a'] = i;
    }
    return nxt;
}
// 判断 t 是否为 s 的子序列：从 pos=0 开始，对每个字符 c：pos = nxt[pos][c] + 1；遇 -1 则否
```

（子序列计数的 DP 用法见 12 章。）

---

## 5.7 经典例题

| 题目 | 考点 | 要点 |
|------|------|------|
| [洛谷 P8306 【模板】字典树](https://www.luogu.com.cn/problem/P8306) | 前缀计数 | 注意输入串数 n/m 与总长；板子直接过 |
| [洛谷 P2580 于是他错误的点名开始了](https://www.luogu.com.cn/problem/P2580) | 点名去重 | `end` 计数 + `visited` 标记重复点名 |
| [洛谷 P4551 最长异或路径](https://www.luogu.com.cn/problem/P4551) | 01-Trie | 树上前缀异或 → 最大异或对 |
| [洛谷 P4735 最大异或和](https://www.luogu.com.cn/problem/P4735) | 可持久化 01-Trie | 前缀异或 + 版本差分区间查询 |
| [洛谷 P3879 [TJOI2010] 阅读理解](https://www.luogu.com.cn/problem/P3879) | Trie + 集合 | 每词存"出现过的文章编号集合" |
| CF 665E *Beautiful Subarrays* | 01-Trie | 子数组异或 ≥ k 计数：前缀异或逐个插查询 |
| [洛谷 P2292 [HNOI2004] L语言](https://www.luogu.com.cn/problem/P2292)（进阶） | Trie + DP | 词典建 Trie，对文本做"从每个可行起点向后匹配单词"的 BFS/DP |

---

## 5.8 小结

```
普通 Trie  ──前缀统计 / 词频 / 字典序 / 补全
    │
    ├─ 01-Trie ──────── 最大异或对、异或路径、子数组异或
    │       └───────── 贪心：高位向反位走
    ├─ 可持久化 Trie ── 区间 [l,r] 内的异或/存在性查询（版本差分）
    │
    └─（接第 06 章）AC 自动机 = Trie + fail 指针 → 多模式匹配
```

> 下一章（06）在 Trie 上添加 **fail 指针**，一次扫描文本同时匹配成千上万个模式串——AC 自动机。
