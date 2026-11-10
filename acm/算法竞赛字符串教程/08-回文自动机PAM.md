# 第 08 章 回文自动机 PAM（Eertree / Palindromic Tree）

> **一句话**：PAM 把字符串中**所有本质不同的回文子串**组织成树，可在 **O(n)** 内求出：
> 本质不同回文子串的个数、**每种回文子串的出现次数**、以及在线（边插入字符边回答）能力。
> 它是回文问题的"后缀自动机"（对回文而言），与 Manacher（07 章）互补：Manacher 管"位置/长度"，PAM 管"种类/次数"。
> 前置：05 章 Trie 的图论直觉。模板与 [OI-wiki 回文树参考实现](https://oi-wiki.org/string/pam/) 一致。

---

## 8.1 关键定理：本质不同回文子串只有 O(n) 个

> **定理**：长度为 n 的串，本质不同的回文子串个数 **≤ n**。

证明要点（归纳）：在 s 后追加字符 c 得到 t，**以 c 结尾**的新回文（即 t 的回文后缀）若多于 1 个，设其中最短的为 P、其余以最长的那个回文后缀覆盖——所有短的回文后缀其实在**旧串 s 中已经出现过**。于是每加一个字符，本质不同回文子串至多增加 1 个。
（详见 OI-wiki 证明：若回文后缀左端点排序为 l1<l2<…，l2 起的回文都被 l1 起的回文"镜像"覆盖，早已出现在 s 中。）

推论：PAM 的状态数（本质不同回文数）+ 2 个根 = **O(n)**，所以可以显式建树。

---

## 8.2 结构

- **每个节点 = 一个本质不同的回文串**（根除外）；
- 转移边：`u --c--> v` 表示在 u 代表的回文串**前后各加一个字符 c**：`str(v) = c + str(u) + c`；
- `fail[u]`：u 的最长**回文后缀**（真后缀）对应的节点；
- 两个根（不表示真实字符串）：
  - 节点 1：**奇根**，`len = -1`（它的"回文"是虚构的，前后加一个字符得到单字符回文）；
  - 节点 0：**偶根**，`len = 0`（真正的空回文串）；
  - `fail[0] = 1`（偶根的最长回文后缀就是空→奇根）；奇根的 fail 无所谓（它不会失配：`-1` 长度的"回文"总能向两边扩展出一个字符）。

这样从偶根扩展得到**偶回文**（len = 2, 4, …），从奇根扩展得到**奇回文**（len = 1, 3, …），结构天然区分奇偶。

---

## 8.3 增量构建

维护 `last` = 当前前缀的**最长回文后缀**节点，字符存入 1 基数组 `s[1..tot]`（`s[0]` 放哨兵）。

插入新字符 `c = s[++tot]`：

1. **找可扩展节点**：从 last 沿 fail 跳，直到节点 u 满足 `s[tot - len[u] - 1] == s[tot]`（u 的回文左边紧邻的字符恰好是 c，于是能扩展成 `c + str(u) + c`）。
   ```cpp
   int getFail(int x) {
       while (s[tot - len[x] - 1] != s[tot]) x = fail[x];
       return x;
   }
   ```
   （当 u 是奇根 len=-1 时，`s[tot+1-1]=s[tot]` 恒成立——这就是"单字符总能扩展"的保证。）

2. **判重 / 建新点**：若 `ch[u][c]` 不存在，新建节点 v（`len[v] = len[u] + 2`），并求 fail：
   ```cpp
   fail[v] = ch[getFail(fail[u])][c];   // 从 u 的 fail 再跳一次得到的“次长回文后缀”扩展
   ch[u][c] = v;
   ```
   若算出来是 0/不存在 → fail[v] = 偶根（空串是最短回文后缀）。
   此时 v 是一个**全新的本质不同回文**（每次插入最多新增 1 个，呼应 8.1）。

3. `last = ch[u][c]; cnt[last]++;`（cnt 先记"作为最长回文后缀出现的次数"）

**复杂度**：每次 getFail 沿 fail 树向上，高度均摊 O(1)；总 O(n)。

---

## 8.4 完整模板（可直接抄，小写字母集）

```cpp
#include <bits/stdc++.h>
using namespace std;

struct PalindromicTree {
    struct Node {
        int ch[26];
        int fail, len;
        long long cnt;          // 出现次数（构建后沿 fail 累加得真实次数）
        Node() { memset(ch, 0, sizeof ch); fail = len = cnt = 0; }
    };
    vector<Node> tr;
    vector<char> s;             // 1 基存储，s[0] = 哨兵
    int tot = 0, last = 0;

    PalindromicTree() { init(); }

    void init() {
        tr.clear();
        tr.resize(2);           // 节点 0 = 偶根(len 0)，节点 1 = 奇根(len -1)
        tr[0].len = 0;          // 偶根
        tr[1].len = -1;         // 奇根
        tr[0].fail = 1;         // 偶根 fail → 奇根
        tr[1].fail = 0;         // 奇根 fail 无所谓（默认 0）
        s.assign(1, '$');       // 哨兵，不在字符集中
        tot = 0;                // 已插入字符数
        last = 0;
    }

    int getFail(int x) {        // 从 x 沿 fail 找第一个能“左右加 c”扩展成回文的节点
        while (s[tot - tr[x].len - 1] != s[tot]) x = tr[x].fail;
        return x;
    }

    // 插入一个字符；返回新形成的（以该位置结尾的）最长回文后缀节点编号
    int insert(char c) {
        s.push_back(c);
        ++tot;
        int u = getFail(last);
        int id = c - 'a';
        if (!tr[u].ch[id]) {
            Node nd;
            tr.push_back(nd);               // 新节点 v
            int v = (int)tr.size() - 1;
            tr[v].len = tr[u].len + 2;      // 两侧各加一个字符
            tr[v].fail = tr[getFail(tr[u].fail)].ch[id];
            tr[u].ch[id] = v;
        }
        last = tr[u].ch[id];
        ++tr[last].cnt;                     // 先记“作为最长回文后缀”的次数
        return last;
    }

    // 建完串后调用：把所有出现次数沿 fail 树累加（节点编号本身是拓扑序，倒序即可）
    void count() {
        for (int i = (int)tr.size() - 1; i >= 2; --i)   // 跳过两个根
            tr[tr[i].fail].cnt += tr[i].cnt;
    }
};
```

**主程序用法示例（[P3649 [APIO2014] 回文串](https://www.luogu.com.cn/problem/P3649)）**：

```cpp
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string s;
    cin >> s;
    PalindromicTree pam;
    for (char c : s) pam.insert(c);
    pam.count();
    long long ans = 0;
    for (int i = 2; i < (int)pam.tr.size(); ++i)        // 跳过 0(偶根) 1(奇根)
        ans = max(ans, pam.tr[i].len * pam.tr[i].cnt);  // 长度 × 出现次数
    cout << ans << '\n';
    return 0;
}
```

> 注意下标约定：**节点 0 是偶根、节点 1 是奇根**（与部分博客相反），抄代码前先看清初始化，否则 `len=-1` 的节点会让你全盘出错。真节点从 2 开始，故"本质不同回文数 = 节点总数 − 2"。

---

## 8.5 核心应用

### 8.5.1 本质不同回文子串个数

```cpp
// 插完整个串后：
int distinctPalCount = (int)pam.tr.size() - 2;
```

### 8.5.2 每种回文的出现次数（核心卖点）

插入过程中每个位置产生的最长回文后缀对应节点被 `cnt++`；某个回文 P 每出现一次，要么它是某个位置的最长回文后缀，要么它是该最长回文后缀沿 fail 链的祖先。因此**把 cnt 沿 fail 树自底向上累加**，`cnt[P 的节点]` 就是 P 的出现次数（模板 `count()` 已实现；节点插入顺序天然满足"子先于父"，倒序扫描即可，无需额外拓扑排序）。

### 8.5.3 以每个位置结尾的回文信息（在线）

`insert` 返回 `last`：若本轮新建了节点，则出现了一个新的本质不同回文；`len[last]` 是以当前位置结尾的最长回文长度。可以边读串边回答，如 [P5496 【模板】回文自动机](https://www.luogu.com.cn/problem/P5496)。

### 8.5.4 最小回文划分（进阶，O(n log n)）

**问题**：把 s 划分成若干回文串，求最少段数 / 方案数。
朴素：`dp[i] = 1 + min dp[j]`（s[j+1..i] 为回文），枚举回文后缀 O(n²)。

**优化（series link / slink）**：回文串后缀有很强的结构（OI-wiki 引理）：记 `diff[u] = len[u] - len[fail[u]]`，`slink[u]` = 沿 fail 一直跳到第一个 `diff` 不同的节点。可以证明：**某个串的所有回文后缀按长度排序后只形成 O(log n) 段等差数列**。于是每个等差段整体转移：

```
g[u] = 该等差数列对 dp 的贡献（维护成段最值/和）
更新：g[u] = g[fail[u]] + (或 min) dp[i - (len[slink[u]] + diff[u])]  （位置换算见 OI-wiki）
dp[i] 用沿 slink 链跳 O(log n) 次得到
```

参考实现与完整引理：OI-wiki《回文树》"最小回文划分"一节；经典题 [CF 932G Palindrome Partition](https://codeforces.com/problemset/problem/932/G)（划分成前后对称的偶数段方案数，先把串交错重排再数偶回文划分）。

### 8.5.5 双端插入等扩展

PAM 还能支持**在串首/串尾双向插入**（维护 lastL、lastR 与最左边界），见 HDU 5421 Victor and String（进阶）。

---

## 8.6 经典例题

| 题目 | 考点 | 要点 |
|------|------|------|
| [洛谷 P5496 【模板】回文自动机（PAM）](https://www.luogu.com.cn/problem/P5496) | 在线插入 | 每次输出以当前位置结尾的新增/回文信息，板子直出 |
| [洛谷 P3649 [APIO2014] 回文串](https://www.luogu.com.cn/problem/P3649) | 次数 × 长度 | `count()` 后扫真节点取 max |
| [洛谷 P1659 [国家集训队] 拉拉队排练](https://www.luogu.com.cn/problem/P1659) | 按长度统计 | 也可用 PAM：每节点有 len 与 cnt，按 len 分组排序 |
| [洛谷 P4287 [SHOI2011] 双倍回文](https://www.luogu.com.cn/problem/P4287) | 回文 + 结构 | 需要"半长仍为回文"的判定，Manacher/PAM 均可做，注意 O(n) 剪枝 |
| [CF 932G Palindrome Partition](https://codeforces.com/problemset/problem/932/G) | 最小回文划分（slink） | 进阶 O(n log n) |
| CodeChef PALPROB *Palindromeness* | 半回文权重 | 利用 fail 树与 len 的递推 |

---

## 8.7 PAM vs Manacher vs 哈希

| 需求 | 推荐工具 |
|------|---------|
| 最长回文子串 / 回文半径 | Manacher O(n) |
| 回文子串按"位置"计数 | Manacher（差分） |
| **本质不同**回文种类 | **PAM** |
| 每种回文的**出现次数** | **PAM** |
| 在线（边读边答） | PAM（哈希也行但碰撞风险） |
| 最小回文划分 | PAM + slink |

> 下一章（09）进入"后缀三剑客"第一剑：**后缀数组 SA**——把 n 个后缀排好序，配上 height，能解决几乎一切"子串即后缀前缀"的计数/查询问题。
