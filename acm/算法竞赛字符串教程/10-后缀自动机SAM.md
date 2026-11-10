# 第 10 章 后缀自动机（Suffix Automaton, SAM）

> **一句话**：SAM 是能接受一个串**所有后缀**的**规模最小**的自动机（状态 ≤ 2n−1、转移 ≤ 3n−4），
> 它能**在线**（每加一个字符 O(1) 均摊）地表示全部子串信息。
> 学完它，你将有：本质不同子串数、任意子串出现次数、字典序第 k 小子串、两串/多串最长公共子串……全部 O(n)。
> 本章构造算法与 [CP-Algorithms](https://cp-algorithms.com/string/suffix-automaton.html) / [OI-wiki](https://oi-wiki.org/string/sam/) 一致。

---

## 10.1 两个直觉

1. **子串 = 后缀的前缀**。若自动机接受所有后缀，则从起点出发的任意路径都对应某个后缀的前缀 = 一个子串；反之每个子串都对应一条路径。因此 **SAM 的所有路径 ↔ s 的所有子串**（一一对应到"路径的字符序列"）。
2. SAM 是"全体后缀的 AC 自动机的**压缩**"：把 endpos 相同的子串合并进同一状态，规模从 O(n²) 降到 O(n)。

| 能做什么 | 复杂度 |
|---|---|
| 本质不同子串个数 / 所有不同子串总长 | O(n) |
| 每个子串的出现次数 | O(n) 预处理 |
| 字典序第 k 小子串 | O(n) 预处理 + O(答案·Σ) 查询 |
| 两串最长公共子串 | O(n+m) |
| 最小循环移位 / 判断子串 | O(n) |
| 在线统计"每加一个字符新增多少子串" | 均摊 O(1)/字符 |

---

## 10.2 核心概念

### 10.2.1 endpos 等价类（SAM 状态 = 等价类）

对子串 t，记 `endpos(t)` = t 在 s 中**所有出现位置的结束下标集合**。例 `s = "abcbc"`：`endpos("bc") = {2,4}`。

- **引理**：两个子串的 endpos 相同 ⇔ 较短者是较长者的后缀、且只在"作为它的后缀"时出现。
- 因此每个状态 v 对应**一族互为后缀的子串**，长度连续取遍 `[minlen(v), len(v)]`；
- **len(v)** = 该状态最长的子串长度；**link(v)**（后缀链接）连到"下一个更短的、endpos 集合不同的最长后缀"所属状态，且满足 **minlen(v) = len(link(v)) + 1**。

### 10.2.2 后缀链接树（parent 树 / link 树）

link 构成一棵以根（状态 0，空串，len=0，link=-1）为根的树，并且：
- `endpos(link(v)) ⊋ endpos(v)`，即 link 树上的祖先的 endpos 是后代 endpos 的超集；
- **子树内节点数/信息聚合** 是统计一切"出现次数"类问题的根基。

### 10.2.3 状态数与转移数

定理：s 长 n，SAM 状态数 **≤ 2n−1**，转移数 **≤ 3n−4**。这是"为什么敢开 2n 数组"的依据。

---

## 10.3 构造算法（增量，逐字符扩展）

维护 `last` = 代表**整个当前串**的状态。加入新字符 c：

1. 新建状态 `cur`，`len[cur] = len[last] + 1`；
2. 从 last 沿 link 上溯，把所有"没有字符 c 转移"的状态补上指向 cur 的转移；
3. 若上溯到 -1（根之上）：`link[cur] = 0`；
   否则停在第一个有 c 转移的状态 p，设 q = 转移目标：
   - 若 `len[p] + 1 == len[q]`：`link[cur] = q`（q 正好是"最长后缀 + c"）；
   - 否则 **分裂**：复制 q 得到 clone（`len[clone] = len[p] + 1`，转移与 link 都复制 q 的），把从 p 开始所有指向 q 的 c 转移改指向 clone，最后 `link[q] = link[cur] = clone`；
4. `last = cur`。

```cpp
#include <bits/stdc++.h>
using namespace std;

struct SuffixAutomaton {
    struct State {
        int next[26];      // 字符集大时换 map/unordered_map（见 10.3.1）
        int link, len;
        State() { memset(next, 0, sizeof next); link = -1; len = 0; }
    };
    vector<State> st;
    int last;

    SuffixAutomaton() { init(); }
    void init() {
        st.clear();
        st.emplace_back();        // 状态 0 = 根（空串），len=0, link=-1
        last = 0;
    }

    void extend(int c) {          // c: 0..25
        int cur = (int)st.size();
        st.emplace_back();
        st[cur].len = st[last].len + 1;

        int p = last;
        while (p != -1 && !st[p].next[c]) {   // 给所有缺 c 转移的后缀状态补转移
            st[p].next[c] = cur;
            p = st[p].link;
        }
        if (p == -1) {
            st[cur].link = 0;                 // 情况一：c 从未出现过
        } else {
            int q = st[p].next[c];
            if (st[p].len + 1 == st[q].len) {
                st[cur].link = q;             // 情况二
            } else {                          // 情况三：分裂 q -> clone
                int clone = (int)st.size();
                st.emplace_back();
                st[clone].len = st[p].len + 1;
                memcpy(st[clone].next, st[q].next, sizeof st[q].next);
                st[clone].link = st[q].link;
                while (p != -1 && st[p].next[c] == q) {
                    st[p].next[c] = clone;    // 重定向
                    p = st[p].link;
                }
                st[q].link = st[cur].link = clone;
            }
        }
        last = cur;
    }
};
```

**复杂度**：均摊 O(1)/字符（严格证明较绕，记住结论即可）；整串构造 O(n)。注意 clone 的 len 可能小于之后创建的状态，所以**创建顺序不是按 len 有序**——需要拓扑序时要按 len 基数排序（见 10.4）。

### 10.3.1 字符集与空间取舍

- 小写字母 + n ≤ 2e5：`int next[26]` 数组最快（约 26×2n×4 字节 ≈ 40MB for n=2e5? 实际 2n×26×4 ≈ 41.6MB for n=2e5, ok; n=1e6 会超 200MB → 换方案）；
- 字符集大或 n 大：每个状态用 `map<int,int>`（O(log) 转移）或把转移存成链表 + 二分（离散化字符）；
- 经典 trick：转移数组先 `vector<array<int,26>>` + `reserve(2n)` 减少内存碎片。

---

## 10.4 拓扑序与"出现次数"（endpos 大小）

**统计口径**：每个"非 clone 状态"对应一个完整前缀的结尾（创建时 endpos 恰好含一个新位置），所以给所有非 clone、非根状态初值 `cnt = 1`；然后把 cnt 沿 link 树**自底向上**累加，`cnt[v] = |endpos(v)|` = 状态 v 里每个子串的出现次数。

自底向上需要按 **len 降序**，用计数排序 O(n)：

```cpp
// 在建好的 SAM 上调 used count 数组：先统计出现次数
vector<int> cnt(st.size(), 0);
for (int v = 1; v < (int)st.size(); ++v)
    if (st[v].link != -1 && /* v 非 clone */ true) cnt[v] = 1;
// 注：clone 状态需要额外标记 is_clone；简单写法是只给"extend 中新建且未被分裂复制"的状态初值 1，
// 下面用最常见的做法：记录每个状态是否 clone。
```

稳妥写法：给 State 加 `bool is_clone`；`extend` 里新建 cur 时 `is_clone=false`，clone 时 `is_clone=true`。然后：

```cpp
int n = (int)st.size();
vector<int> order(n), bucket(n + 1, 0);
for (auto& s : st) ++bucket[s.len];
for (int i = 1; i <= n; ++i) bucket[i] += bucket[i - 1];
for (int i = n - 1; i >= 0; --i) order[--bucket[st[i].len]] = i;   // 按 len 升序

vector<int> cnt(n, 0);
for (int v = 1; v < n; ++v) if (!st[v].is_clone) cnt[v] = 1;       // 每个位置贡献 1
for (int i = n - 1; i >= 0; --i) {          // len 降序 = 拓扑逆序
    int v = order[i];
    if (st[v].link != -1) cnt[st[v].link] += cnt[v];
}
// cnt[v] = 状态 v 中任意子串（如 longest(v)）在 s 中出现的次数
```

> [P3804 【模板】后缀自动机](https://www.luogu.com.cn/problem/P3804)：求 `max{ len[v] × cnt[v] | cnt[v] > 1 }`。

**为什么 clone 初值为 0**：clone 不"拥有"新的结束位置，它只是把 q 的一部分位置拆出来；那些位置已经由它复制到的祖先/子孙统计了。

---

## 10.5 应用大礼包

### ① 本质不同子串个数 / 总长度

状态 v 新增的子串数为 `len[v] − len[link[v]]`（长度区间 `[len(link)+1, len]` 里的子串首次在此状态出现）：

```cpp
long long distinctCount() {
    long long ans = 0;
    for (int v = 1; v < (int)st.size(); ++v)
        ans += st[v].len - st[st[v].link].len;   // link[v] 的 len
    return ans;
}
// 总长度 = Σ len[v]*(len[v]+1)/2 - len[link[v]]*(len[link[v]]+1)/2（等差数列，见 OI-wiki）
```

**在线版**（[P4070 [SDOI2016] 生成魔咒](https://www.luogu.com.cn/problem/P4070)，字符集大、next 需用 map）：每 extend 一个字符，**新增的本质不同子串数 = `len[cur] − len[link[cur]]`**。原因：新的本质不同子串只能是"以新位置结尾且此前未出现过"的后缀，长度恰好取 `(len[link[cur]], len[cur]]`。逐字符累加即为全程答案（分裂产生的 clone 不增加总量：它代表的子串此前已被其它状态统计过）。

### ② 字典序第 k 小子串（[P3975 [TJOI2015] 弦论](https://www.luogu.com.cn/problem/P3975)）

两种口径：
- **t = 0**：相同的子串算一种 → 每个状态子串数 = `len[v] - len[link[v]]`；
- **t = 1**：按出现位置计数（相同内容多次出现算多次）→ 每个状态子串权重 = `cnt[v]`（endpos 大小）。

做法：先按 len 降序（拓扑逆序）计算 `dp[v]` = 从 v 出发（含空串）的路径条数（按口径给"停在 v"加权），然后从根出发贪心：从小到大试每个字符 ch，若 `dp[目标] ≥ k` 就走并输出 ch，否则 `k -= dp[目标]`。

```cpp
// 伪代码（t=1 口径）：ways[v] = 从 v 开始还能形成的子串数（含空路径）
// ways[v] = cnt[v](停在v, t=1)/1(t=0) + Σ ways[next[v][c]]
// 查询：cur=0; while(k>0){ for c: 若 cur 有转移 c 且 ways[next]>=k → 输出c,cur=next,k-=cnt? 见题解细粒度 }
```

> 此题经典提示：先在 SAM 上跑 ① 的 DP（每状态贡献 cnt[v]），再按字典序贪心。O(n·Σ)。

### ③ 两串最长公共子串（SPOJ LCS）

对 S 建 SAM，用 T 在上面"尽量匹配"：维护 `(v, l)` = 当前匹配到的状态与长度。加字符 c：
- 有转移 → `v = next[v][c], l++`；
- 没转移 → 沿 link 回退缩短（`v = link[v], l = len[v]`），直到有转移或到根；
答案 = 过程中 l 的最大值。总 O(|S| + |T|)。

```cpp
int lcsLength(const string& S, const string& T) {
    SuffixAutomaton sam;
    for (char ch : S) sam.extend(ch - 'a');
    int v = 0, l = 0, best = 0;
    for (char ch : T) {
        int c = ch - 'a';
        while (v != -1 && !sam.st[v].next[c]) {     // 注意根(0)无转移时也退
            v = sam.st[v].link;
            if (v != -1) l = sam.st[v].len;
            else { v = 0; l = 0; break; }
        }
        if (sam.st[v].next[c]) { v = sam.st[v].next[c]; ++l; }
        best = max(best, l);
    }
    return best;
}
```

### ④ 多串最长公共子串（SPOJ LCS2）

对**最短**的串建 SAM；对其它每个串跑 ③ 的匹配，同时记录"每个状态被当前串匹配到的最长长度 mx[v]"；跑完一个串后按 len 降序**把子节点的 mx 上传到 link 祖先**（`mx[link[v]] = max(mx[link[v]], min(mx[v], len[link[v]]))`）。全部串处理完后，答案 = `max_v min_v(mx over all strings)`。复杂度 O(Σ|Sᵢ|)。（详见 OI-wiki LCS2 解法二。）

### ⑤ 最小循环移位 / 最小表示

对 `S+S` 建 SAM，从根出发贪心走**最小字符**转移共 n 步，得到的就是字典序最小的长度 n 的循环移位。（也可用第 11 章的 Lyndon/最小表示法，更省空间。）

### ⑥ 子串出现位置 / 判断子串

- **判断 P 是否是 s 的子串**：沿转移走，走不完则不是；O(|P|)。
- **找 P 出现次数**：走完得到状态 v，答案 = cnt[v]（10.4）。
- **第一次出现位置**：extend 时给 cur 记 `firstpos = len-1`，clone 继承 q 的 firstpos；答案 = firstpos[v] − |P| + 1。
- **所有出现位置**：对 link 树建反图（children 列表），从 P 所在状态 v 沿反图 DFS/BFS，输出所有非 clone 节点的 firstpos。

### ⑦ 广义 SAM（多串共用）

把多串依次插入，**每串开头把 last 重置为 0**。插入前做一次"如果 last 已有字符 c 的转移则直接走过去"的判重（经典写法），避免重复计数。例题：P6139【模板】广义后缀自动机（广义 SAM）。进阶可研究"SAM on Trie"。

### ⑧ SAM + 其它数据结构的进阶方向（了解即可）

- SAM 的 link 树 + 线段树合并：求"某个子串在区间 [l,r] 内是否出现 / 出现次数"；
- SAM + DP / 二分（[P4022 [CTSC2012] 熟悉的文章](https://www.luogu.com.cn/problem/P4022)）；
- [NOI2018 你的名字](https://uoj.ac/problem/395) 级别的综合（SAM 套 SAM / 可持久化线段树）。

---

## 10.6 经典例题

| 题目 | 考点 | 要点 |
|------|------|------|
| [洛谷 P3804 【模板】后缀自动机](https://www.luogu.com.cn/problem/P3804) | 出现次数 | cnt 拓扑累加，max(len×cnt) |
| [洛谷 P3975 [TJOI2015] 弦论](https://www.luogu.com.cn/problem/P3975) | 字典序第 k 小 | 两种口径 DP + 贪心走边 |
| [洛谷 P4070 [SDOI2016] 生成魔咒](https://www.luogu.com.cn/problem/P4070) | 在线本质不同 | 每步新增 = Σ(新建状态 len − len[link]) |
| SPOJ LCS / LCS2 | 最长公共子串（2/多串） | ③④ 的方法 |
| [洛谷 P6139 【模板】广义后缀自动机](https://www.luogu.com.cn/problem/P6139) | 广义 SAM | 多串共用 + 判重 |
| CF 235C *Cyclical Quest* | 循环同构查询 | SAM 上匹配所有旋转并去重 |
| CF 802I / CF 873F | 出现次数 × 长度 | cnt + DP |
| [洛谷 P3346 诸神眷顾的幻想乡](https://www.luogu.com.cn/problem/P3346)（进阶） | Trie 上 SAM | 广义 SAM 应用题 |

---

## 10.7 SAM vs SA（选谁？）

| | 后缀数组 SA | 后缀自动机 SAM |
|---|---|---|
| 构建 | O(n log n)，常数小 | O(n)，常数略大、内存略大 |
| 在线加字符 | 不行（重排） | **可以** |
| 任意子串出现次数 | 二分 + 区间 | cnt[v] 直出 |
| 本质不同子串数 | Σ(n−sa+1)−Σheight | Σ(len−len[link]) |
| 字典序第 k 小 | 很好写 | 可以（DP） |
| 多模式串 | 拼接分隔符 | 广义 SAM |
| 直观程度 | 较直观 | 概念门槛高 |

> **建议**：两道都背熟。SA 代码简单粗暴适合"查区间/字典序"；SAM 适合"出现次数/增量/在线"。字符串大题往往是 SA/SAM 与数据结构、DP 的联合（如 NOI2015 品酒大会、NOI2018 你的名字）。

> 下一章（11）回到相对轻量的主题：**Lyndon 分解与最小表示法**——循环同构问题的最优解法。
