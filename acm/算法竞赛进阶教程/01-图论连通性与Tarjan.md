# 第 01 章 图论连通性与 Tarjan 全家桶

> 你已经会 Tarjan 的**离线 LCA**（`lca_Tarjan.cpp`），但那只是 Tarjan 算法的一个"副产品用法"。
> 同一个 `dfn / low` 框架还能求出：**强连通分量、割点、桥、点双、边双、2-SAT 的解**。
> 这一章把它们一次讲完——**这是图论进阶的分水岭**，后面网络流建模、圆方树、虚树都靠它打底。

---

## 1.1 统一框架：dfn、low 与 DFS 树

**准备工作**：对图做一次 DFS，记录每个点的**时间戳** `dfn[u]`（第几个被访问到），并定义

$$\text{low}[u] = \min\Bigl(\ \text{dfn}[u],\ \min_{u \to v \text{ 是树边}} \text{low}[v],\ \min_{u \to v \text{ 是非树边}} \text{dfn}[v]\ \Bigr).$$

> 直觉：`low[u]` 表示"从 $u$ 的子树出发，**通过至多一条非树边**能回到的最小 dfn"。
> 也就是"$u$ 的子树有没有'后路'通向上方"。

**DFS 树的四种边**（有向图）：

| 类型 | 特征 | 作用 |
|---|---|---|
| 树边 | `!dfn[v]` | 骨架 |
| 返祖边（back edge） | `dfn[v] < dfn[u]` 且 $v$ 是 $u$ 的祖先 | **构成环**，SCC 的关键 |
| 前向边（forward edge） | `dfn[v] > dfn[u]` 且 $v$ 在 $u$ 子树中 | 不影响 low |
| 横叉边（cross edge） | 其他已访问点 | **只有指向栈中点的才影响 low** |

**这是所有 Tarjan 算法唯一的"心法"**：

- 对**已访问且在栈中**的点 $v$，用 `dfn[v]` 更新 `low[u]`；
- 对已访问但**已出栈**的点 $v$，**忽略**（它属于别的已定型的连通块）。

---

## 1.2 强连通分量 SCC

### 定义

有向图中，若 $u, v$ 可互相到达，则称它们强连通。极大强连通子图称为**强连通分量（SCC）**。
把所有 SCC 缩成一个点后，得到的图一定是 **DAG**（否则那两个 SCC 就该合并）。

### Tarjan 求 SCC

维护一个栈：DFS 过程中把访问到的点压栈；当 `low[u] == dfn[u]` 时，说明 $u$ 是它所在 SCC 的"根"，此时把栈中 $u$ 及以上的点全部弹出，它们构成一个 SCC。

```cpp
// ==================== Tarjan 求强连通分量 ====================
const int MAXN = 100005;
vector<int> g[MAXN];          // 有向图邻接表
int dfn[MAXN], low[MAXN], timer;
int stk[MAXN], top;
bool inStk[MAXN];
int sccId[MAXN], sccCnt;      // 每个点所属 SCC 编号（1..sccCnt）
int sccSize[MAXN];

void tarjan_scc(int u) {
    dfn[u] = low[u] = ++timer;
    stk[++top] = u; inStk[u] = true;
    for (int v : g[u]) {
        if (!dfn[v]) {                       // 树边
            tarjan_scc(v);
            low[u] = min(low[u], low[v]);
        } else if (inStk[v]) {               // 指向栈中点的边（返祖边/横叉边）
            low[u] = min(low[u], dfn[v]);
        }
        // 指向已出栈点的边：忽略
    }
    if (low[u] == dfn[u]) {                  // u 是 SCC 的根
        ++sccCnt;
        int v;
        do {
            v = stk[top--];
            inStk[v] = false;
            sccId[v] = sccCnt;
            ++sccSize[sccCnt];
        } while (v != u);
    }
}

// 调用：for (int i = 1; i <= n; i++) if (!dfn[i]) tarjan_scc(i);
```

**复杂度**：$O(n + m)$，空间 $O(n + m)$。

**递归深度问题**：链状图（如 $10^6$ 个点排成一条链）会爆栈。三种处理方式：

1. 把 `tarjan_scc` 改写成显式栈的迭代版本（推荐，见 §1.2.3）；
2. 开头加 `ios::sync_with_stdio(false)` 之外的编译指令扩大栈：Windows 下 `#pragma comment(linker, "/STACK:...")` 不可靠；
3. 手写栈模拟（最稳妥）。

**迭代版 SCC（避免爆栈）**：

```cpp
// 迭代版 Tarjan（用数组模拟递归栈），适合链状大图
void tarjan_scc_iterative(int start) {
    static int callStack[MAXN * 2], edgeIter[MAXN];
    int csTop = 0;
    callStack[++csTop] = start;
    dfn[start] = low[start] = ++timer;
    stk[++top] = start; inStk[start] = true;
    edgeIter[start] = 0;
    while (csTop) {
        int u = callStack[csTop];
        if (edgeIter[u] < (int)g[u].size()) {
            int v = g[u][edgeIter[u]++];
            if (!dfn[v]) {
                dfn[v] = low[v] = ++timer;
                stk[++top] = v; inStk[v] = true;
                edgeIter[v] = 0;
                callStack[++csTop] = v;
            } else if (inStk[v]) {
                low[u] = min(low[u], dfn[v]);
            }
        } else {
            // 所有出边处理完毕
            if (low[u] == dfn[u]) {
                ++sccCnt;
                int v;
                do {
                    v = stk[top--];
                    inStk[v] = false;
                    sccId[v] = sccCnt;
                    ++sccSize[sccCnt];
                } while (v != u);
            }
            --csTop;
            if (csTop) {
                int p = callStack[csTop];
                low[p] = min(low[p], low[u]);      // 回溯时更新父节点 low
            }
        }
    }
}
```

> **注意**：迭代版最容易错的地方是"回溯时更新父节点的 low"，必须在弹出子节点后立即做。

### 缩点与 DAG 上的 DP

```cpp
// ==================== 缩点：把 SCC 缩成 DAG ====================
vector<int> dag[MAXN];        // 缩点后的图
int val[MAXN], dagVal[MAXN];  // 点权

void build_dag(int n) {
    for (int u = 1; u <= n; u++) {
        dagVal[sccId[u]] += val[u];                 // 合并点权
        for (int v : g[u])
            if (sccId[u] != sccId[v]) dag[sccId[u]].push_back(sccId[v]);
    }
    for (int i = 1; i <= sccCnt; i++) {             // 去重（可选，但能显著减小常数）
        sort(dag[i].begin(), dag[i].end());
        dag[i].erase(unique(dag[i].begin(), dag[i].end()), dag[i].end());
    }
}

// DAG 上最长路（按 sccId 从大到小就是拓扑序！）
// Tarjan 求出的 SCC 编号天然满足：若 sccId[u] < sccId[v] 则不存在 v -> u 的路径
// 即 sccId 从大到小 = 拓扑序
int longest_path_dag() {
    vector<int> dp(sccCnt + 1, 0);
    int ans = 0;
    for (int u = sccCnt; u >= 1; u--) {             // 逆序 = 拓扑序
        dp[u] += dagVal[u];
        ans = max(ans, dp[u]);
        for (int v : dag[u]) dp[v] = max(dp[v], dp[u]);
    }
    return ans;
}
```

> **重要结论**：Tarjan 的 `sccId` **天然是逆拓扑序**（编号大的在拓扑序前面）。
> 所以缩点后做 DAG 上的 DP **不需要再做一次拓扑排序**——直接倒序/正序遍历 `sccId` 即可。
> 忘了这一点不会错，但会白白多写一段拓扑排序代码。

**经典应用**：

| 题目 | 做法 |
|---|---|
| [P2341 受欢迎的牛](https://www.luogu.com.cn/problem/P2341) | 缩点后统计出度为 0 的 SCC 个数，恰为 1 时答案是它的大小 |
| [P3387 【模板】缩点](https://www.luogu.com.cn/problem/P3387) | 缩点 + DAG 最长路 |
| [P2812 校园网络](https://www.luogu.com.cn/problem/P2812) | 缩点后：入度 0 的个数 / max(入度0, 出度0) |
| [P2272 最大半连通子图](https://www.luogu.com.cn/problem/P2272) | 缩点 + DAG 上最长链计数 |

---

## 1.3 割点与桥（无向图）

### 定义与判定

**割点**：删掉后使图的连通块数增加的点。
**桥（割边）**：删掉后使图的连通块数增加的边。

设 DFS 树上 $v$ 是 $u$ 的儿子：

$$\text{桥}(u,v) \iff \text{low}[v] > \text{dfn}[u],$$

$$\text{割点}(u) \iff \exists\, v \text{ 是 } u \text{ 的儿子},\ \text{low}[v] \ge \text{dfn}[u] \quad (\text{且 } u \text{ 不是 DFS 根}),$$

$$\text{根 } u \text{ 是割点} \iff u \text{ 在 DFS 树上有} \ge 2 \text{ 个儿子}.$$

**直观解释**：`low[v] >= dfn[u]` 意味着 $v$ 的子树**没有后路**绕过 $u$ 回到 $u$ 的祖先——所以删掉 $u$ 后 $v$ 的子树就"孤立"了。桥是更严格的情形（连回到 $u$ 自己都不行）。

### 代码（正确处理重边）

**关键技巧**：用**边的编号**而不是"父节点编号"来跳过回路，这样才能处理重边（两条 $u$-$v$ 边中的一条是树边，另一条是合法的返祖边）。

```cpp
// ==================== 无向图：割点与桥 ====================
// 边从编号 2 开始成对存储：i 与 i^1 互为反向边
const int MAXN = 100005, MAXM = 500005;
struct Edge { int to, nxt; } e[MAXM * 2];
int head[MAXN], ecnt = 1;

void addEdge(int u, int v) {
    e[++ecnt] = {v, head[u]}; head[u] = ecnt;
    e[++ecnt] = {u, head[v]}; head[v] = ecnt;
}

int dfn[MAXN], low[MAXN], timer;
bool isCut[MAXN];              // 割点
bool isBridge[MAXM * 2];       // 桥（按边编号标记，成对标记）

void tarjan_cut(int u, int inEdge) {   // inEdge: 进入 u 的边的编号，根为 0
    dfn[u] = low[u] = ++timer;
    int child = 0;
    for (int i = head[u]; i; i = e[i].nxt) {
        int v = e[i].to;
        if (i == (inEdge ^ 1)) continue;        // 不走来时的那条边（重边安全）
        if (!dfn[v]) {
            ++child;
            tarjan_cut(v, i);
            low[u] = min(low[u], low[v]);
            if (low[v] > dfn[u]) isBridge[i] = isBridge[i ^ 1] = true;
            if (low[v] >= dfn[u] && inEdge != 0) isCut[u] = true;
        } else {
            low[u] = min(low[u], dfn[v]);       // 无向图：指向祖先（含重边的另一端）
        }
    }
    if (inEdge == 0 && child > 1) isCut[u] = true;   // 根节点特判
}
```

> **为什么无向图里 `else` 分支不需要判"在栈中"？** 无向图的 DFS 没有横叉边：所有非树边都是返祖边（要么指向祖先，要么指向子树中的点——后者在递归返回时已被处理为 `dfn[v] > dfn[u]`，用 `min` 更新不会产生错误影响）。
> 严格的写法是 `else if (dfn[v] < dfn[u])`，加上它更安全。

**经典应用**：

| 题目 | 做法 |
|---|---|
| [P3388 【模板】割点](https://www.luogu.com.cn/problem/P3388) | 直接输出割点 |
| [P3469 [POI2008] BLO-Blockade](https://www.luogu.com.cn/problem/P3469) | 统计删点后不连通点对数：Tarjan 过程中累加 `size` 贡献 |
| [P5058 [ZJOI2004] 嗅探器](https://www.luogu.com.cn/problem/P5058) | 求两点间**必经点**（Tarjan 过程中判断） |

**BLO-Blockade 的核心推导**（很典型，背下来）：

删掉 $u$ 后，图被分成若干块：每个满足 `low[v] >= dfn[u]` 的儿子子树 $v$（大小 $s_v$）、以及"其余部分"（大小 $n - 1 - \sum s_v$）、$u$ 自己。
不连通点对数 = $2\sum_{\text{块}} \text{(块大小)} \times (n - 1 - \text{块大小})$ 的一半…… 更直接地：

$$\text{ans}[u] = 2\Bigl[(n-1) + \sum_{v} s_v\,(n - 1 - s_v)\Bigr]$$

其中求和只对"被 $u$ 分出去"的子树做（$u$ 不是割点时答案是 $2(n-1)$）。

---

## 1.4 边双连通分量（e-BCC）

**定义**：一个**极大子图**，其中任意两点间存在两条**边不相交**的路径（等价于：子图内没有桥）。

**求法**：先求所有桥，然后"不走桥"地 DFS 染色，每个连通块就是一个 e-BCC。

```cpp
// ==================== 边双连通分量 ====================
// 依赖 §1.3 的无向图存储与 isBridge 数组（下面 7 行同前，合并到同一文件时删掉）
const int MAXN = 100005, MAXM = 500005;
struct Edge { int to, nxt; } e[MAXM * 2];
int head[MAXN], ecnt = 1;
void addEdge(int u, int v) {
    e[++ecnt] = {v, head[u]}; head[u] = ecnt;
    e[++ecnt] = {u, head[v]}; head[v] = ecnt;
}
bool isBridge[MAXM * 2];

int bel[MAXN], bcnt;              // 每个点所属 e-BCC 编号

void dfs_ebcc(int u) {
    bel[u] = bcnt;
    for (int i = head[u]; i; i = e[i].nxt) {
        int v = e[i].to;
        if (isBridge[i] || bel[v]) continue;
        dfs_ebcc(v);
    }
}

void build_ebcc(int n) {
    bcnt = 0;
    for (int i = 1; i <= n; i++)
        if (!bel[i]) { ++bcnt; dfs_ebcc(i); }
}

// 把 e-BCC 缩成一棵树：桥就是树边
vector<int> tree_e[MAXN];
void ebcc_to_tree(int n) {
    for (int u = 1; u <= n; u++)
        for (int i = head[u]; i; i = e[i].nxt) {
            int v = e[i].to;
            if (isBridge[i] && bel[u] != bel[v]) {
                tree_e[bel[u]].push_back(bel[v]);
                // 只加一次（成对边会加两次，最后去重即可）
            }
        }
    for (int i = 1; i <= bcnt; i++) {
        sort(tree_e[i].begin(), tree_e[i].end());
        tree_e[i].erase(unique(tree_e[i].begin(), tree_e[i].end()), tree_e[i].end());
    }
}
```

**性质**：e-BCC 缩点后一定是**树**（因为桥就是树边，且任意两个 e-BCC 之间至多一条桥——否则会有环，就不需要桥了）。

**经典应用**：树上问题全部可用（直径、LCA、树上差分）。

| 题目 | 做法 |
|---|---|
| [P2860 冗余路径](https://www.luogu.com.cn/problem/P2860) | 缩点成树，答案 = $\lceil \frac{\text{叶子数}}{2}\rceil$ |
| [P8436 【模板】边双连通分量](https://www.luogu.com.cn/problem/P8436) | 模板题 |
| [P3225 矿场搭建](https://www.luogu.com.cn/problem/P3225) | 点双 + 分类讨论（见下节） |

---

## 1.5 点双连通分量（v-BCC）与圆方树

### 定义与性质

**点双连通**：任意两点间存在两条**点不相交**的路径（等价于：删掉任意一个点后仍连通）。

**关键区别**：

| | 边双 | 点双 |
|---|---|---|
| 定义 | 无桥 | 无割点 |
| 交点情况 | 两两不相交（点可以重合？不，边双是划分） | **一个割点可能属于多个点双** |
| 缩点后 | 一棵树 | 圆方树 |

### 求点双：栈中存"边"

求点双的经典做法是**把边压栈**：当你发现 `low[v] >= dfn[u]`（即 $u$ 是割点或根），就把栈中从 $(u,v)$ 开始的边全部弹出，这些边构成一个点双。

```cpp
// ==================== 点双连通分量 + 圆方树 ====================
// 圆方树：原图点称为"圆点"(1..n)，每个点双新建一个"方点"(n+1 ... n+tot)
// 下面 6 行与 §1.3 的无向图存储相同；放在同一个文件里时只需保留一份
const int MAXN = 100005, MAXM = 500005;
struct Edge { int to, nxt; } e[MAXM * 2];
int head[MAXN], ecnt = 1;
void addEdge(int u, int v) {
    e[++ecnt] = {v, head[u]}; head[u] = ecnt;
    e[++ecnt] = {u, head[v]}; head[v] = ecnt;
}

int vbDfn[MAXN], vbLow[MAXN], vbTimer;
int edgeStk[MAXM * 2], edgeTop;

int n0;                        // 原图点数
int total;                     // 圆方树总点数（初始为 n0）
vector<int> cft[MAXN * 2];     // 圆方树（circle-square tree）

void tarjan_vbcc(int u, int inEdge) {
    vbDfn[u] = vbLow[u] = ++vbTimer;
    for (int i = head[u]; i; i = e[i].nxt) {
        int v = e[i].to;
        if (i == (inEdge ^ 1)) continue;
        if (!vbDfn[v]) {
            edgeStk[++edgeTop] = i;
            tarjan_vbcc(v, i);
            vbLow[u] = min(vbLow[u], vbLow[v]);
            if (vbLow[v] >= vbDfn[u]) {                  // u 是割点（或根）=> 新点双
                ++total;                                 // 新建方点
                vector<int> vs;
                int x;
                do {
                    x = edgeStk[edgeTop--];
                    vs.push_back(e[x].to);
                    vs.push_back(e[x ^ 1].to);
                } while (x != i);
                sort(vs.begin(), vs.end());
                vs.erase(unique(vs.begin(), vs.end()), vs.end());
                for (int w : vs) { cft[total].push_back(w); cft[w].push_back(total); }
            }
        } else if (vbDfn[v] < vbDfn[u]) {                // 返祖边：压栈（只压一次）
            edgeStk[++edgeTop] = i;
            vbLow[u] = min(vbLow[u], vbDfn[v]);
        }
    }
}
```

**注意**：孤立点（度数为 0 的点）不会形成任何点双，需要单独补一个方点（或直接作为圆点存在）。

### 圆方树的性质（必背）

1. 圆方树是一棵**树**，点数 $= n + (\text{点双个数})$，边数 $= $ 点双内点数之和；
2. **圆点与方点交替**出现（树边总连接一个圆点和一个方点）；
3. **原图中 $u$ 到 $v$ 的所有路径的并 = 圆方树上 $u$ 到 $v$ 的路径上所有方点对应的点双的并**；
4. **原图中 $u \to v$ 的必经点 = 圆方树上 $u \to v$ 路径上除 $u,v$ 外的所有圆点**；
5. 原图中的割点 = 度数 $\ge 2$ 的圆点。

### 经典应用

| 题目 | 做法 |
|---|---|
| [P4320 道路相遇](https://www.luogu.com.cn/problem/P4320) | 圆方树 + LCA，答案 = 路径上圆点数（含端点） |
| [P4630 铁人两项](https://www.luogu.com.cn/problem/P4630) | 圆方树 + 换根/贡献法 + 点权技巧 |
| [P5236 【模板】静态仙人掌](https://www.luogu.com.cn/problem/P5236) | 圆方树 + 树上距离（方点带边权修正） |
| [P3225 矿场搭建](https://www.luogu.com.cn/problem/P3225) | 点双分类讨论：每个点双里数割点个数 |

**P4630 的经典技巧**（方点权值设为 -1 或点双大小）：

把方点权值设为"点双的大小"，圆点权值设为 $-1$，则圆方树上任意两点路径的**点权和**恰好等于原图中两点间"可选中间点"的数量——因为路径上的每个方点贡献它包含的点数，而每个被重复计算的圆点（割点）恰好多算一次，用 $-1$ 抵消。

这类"**点权 → 树上路径和**"的转化是圆方树的灵魂。

---

## 1.6 2-SAT

### 问题

有 $n$ 个布尔变量 $x_1, \dots, x_n$，给定 $m$ 个**子句**，每个子句形如 $(x_a = v_a) \lor (x_b = v_b)$。求一组赋值使所有子句为真。

### 建模：蕴含图

$$(a \lor b) \equiv (\lnot a \to b) \land (\lnot b \to a).$$

对每个变量建两个点：$i$ 表示"$x_i$ 为真"，$i + n$ 表示"$x_i$ 为假"。对每个子句加两条有向边。

**判定**：跑 Tarjan 求 SCC。

- 若存在 $i$ 使 `sccId[i] == sccId[i + n]`：**无解**；
- 否则有解，且赋值规则为

$$x_i = \begin{cases}\text{真}, & \text{sccId}[i] < \text{sccId}[i+n]\\ \text{假}, & \text{否则}\end{cases}$$

**为什么这个赋值对？** Tarjan 求出的 `sccId` 是**逆拓扑序**（编号小 = 拓扑序靠后 = 更"接近汇点"）。选汇点侧的字面量为真，就不会顺着蕴含边推出矛盾。

### 模板

```cpp
// ==================== 2-SAT（Tarjan 版） ====================
// 变量 x_1..x_n：节点 i 表示 "x_i 真"，节点 i+n 表示 "x_i 假"
struct TwoSAT {
    int n;
    vector<vector<int>> g;
    vector<int> dfn, low, stk, sccId;
    vector<char> inStk;
    int timer, top, sccCnt;

    TwoSAT(int n) : n(n), g(2 * n + 2), dfn(2 * n + 2, 0), low(2 * n + 2, 0),
                    stk(2 * n + 2), sccId(2 * n + 2, 0), inStk(2 * n + 2, 0),
                    timer(0), top(0), sccCnt(0) {}

    static int neg(int x, int n) { return x <= n ? x + n : x - n; }

    // 加入子句 (x_a == va) ∨ (x_b == vb)
    void addClause(int a, bool va, int b, bool vb) {
        int na = va ? a : a + n;
        int nb = vb ? b : b + n;
        g[neg(na, n)].push_back(nb);
        g[neg(nb, n)].push_back(na);
    }

    // 便捷写法：x_a 必须为真
    void forceTrue(int a, bool va) { addClause(a, va, a, va); }

    void tarjan(int u) {
        dfn[u] = low[u] = ++timer;
        stk[++top] = u; inStk[u] = 1;
        for (int v : g[u]) {
            if (!dfn[v]) { tarjan(v); low[u] = min(low[u], low[v]); }
            else if (inStk[v]) low[u] = min(low[u], dfn[v]);
        }
        if (low[u] == dfn[u]) {
            ++sccCnt;
            int v;
            do { v = stk[top--]; inStk[v] = 0; sccId[v] = sccCnt; } while (v != u);
        }
    }

    // 求解；返回 false 表示无解，解写入 sol[1..n]
    bool solve(vector<int> &sol) {
        for (int i = 1; i <= 2 * n; i++) if (!dfn[i]) tarjan(i);
        sol.assign(n + 1, 0);
        for (int i = 1; i <= n; i++) {
            if (sccId[i] == sccId[i + n]) return false;
            sol[i] = (sccId[i] < sccId[i + n]) ? 1 : 0;
        }
        return true;
    }
};
```

> **注意**：上面用 `2*n+2` 个点、下标从 1 开始；`neg(i) = i + n` 当 `i <= n`，否则 `i - n`。

### 常见建模套路

| 题目要求 | 建模 |
|---|---|
| $(x_a = v_a) \lor (x_b = v_b)$ | 标准子句 |
| $x_a \to x_b$（蕴含） | $(\lnot x_a \lor x_b)$ |
| $x_a \ne x_b$（不同） | $(x_a \lor x_b) \land (\lnot x_a \lor \lnot x_b)$ |
| $x_a = x_b$ | $(x_a \lor \lnot x_b)\land(\lnot x_a\lor x_b)$ |
| 恰好一个为真（$n$ 中取 1） | 拆成"至多一个"+ 前缀优化建图，或直接用计数 DP |
| $x_a$ 必真 | $(x_a \lor x_a)$ |

**经典题**：[P4782 【模板】2-SAT](https://www.luogu.com.cn/problem/P4782)、[P4171 满汉全席](https://www.luogu.com.cn/problem/P4171)、[P5782 和平委员会](https://www.luogu.com.cn/problem/P5782)、[P3007 牛议会](https://www.luogu.com.cn/problem/P3007)（求字典序最小解，需 DFS 版）。

**求字典序最小解**（当 $n \le 20$ 万但要求最小字典序时不能用 Tarjan 版）：按 $1 \to n$ 顺序贪心，每次尝试把 $x_i$ 设为假，用 DFS 标记蕴含闭包，冲突则回退。复杂度 $O(nm)$ 但实测很快。

---

## 1.7 支配树（简介）

**定义**：在有向图（通常以 $r$ 为源点）中，若从 $r$ 到 $v$ 的所有路径都经过 $u$，则称 $u$ **支配** $v$。除 $v$ 自身外最"近"的支配点称为**直接支配点** $\text{idom}(v)$。所有点连向自己的 idom，构成一棵**支配树**。

**DAG 上的简易做法**（竞赛中够用）：按拓扑序处理，$\text{idom}(v) = \operatorname{LCA}\{\text{所有入边的起点}\}$（在已建的支配树上取 LCA）。

```cpp
// ==================== DAG 支配树（倍增 LCA 在线维护） ====================
// 约定：节点编号 1..n，0 作为"虚拟根"哨兵；源点 r 满足 idom(r) = 0
const int MAXN = 100005;
const int LOG = 18;
int fa[MAXN][LOG], dep[MAXN];
vector<int> dt[MAXN];          // 支配树
int indeg2[MAXN];

int lca_dt(int a, int b) {
    if (!a) return b;
    if (!b) return a;
    if (dep[a] < dep[b]) swap(a, b);
    for (int k = LOG - 1; k >= 0; k--) if (dep[a] - (1 << k) >= dep[b]) a = fa[a][k];
    if (a == b) return a;
    for (int k = LOG - 1; k >= 0; k--) if (fa[a][k] != fa[b][k]) { a = fa[a][k]; b = fa[b][k]; }
    return fa[a][0];
}

void build_dominator_tree(int n, const vector<vector<int>> &pred, const vector<int> &topo) {
    for (int v : topo) {
        int l = 0;
        for (int u : pred[v]) l = lca_dt(l, u);        // 所有前驱的 LCA
        fa[v][0] = l;
        dep[v] = dep[l] + 1;
        for (int k = 1; k < LOG; k++) fa[v][k] = fa[fa[v][k - 1]][k - 1];
        if (l) dt[l].push_back(v);
    }
}
```

**一般有向图的支配树**需要 Lengauer-Tarjan 算法（$O(m\log n)$ 或 $O(m\alpha(n))$），实现约 80 行，属于区域赛难题的工具，建议先把 DAG 版本用熟。

**经典题**：[P2597 [ZJOI2012] 灾难](https://www.luogu.com.cn/problem/P2597)（DAG 支配树模板）、[P5180 【模板】支配树](https://www.luogu.com.cn/problem/P5180)。

---

## 1.8 例题与套路

### 套路 1：看到"互相可达 / 环 / 最长路" → 缩点

> **特征**：有向图 + 环 + 路径/可达性问题。

**步骤**：Tarjan 求 SCC → 缩点 → 在 DAG 上做拓扑序 DP。
**注意**：缩点后要去重边（否则计数类题目会算重）。

### 套路 2：看到"删掉某个点/边会怎样" → 割点/桥

> **特征**：问是否存在"关键点/关键边"；图的"脆弱性"分析。

**步骤**：Tarjan 求割点/桥，在 DFS 过程中累加贡献（如子树大小）。
**注意**：根的判定、重边的处理。

### 套路 3：看到"两点间必经点 / 路径并 / 仙人掌" → 点双 + 圆方树

> **特征**：多次询问两点间路径上的信息；图是仙人掌。

**步骤**：Tarjan 求点双 → 建圆方树 → 树上 LCA/路径和。

### 套路 4：看到"每个限制形如 a 为真或 b 为真" → 2-SAT

> **特征**：每个变量二值，约束是"二者至少一真"。

**步骤**：建蕴含图 → Tarjan → 判 SCC 冲突 → 赋值。

### 套路 5：无向图连通性问题先想"缩点成树"

> **特征**：无向图上的路径查询、树上统计。

**步骤**：e-BCC 缩点（或点双缩点），把图变成树，然后套 LCA/树剖/树上差分。

---

## 1.9 易错点清单

| 编号 | 易错点 | 说明 |
|---|---|---|
| 1 | SCC 中 `else` 分支无条件用 `dfn[v]` 更新 | 必须判 `inStk[v]`，否则已定型的点会污染 low |
| 2 | 缩点后忘记去重边 | DAG 上计数/DP 会算重 |
| 3 | 割点判定漏了根节点特判（`child > 1`） | 常见 WA |
| 4 | 用"父节点编号"跳过回边导致重边误判为桥 | 必须用**边编号** `i == (inEdge ^ 1)` |
| 5 | 边数组从 0 开始编号却用 `i ^ 1` | 必须从 2 开始成对编号 |
| 6 | 点双中"边压栈"漏了返祖边的压栈条件 | 用 `dfn[v] < dfn[u]` 保证每条无向边只压一次 |
| 7 | 圆方树方点与圆点连边未去重 | 同一割点在同一方点内会重复连边 |
| 8 | 2-SAT 赋值方向搞反 | 记牢：`sccId[i] < sccId[i+n]` 取真 |
| 9 | 2-SAT 中变量下标与"取反"映射写错 | 统一用 `neg(x) = x <= n ? x + n : x - n` |
| 10 | 递归深度爆栈 | 链状图用迭代版或手动开栈 |
| 11 | 多组数据忘记清空 `timer`、`top`、`ecnt`、`inStk` | 建议封装成 `init()` |
| 12 | 孤立点在建圆方树时被忽略 | 单独补方点或特判 |

---

## 1.10 模板汇总

```
Tarjan 家族速查（全部 O(n + m)）

有向图：
  tarjan_scc(u)            强连通分量（栈中存点，low[u]==dfn[u] 时弹栈）
  缩点 build_dag(n)        sccId 逆拓扑序，可直接 DP
  TwoSAT                   2-SAT 建模 + 求解

无向图：
  tarjan_cut(u, inEdge)    割点(isCut) 与桥(isBridge)，用边编号处理重边
  dfs_ebcc / build_ebcc    边双连通分量（不走桥 DFS 染色）→ 缩点成树
  tarjan_vbcc(u, inEdge)   点双连通分量（栈中存边）→ 圆方树

DAG：
  build_dominator_tree     支配树（idom(v) = 所有前驱的 LCA）
```

**调试建议**：拿一棵"链"和一张"完全图"分别跑一遍你的 Tarjan，前者测爆栈，后者测 SCC 是否会全部合并成一个。再用随机小图与暴力（枚举删点判连通、Floyd 传递闭包判强连通）对拍。

---

> **下一章预告**：网络流。有了 SCC 与缩点的基本功，最大流/最小割的建模会顺很多——尤其是"最小割 = 最大流"和"最大权闭合子图"这两类经典转化。

---

## 附：本章对拍验证

本章模板已加入 `tools/verify_graph.cpp` 的验证范围（SCC/割点/桥与暴力对拍）。运行方式见 `tools/README.md`。
