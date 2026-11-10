// 第 01 章 模板验证：SCC / 割点 / 桥 / 边双 / 点双+圆方树 / 2-SAT / DAG 支配树
// 编译: g++ -std=c++17 -O2 -o verify_graph.exe tools/verify_graph.cpp
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <queue>
#include <random>
#include <functional>

using namespace std;
using ll = long long;

static int g_fail = 0, g_pass = 0;
static void check(bool cond, const char *name) {
    if (cond) { ++g_pass; printf("  [ok]   %s\n", name); }
    else      { ++g_fail; printf("  [FAIL] %s\n", name); }
}
static mt19937 rng(20240607);

const int MAXN = 14, MAXM = 80;

// ==================================================================
// 教程模板 1：有向图 SCC
// ==================================================================
vector<int> gDir[MAXN];
int dfn[MAXN], low[MAXN], timer_;
int stk_[MAXN], top_;
bool inStk[MAXN];
int sccId[MAXN], sccCnt, sccSize[MAXN];

void tarjan_scc(int u) {
    dfn[u] = low[u] = ++timer_;
    stk_[++top_] = u; inStk[u] = true;
    for (int v : gDir[u]) {
        if (!dfn[v]) {
            tarjan_scc(v);
            low[u] = min(low[u], low[v]);
        } else if (inStk[v]) {
            low[u] = min(low[u], dfn[v]);
        }
    }
    if (low[u] == dfn[u]) {
        ++sccCnt;
        int v;
        do {
            v = stk_[top_--];
            inStk[v] = false;
            sccId[v] = sccCnt;
            ++sccSize[sccCnt];
        } while (v != u);
    }
}

// ==================================================================
// 教程模板 2：无向图 割点 / 桥
// ==================================================================
struct Edge { int to, nxt; } e[MAXM * 2];
int head[MAXN], ecnt = 1;
void addEdge(int u, int v) {
    e[++ecnt] = {v, head[u]}; head[u] = ecnt;
    e[++ecnt] = {u, head[v]}; head[v] = ecnt;
}
int dfnU[MAXN], lowU[MAXN], timerU;
bool isCut[MAXN], isBridge[MAXM * 2];

void tarjan_cut(int u, int inEdge) {
    dfnU[u] = lowU[u] = ++timerU;
    int child = 0;
    for (int i = head[u]; i; i = e[i].nxt) {
        int v = e[i].to;
        if (i == (inEdge ^ 1)) continue;
        if (!dfnU[v]) {
            ++child;
            tarjan_cut(v, i);
            lowU[u] = min(lowU[u], lowU[v]);
            if (lowU[v] > dfnU[u]) isBridge[i] = isBridge[i ^ 1] = true;
            if (lowU[v] >= dfnU[u] && inEdge != 0) isCut[u] = true;
        } else if (dfnU[v] < dfnU[u]) {
            lowU[u] = min(lowU[u], dfnU[v]);
        }
    }
    if (inEdge == 0 && child > 1) isCut[u] = true;
}

// ==================================================================
// 教程模板 3：点双 + 圆方树（用独立的存储，避免与上面冲突）
// ==================================================================
struct Edge2 { int to, nxt; } e2[MAXM * 2];
int head2[MAXN], ecnt2 = 1;
void addEdge2(int u, int v) {
    e2[++ecnt2] = {v, head2[u]}; head2[u] = ecnt2;
    e2[++ecnt2] = {u, head2[v]}; head2[v] = ecnt2;
}
int vbDfn[MAXN], vbLow[MAXN], vbTimer;
int edgeStk[MAXM * 2], edgeTop;
int n0_;
int total_;
vector<int> cft[MAXN * 2 + 5];

void tarjan_vbcc(int u, int inEdge) {
    vbDfn[u] = vbLow[u] = ++vbTimer;
    for (int i = head2[u]; i; i = e2[i].nxt) {
        int v = e2[i].to;
        if (i == (inEdge ^ 1)) continue;
        if (!vbDfn[v]) {
            edgeStk[++edgeTop] = i;
            tarjan_vbcc(v, i);
            vbLow[u] = min(vbLow[u], vbLow[v]);
            if (vbLow[v] >= vbDfn[u]) {
                ++total_;
                vector<int> vs;
                int x;
                do {
                    x = edgeStk[edgeTop--];
                    vs.push_back(e2[x].to);
                    vs.push_back(e2[x ^ 1].to);
                } while (x != i);
                sort(vs.begin(), vs.end());
                vs.erase(unique(vs.begin(), vs.end()), vs.end());
                for (int w : vs) { cft[total_].push_back(w); cft[w].push_back(total_); }
            }
        } else if (vbDfn[v] < vbDfn[u]) {
            edgeStk[++edgeTop] = i;
            vbLow[u] = min(vbLow[u], vbDfn[v]);
        }
    }
}

// ==================================================================
// 教程模板 4：2-SAT
// ==================================================================
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

    void addClause(int a, bool va, int b, bool vb) {
        int na = va ? a : a + n;
        int nb = vb ? b : b + n;
        g[neg(na, n)].push_back(nb);
        g[neg(nb, n)].push_back(na);
    }
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

// ==================================================================
// 教程模板 5：DAG 支配树
// ==================================================================
const int LOGD = 8;
int fa[MAXN][LOGD], dep[MAXN];
vector<int> dt[MAXN];

int lca_dt(int a, int b) {
    if (!a) return b;
    if (!b) return a;
    if (dep[a] < dep[b]) swap(a, b);
    for (int k = LOGD - 1; k >= 0; k--) if (dep[a] - (1 << k) >= dep[b]) a = fa[a][k];
    if (a == b) return a;
    for (int k = LOGD - 1; k >= 0; k--) if (fa[a][k] != fa[b][k]) { a = fa[a][k]; b = fa[b][k]; }
    return fa[a][0];
}

void build_dominator_tree(int n, const vector<vector<int>> &pred, const vector<int> &topo) {
    for (int i = 0; i <= n; i++) {   // 清空
        for (int k = 0; k < LOGD; k++) fa[i][k] = 0;
        dep[i] = 0; dt[i].clear();
    }
    for (int v : topo) {
        int l = 0;
        for (int u : pred[v]) l = lca_dt(l, u);
        fa[v][0] = l;
        dep[v] = dep[l] + 1;
        for (int k = 1; k < LOGD; k++) fa[v][k] = fa[fa[v][k - 1]][k - 1];
        if (l) dt[l].push_back(v);
    }
}

// ==================================================================
// 暴力参照
// ==================================================================
int n, m;
vector<pair<int,int>> edges;
bool reachM[MAXN][MAXN];

// 无向图：计算连通块数（可禁用某个点 / 某条边）
int countComp(int banNode, int banEdge) {
    vector<int> pa(n + 1);
    for (int i = 1; i <= n; i++) pa[i] = i;
    function<int(int)> find = [&](int x) { return pa[x] == x ? x : pa[x] = find(pa[x]); };
    for (int i = 0; i < (int)edges.size(); i++) {
        if (i == banEdge) continue;
        int u = edges[i].first, v = edges[i].second;
        if (u == banNode || v == banNode) continue;
        pa[find(u)] = find(v);
    }
    int c = 0;
    for (int i = 1; i <= n; i++) if (i != banNode && find(i) == i) ++c;
    return c;
}

// 有向图：传递闭包（用于 SCC 判定）
void floydReach() {
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++) reachM[i][j] = false;
    for (int u = 1; u <= n; u++) for (int v : gDir[u]) reachM[u][v] = true;
    for (int i = 1; i <= n; i++) reachM[i][i] = true;
    for (int k = 1; k <= n; k++)
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                if (reachM[i][k] && reachM[k][j]) reachM[i][j] = true;
}

// 划分比较工具：把 comp 编号转成"划分的规范形式"，两两比较是否等价
string canonicalPartition(const int *id, int nn) {
    vector<int> firstSeen;
    string s;
    for (int i = 1; i <= nn; i++) {
        int c = id[i];
        bool found = false;
        for (size_t j = 0; j < firstSeen.size(); j++)
            if (firstSeen[j] == c) { s += char('a' + (int)j); found = true; break; }
        if (!found) { firstSeen.push_back(c); s += char('a' + (int)firstSeen.size() - 1); }
    }
    return s;
}

// ==================================================================
int main() {
    printf("=== 1. 强连通分量 SCC ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 400; iter++) {
            n = (int)(rng() % 8) + 1;
            for (int i = 1; i <= n; i++) gDir[i].clear();
            m = (int)(rng() % 14);
            for (int i = 0; i < m; i++) {
                int u = (int)(rng() % n) + 1, v = (int)(rng() % n) + 1;
                gDir[u].push_back(v);
            }
            // 我的 Tarjan
            timer_ = top_ = sccCnt = 0;
            for (int i = 1; i <= n; i++) { dfn[i] = low[i] = 0; inStk[i] = false; sccId[i] = 0; sccSize[i] = 0; }
            for (int i = 1; i <= n; i++) if (!dfn[i]) tarjan_scc(i);
            // 暴力：传递闭包
            floydReach();
            vector<int> bruteId(n + 1, 0);
            int bc = 0;
            vector<int> assigned(n + 1, 0);
            for (int i = 1; i <= n; i++) {
                if (assigned[i]) continue;
                ++bc;
                for (int j = 1; j <= n; j++)
                    if (reachM[i][j] && reachM[j][i]) { assigned[j] = 1; bruteId[j] = bc; }
            }
            // 比较划分是否一致（SCC 的编号顺序无意义）
            if (canonicalPartition(sccId, n) != canonicalPartition(bruteId.data(), n)) {
                ok = false;
                printf("    n=%d m=%d 划分不一致\n", n, m);
                break;
            }
            // 检查 SCC 大小之和
            int sum = 0;
            for (int i = 1; i <= sccCnt; i++) sum += sccSize[i];
            if (sum != n) { ok = false; break; }
        }
        check(ok, "Tarjan SCC 划分与传递闭包暴力一致（400 组随机图）");

        // 缩点后 sccId 是逆拓扑序：若存在 u->v 的边且跨 SCC，则 sccId[u] > sccId[v]
        bool ok2 = true;
        for (int iter = 0; iter < 200; iter++) {
            n = (int)(rng() % 8) + 1;
            for (int i = 1; i <= n; i++) gDir[i].clear();
            m = (int)(rng() % 14);
            for (int i = 0; i < m; i++) {
                int u = (int)(rng() % n) + 1, v = (int)(rng() % n) + 1;
                gDir[u].push_back(v);
            }
            timer_ = top_ = sccCnt = 0;
            for (int i = 1; i <= n; i++) { dfn[i] = low[i] = 0; inStk[i] = false; sccId[i] = 0; }
            for (int i = 1; i <= n; i++) if (!dfn[i]) tarjan_scc(i);
            for (int u = 1; u <= n; u++)
                for (int v : gDir[u])
                    if (sccId[u] != sccId[v] && sccId[u] < sccId[v]) ok2 = false;
        }
        check(ok2, "sccId 满足「跨 SCC 的边必然从大编号指向小编号」（逆拓扑序）");
    }

    printf("=== 2. 割点与桥 ===\n");
    {
        bool okCut = true, okBridge = true;
        for (int iter = 0; iter < 400; iter++) {
            n = (int)(rng() % 8) + 1;
            // 先造一棵随机树保证连通性，再加随机边（允许重边）
            edges.clear();
            for (int i = 2; i <= n; i++) edges.push_back({(int)(rng() % (i - 1)) + 1, i});
            int extra = (int)(rng() % 6);
            for (int i = 0; i < extra; i++)
                edges.push_back({(int)(rng() % n) + 1, (int)(rng() % n) + 1});
            // 加自环测试
            if (rng() % 4 == 0 && n >= 1) edges.push_back({1, 1});

            // 建图
            ecnt = 1;
            for (int i = 1; i <= n; i++) head[i] = 0;
            for (int i = 0; i < MAXM * 2; i++) { e[i].to = 0; e[i].nxt = 0; }
            for (int i = 0; i < (int)edges.size(); i++)
                if (edges[i].first != edges[i].second) addEdge(edges[i].first, edges[i].second);

            for (int i = 0; i < MAXN; i++) { dfnU[i] = lowU[i] = 0; isCut[i] = false; }
            for (int i = 0; i < MAXM * 2; i++) isBridge[i] = false;
            timerU = 0;
            for (int i = 1; i <= n; i++) if (!dfnU[i]) tarjan_cut(i, 0);

            int base = countComp(0, -1);
            for (int u = 1; u <= n; u++) {
                bool bruteCut = countComp(u, -1) > base;
                if (bruteCut != isCut[u]) { okCut = false; printf("    割点 n=%d u=%d 算法=%d 暴力=%d\n", n, u, (int)isCut[u], (int)bruteCut); break; }
            }
            if (!okCut) break;
            // 桥：把链式前向星里成对的边映射回 edges
            for (int idx = 0; idx < (int)edges.size() && okBridge; idx++) {
                if (edges[idx].first == edges[idx].second) continue;   // 自环不是桥
                bool bruteBridge = countComp(-1, idx) > base;
                // 找算法结果：遍历 u 的边找到这条边
                bool algBridge = false;
                for (int i = 2; i <= ecnt; i++) {
                    int a = e[i].to, b = e[i ^ 1].to;
                    if ((a == edges[idx].first && b == edges[idx].second) ||
                        (a == edges[idx].second && b == edges[idx].first)) {
                        // 重边情形：只要有一条同端点的边不是桥，就都不是桥
                        if (!isBridge[i]) { algBridge = false; break; }
                        algBridge = true;
                    }
                }
                if (bruteBridge != algBridge) {
                    okBridge = false;
                    printf("    桥 n=%d edge(%d,%d) 算法=%d 暴力=%d\n", n, edges[idx].first, edges[idx].second, (int)algBridge, (int)bruteBridge);
                }
            }
            if (!okBridge) break;
        }
        check(okCut, "割点判定与「删点后连通块数增加」暴力一致（400 组随机图）");
        check(okBridge, "桥判定与「删边后连通块数增加」暴力一致（含重边/自环）");
    }

    printf("=== 3. 点双 + 圆方树 ===\n");
    {
        bool okTree = true, okAlt = true, okCutEq = true;
        for (int iter = 0; iter < 300; iter++) {
            n = (int)(rng() % 8) + 2;
            edges.clear();
            for (int i = 2; i <= n; i++) edges.push_back({(int)(rng() % (i - 1)) + 1, i});
            int extra = (int)(rng() % 4);
            for (int i = 0; i < extra; i++)
                edges.push_back({(int)(rng() % n) + 1, (int)(rng() % n) + 1});

            ecnt2 = 1;
            for (int i = 1; i <= n; i++) head2[i] = 0;
            for (int i = 0; i < MAXN; i++) { vbDfn[i] = vbLow[i] = 0; }
            vbTimer = edgeTop = 0;
            n0_ = n; total_ = n;
            for (int i = 1; i <= 2 * MAXN + 4; i++) cft[i].clear();
            for (auto &pr : edges)
                if (pr.first != pr.second) addEdge2(pr.first, pr.second);
            for (int i = 1; i <= n; i++) if (!vbDfn[i]) tarjan_vbcc(i, 0);

            // 1) 圆方树是树：边数 = 点数 - 1（连通图）
            int treeEdges = 0, treeNodes = total_;
            for (int i = 1; i <= total_; i++) treeEdges += (int)cft[i].size();
            treeEdges /= 2;
            if (treeEdges != treeNodes - 1) {
                okTree = false;
                printf("    n=%d 圆方树边数=%d 点数=%d\n", n, treeEdges, treeNodes);
                break;
            }
            // 2) 圆点与方点交替
            for (int u = 1; u <= total_; u++)
                for (int v : cft[u])
                    if (!((u <= n) != (v <= n))) { okAlt = false; break; }
            if (!okAlt) break;

            // 3) 割点 = 圆方树中度数 >= 2 的圆点
            // 先用暴力求割点
            ecnt = 1;
            for (int i = 1; i <= n; i++) head[i] = 0;
            for (auto &pr : edges) if (pr.first != pr.second) addEdge(pr.first, pr.second);
            for (int i = 0; i < MAXN; i++) { dfnU[i] = lowU[i] = 0; isCut[i] = false; }
            timerU = 0;
            for (int i = 1; i <= n; i++) if (!dfnU[i]) tarjan_cut(i, 0);
            int base = countComp(0, -1);
            for (int u = 1; u <= n; u++) {
                bool bruteCut = countComp(u, -1) > base;
                bool fromTree = ((int)cft[u].size() >= 2);
                if (bruteCut != fromTree) {
                    okCutEq = false;
                    printf("    n=%d u=%d 暴力割点=%d 圆方树度=%d\n", n, u, (int)bruteCut, (int)cft[u].size());
                    break;
                }
            }
            if (!okCutEq) break;
        }
        check(okTree, "圆方树是一棵树（边数 = 点数 - 1）");
        check(okAlt, "圆方树中圆点与方点交替出现");
        check(okCutEq, "原图割点 == 圆方树中度数 ≥ 2 的圆点");
    }

    printf("=== 4. 2-SAT ===\n");
    {
        bool ok = true;
        int tested = 0, satCnt = 0;
        for (int iter = 0; iter < 600; iter++) {
            int nn = (int)(rng() % 5) + 1;
            int mm = (int)(rng() % 9);
            vector<array<int,4>> cls;
            for (int i = 0; i < mm; i++) {
                int a = (int)(rng() % nn) + 1, b = (int)(rng() % nn) + 1;
                int va = (int)(rng() % 2), vb = (int)(rng() % 2);
                cls.push_back({a, va, b, vb});
            }
            // 暴力枚举
            bool bruteSat = false;
            vector<int> bruteSol(nn + 1, 0);
            for (int mask = 0; mask < (1 << nn); mask++) {
                bool good = true;
                for (auto &c : cls) {
                    int xa = (mask >> (c[0] - 1)) & 1;
                    int xb = (mask >> (c[2] - 1)) & 1;
                    if (!((xa == c[1]) || (xb == c[3]))) { good = false; break; }
                }
                if (good) {
                    bruteSat = true;
                    for (int i = 1; i <= nn; i++) bruteSol[i] = (mask >> (i - 1)) & 1;
                    break;
                }
            }
            // 我的求解器
            TwoSAT ts(nn);
            for (auto &c : cls) ts.addClause(c[0], c[1] == 1, c[2], c[3] == 1);
            vector<int> sol;
            bool mySat = ts.solve(sol);
            tested++;
            if (bruteSat) satCnt++;
            if (bruteSat != mySat) {
                ok = false;
                printf("    n=%d m=%d 算法=%d 暴力=%d\n", nn, mm, (int)mySat, (int)bruteSat);
                break;
            }
            if (mySat) {
                for (auto &c : cls) {
                    if (!((sol[c[0]] == c[1]) || (sol[c[2]] == c[3]))) { ok = false; break; }
                }
                if (!ok) { printf("    解不满足子句\n"); break; }
            }
        }
        check(ok, "2-SAT 判定与赋值均与暴力枚举一致（含随机可满足/不可满足实例）");
        printf("    （%d 组实例，其中 %d 组可满足）\n", tested, satCnt);
    }

    printf("=== 5. DAG 支配树 ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 300; iter++) {
            int nn = (int)(rng() % 7) + 2;          // 点 1..nn，1 为源点，0 是虚拟根哨兵
            vector<vector<int>> pred(nn + 1), succ(nn + 1);
            for (int v = 2; v <= nn; v++) {          // 保证源点 1 能到达所有点
                int u = (int)(rng() % (v - 1)) + 1;
                pred[v].push_back(u);
                succ[u].push_back(v);
            }
            int extra = (int)(rng() % 6);
            for (int i = 0; i < extra; i++) {
                int u = (int)(rng() % nn) + 1, v = (int)(rng() % nn) + 1;
                if (u < v) { pred[v].push_back(u); succ[u].push_back(v); }
            }
            vector<int> topo;
            for (int i = 1; i <= nn; i++) topo.push_back(i);
            build_dominator_tree(nn, pred, topo);

            // 暴力：对每个 v，D(v) = {u : 删除 u 后 1 无法到达 v}
            for (int v = 2; v <= nn; v++) {
                int ancMask = 0;                     // 支配树上除 v 外的祖先（含源点 1）
                for (int x = fa[v][0]; x; x = fa[x][0]) ancMask |= (1 << x);
                int bruteMask = 0;
                for (int ban = 1; ban <= nn; ban++) {
                    if (ban == v) continue;          // 只关心除自身外的支配点
                    vector<char> vis(nn + 1, 0);
                    queue<int> q;
                    if (ban != 1) { vis[1] = 1; q.push(1); }
                    while (!q.empty()) {
                        int u = q.front(); q.pop();
                        for (int w : succ[u])
                            if (w != ban && !vis[w]) { vis[w] = 1; q.push(w); }
                    }
                    if (!vis[v]) bruteMask |= (1 << ban);
                }
                if (ancMask != bruteMask) {
                    ok = false;
                    printf("    nn=%d v=%d 支配树祖先=%d 暴力=%d\n", nn, v, ancMask, bruteMask);
                    break;
                }
            }
            if (!ok) break;
        }
        check(ok, "DAG 支配树的祖先集合 == 暴力定义的支配点集合");
    }

    printf("\n===== 通过 %d 项，失败 %d 项 =====\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
