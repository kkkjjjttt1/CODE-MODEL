// 第 04 章 模板验证：线段树合并 / 李超树 / 可撤销并查集 / 吉司机线段树 / 左偏树 / 笛卡尔树 / 分块 / 莫队 / 三元环
// 编译: g++ -std=c++17 -O2 -o verify_ds2.exe tools/verify_ds2.cpp
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <algorithm>
#include <numeric>
#include <set>
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

// ==================================================================
// 1. 线段树合并
// ==================================================================
const int SMAXN = 60, SVAL = 16;
int slc[SMAXN * 12], src[SMAXN * 12], ssum[SMAXN * 12];
int sroot[SMAXN], stot;

int newNodeS() { ++stot; slc[stot] = src[stot] = ssum[stot] = 0; return stot; }

int insertS(int &node, int l, int r, int pos) {
    if (!node) node = newNodeS();
    ++ssum[node];
    if (l == r) return node;
    int mid = (l + r) >> 1;
    if (pos <= mid) insertS(slc[node], l, mid, pos);
    else insertS(src[node], mid + 1, r, pos);
    return node;
}
int mergeTreeS(int a, int b, int l, int r) {
    if (!a) return b;
    if (!b) return a;
    if (l == r) { ssum[a] += ssum[b]; return a; }
    int mid = (l + r) >> 1;
    slc[a] = mergeTreeS(slc[a], slc[b], l, mid);
    src[a] = mergeTreeS(src[a], src[b], mid + 1, r);
    ssum[a] = ssum[slc[a]] + ssum[src[a]];
    return a;
}
int queryS(int node, int l, int r, int ql, int qr) {
    if (!node || qr < l || r < ql) return 0;
    if (ql <= l && r <= qr) return ssum[node];
    int mid = (l + r) >> 1;
    return queryS(slc[node], l, mid, ql, qr) + queryS(src[node], mid + 1, r, ql, qr);
}

// ==================================================================
// 2. 李超线段树（维护最大值）
// ==================================================================
const int LCHN = 260;
const double INFDBL = 1e18;

struct Line { double k, b; Line(double k = 0, double b = -INFDBL) : k(k), b(b) {} 
    double at(double x) const { return k * x + b; } };

Line lchTree[LCHN * 4];

void lchInsert(Line nw, int node, int l, int r) {
    int mid = (l + r) >> 1;
    bool leftBetter = nw.at(l) > lchTree[node].at(l);
    bool midBetter = nw.at(mid) > lchTree[node].at(mid);
    if (midBetter) std::swap(nw, lchTree[node]);
    if (l == r) return;
    if (leftBetter != midBetter) lchInsert(nw, node << 1, l, mid);
    else lchInsert(nw, node << 1 | 1, mid + 1, r);
}
double lchQuery(int node, int l, int r, int x) {
    double res = lchTree[node].at(x);
    if (l == r) return res;
    int mid = (l + r) >> 1;
    if (x <= mid) res = std::max(res, lchQuery(node << 1, l, mid, x));
    else res = std::max(res, lchQuery(node << 1 | 1, mid + 1, r, x));
    return res;
}
void lchReset(int node, int l, int r) {
    lchTree[node] = Line();
    if (l == r) return;
    int mid = (l + r) >> 1;
    lchReset(node << 1, l, mid);
    lchReset(node << 1 | 1, mid + 1, r);
}

// ==================================================================
// 3. 可撤销并查集
// ==================================================================
const int UFMAXN = 40;
int ufFa[UFMAXN], ufRk[UFMAXN];
int histU[UFMAXN * 20], histV[UFMAXN * 20], histTop, histVal[UFMAXN * 20];

void initUF(int n) {
    for (int i = 1; i <= n; i++) { ufFa[i] = i; ufRk[i] = 1; }
    histTop = 0;
}
int findUF(int x) { while (ufFa[x] != x) x = ufFa[x]; return x; }
bool uniteUF(int a, int b) {
    int ra = findUF(a), rb = findUF(b);
    if (ra == rb) { histU[++histTop] = 0; histV[histTop] = 0; histVal[histTop] = 0; return false; }
    if (ufRk[ra] < ufRk[rb]) std::swap(ra, rb);
    int add = (ufRk[ra] == ufRk[rb]) ? 1 : 0;
    ufFa[rb] = ra;
    ufRk[ra] += add;
    histU[++histTop] = ra; histV[histTop] = rb; histVal[histTop] = add;
    return true;
}
void undoUF(int stamp) {
    while (histTop > stamp) {
        int ra = histU[histTop], rb = histV[histTop], add = histVal[histTop];
        --histTop;
        if (ra) { ufFa[rb] = rb; ufRk[ra] -= add; }
    }
}

// ==================================================================
// 4. 吉司机线段树（区间 chmin + 区间求和）
// ==================================================================
const int GSMAXN = 70;
const int GINF = 0x3f3f3f3f;
int gmx[GSMAXN * 4], gse[GSMAXN * 4], gcnt[GSMAXN * 4];
ll gsum[GSMAXN * 4];

void gPushUp(int node) {
    int l = node << 1, r = node << 1 | 1;
    gsum[node] = gsum[l] + gsum[r];
    if (gmx[l] == gmx[r]) {
        gmx[node] = gmx[l]; gcnt[node] = gcnt[l] + gcnt[r];
        gse[node] = std::max(gse[l], gse[r]);
    } else if (gmx[l] > gmx[r]) {
        gmx[node] = gmx[l]; gcnt[node] = gcnt[l];
        gse[node] = std::max(gse[l], gmx[r]);
    } else {
        gmx[node] = gmx[r]; gcnt[node] = gcnt[r];
        gse[node] = std::max(gmx[l], gse[r]);
    }
}
void gBuild(int node, int l, int r, const int *a) {
    if (l == r) { gmx[node] = a[l]; gsum[node] = a[l]; gse[node] = -GINF; gcnt[node] = 1; return; }
    int mid = (l + r) >> 1;
    gBuild(node << 1, l, mid, a);
    gBuild(node << 1 | 1, mid + 1, r, a);
    gPushUp(node);
}
void gUpdateMin(int node, int x) {
    if (x >= gmx[node]) return;
    gsum[node] -= (ll)gcnt[node] * (gmx[node] - x);
    gmx[node] = x;
}
void gPushDown(int node) {
    gUpdateMin(node << 1, gmx[node]);
    gUpdateMin(node << 1 | 1, gmx[node]);
}
void gChmin(int node, int l, int r, int ql, int qr, int x) {
    if (qr < l || r < ql || x >= gmx[node]) return;
    if (ql <= l && r <= qr && gse[node] < x) { gUpdateMin(node, x); return; }
    int mid = (l + r) >> 1;
    gPushDown(node);
    gChmin(node << 1, l, mid, ql, qr, x);
    gChmin(node << 1 | 1, mid + 1, r, ql, qr, x);
    gPushUp(node);
}
ll gQuerySum(int node, int l, int r, int ql, int qr) {
    if (qr < l || r < ql) return 0;
    if (ql <= l && r <= qr) return gsum[node];
    int mid = (l + r) >> 1;
    gPushDown(node);
    return gQuerySum(node << 1, l, mid, ql, qr) + gQuerySum(node << 1 | 1, mid + 1, r, ql, qr);
}

// ==================================================================
// 5. 左偏树
// ==================================================================
const int LTMAXN = 60;
int ltL[LTMAXN], ltR[LTMAXN], ltDis[LTMAXN], ltVal[LTMAXN], ltFa[LTMAXN];

void initLT(int n) {
    for (int i = 1; i <= n; i++) { ltL[i] = ltR[i] = ltFa[i] = 0; ltDis[i] = 1; }
}
int findLTRoot(int x) { while (ltFa[x]) x = ltFa[x]; return x; }
int mergeLT(int x, int y) {
    if (!x || !y) return x | y;
    if (ltVal[x] > ltVal[y]) std::swap(x, y);
    ltR[x] = mergeLT(ltR[x], y);
    if (ltDis[ltL[x]] < ltDis[ltR[x]]) std::swap(ltL[x], ltR[x]);
    ltDis[x] = ltDis[ltR[x]] + 1;
    ltFa[ltL[x]] = ltFa[ltR[x]] = x;
    return x;
}
int popLT(int x) {
    int l = ltL[x], r = ltR[x];
    ltVal[x] = -1;
    ltFa[l] = ltFa[r] = 0;
    ltL[x] = ltR[x] = 0; ltDis[x] = 1;
    return mergeLT(l, r);
}

// ==================================================================
// 6. 笛卡尔树
// ==================================================================
const int CTMAXN = 70;
int ctA[CTMAXN], ctL[CTMAXN], ctR[CTMAXN], ctStk[CTMAXN];

int buildCartesian(int n) {
    int top = 0;
    for (int i = 1; i <= n; i++) {
        ctL[i] = ctR[i] = 0;
        int last = 0;
        while (top && ctA[ctStk[top]] > ctA[i]) last = ctStk[top--];
        ctL[i] = last;
        if (top) ctR[ctStk[top]] = i;
        ctStk[++top] = i;
    }
    return ctStk[1];
}

// ==================================================================
// 7. 分块（区间加 + 区间求和）
// ==================================================================
const int BKMAXN = 300;
int blkSize, blkCnt;
int bel[BKMAXN], blkL[BKMAXN], blkR[BKMAXN];
ll bkA[BKMAXN], addTag[BKMAXN], blkSum[BKMAXN];

void buildBlock(int n) {
    blkSize = (int)sqrt((double)n) + 1;
    blkCnt = (n + blkSize - 1) / blkSize;
    for (int i = 1; i <= n; i++) bel[i] = (i - 1) / blkSize + 1;
    for (int b = 1; b <= blkCnt; b++) {
        blkL[b] = (b - 1) * blkSize + 1;
        blkR[b] = std::min(n, b * blkSize);
        addTag[b] = 0; blkSum[b] = 0;
        for (int i = blkL[b]; i <= blkR[b]; i++) blkSum[b] += bkA[i];
    }
}
void rangeAdd(int l, int r, ll v) {
    if (bel[l] == bel[r]) {
        for (int i = l; i <= r; i++) { bkA[i] += v; blkSum[bel[i]] += v; }
        return;
    }
    for (int i = l; i <= blkR[bel[l]]; i++) { bkA[i] += v; blkSum[bel[i]] += v; }
    for (int i = blkL[bel[r]]; i <= r; i++) { bkA[i] += v; blkSum[bel[i]] += v; }
    for (int b = bel[l] + 1; b <= bel[r] - 1; b++) { addTag[b] += v; blkSum[b] += v * (blkR[b] - blkL[b] + 1); }
}
ll rangeSum(int l, int r) {
    ll res = 0;
    if (bel[l] == bel[r]) {
        for (int i = l; i <= r; i++) res += bkA[i] + addTag[bel[i]];
        return res;
    }
    for (int i = l; i <= blkR[bel[l]]; i++) res += bkA[i] + addTag[bel[i]];
    for (int i = blkL[bel[r]]; i <= r; i++) res += bkA[i] + addTag[bel[i]];
    for (int b = bel[l] + 1; b <= bel[r] - 1; b++) res += blkSum[b];
    return res;
}

// ==================================================================
// 8. 莫队（区间不同颜色数）
// ==================================================================
const int MOMAXN = 300;
int moBlockSize, moA[MOMAXN], moCnt[MOMAXN * 4], moDistinct;
struct MoQuery { int l, r, id; };
bool cmpMo(const MoQuery &a, const MoQuery &b) {
    int ba = a.l / moBlockSize, bb = b.l / moBlockSize;
    if (ba != bb) return ba < bb;
    return (ba & 1) ? a.r > b.r : a.r < b.r;
}
void moAdd(int p) { if (++moCnt[moA[p]] == 1) ++moDistinct; }
void moDel(int p) { if (--moCnt[moA[p]] == 0) --moDistinct; }

// ==================================================================
// 9. 三元环计数
// ==================================================================
long long countTriangles(int n, const std::vector<std::pair<int,int>> &edges) {
    std::vector<int> deg(n + 1, 0);
    for (auto &e2 : edges) { ++deg[e2.first]; ++deg[e2.second]; }
    std::vector<std::vector<int>> out(n + 1);
    for (auto &e2 : edges) {
        int u = e2.first, v = e2.second;
        if (deg[u] > deg[v] || (deg[u] == deg[v] && u > v)) std::swap(u, v);
        out[u].push_back(v);
    }
    std::vector<char> mark(n + 1, 0);
    long long cnt = 0;
    for (int u = 1; u <= n; u++) {
        for (int v : out[u]) mark[v] = 1;
        for (int v : out[u])
            for (int w : out[v]) if (mark[w]) ++cnt;
        for (int v : out[u]) mark[v] = 0;
    }
    return cnt;
}

// ==================================================================
int main() {
    setvbuf(stdout, NULL, _IONBF, 0);

    printf("=== 1. 线段树合并 ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 200; iter++) {
            int n = (int)(rng() % 15) + 2;
            vector<vector<int>> tree(n + 1);
            for (int i = 2; i <= n; i++) {
                int p = (int)(rng() % (i - 1)) + 1;
                tree[p].push_back(i); tree[i].push_back(p);
            }
            vector<int> col(n + 1);
            for (int i = 1; i <= n; i++) col[i] = (int)(rng() % SVAL) + 1;

            stot = 0;
            for (int i = 1; i <= n; i++) sroot[i] = 0;
            for (int i = 1; i <= n; i++) insertS(sroot[i], 1, SVAL, col[i]);

            // 以 1 为根 DFS 合并
            vector<int> par(n + 1, 0), ord;
            function<void(int,int)> dfs = [&](int u, int f) {
                par[u] = f;
                for (int v : tree[u]) if (v != f) dfs(v, u);
                ord.push_back(u);
            };
            dfs(1, 0);

            // 关键：必须在"合并进父节点之前"记录子树统计
            // （合并是原地复用节点，被合并掉的树之后不能再查询）
            vector<vector<int>> rec(n + 1, vector<int>(SVAL + 1, 0));
            for (int u : ord) {
                for (int c = 1; c <= SVAL; c++) rec[u][c] = queryS(sroot[u], 1, SVAL, c, c);
                if (par[u]) sroot[par[u]] = mergeTreeS(sroot[par[u]], sroot[u], 1, SVAL);
            }
            // 验证：每个子树内各颜色的计数
            for (int u = 1; u <= n; u++) {
                vector<int> sub;
                function<void(int,int)> collect = [&](int x, int f) {
                    sub.push_back(x);
                    for (int v : tree[x]) if (v != f) collect(v, x);
                };
                collect(u, par[u]);
                for (int c = 1; c <= SVAL; c++) {
                    int expect = 0;
                    for (int x : sub) if (col[x] == c) ++expect;
                    int got = rec[u][c];
                    if (got != expect) { ok = false; printf("    u=%d c=%d got=%d expect=%d\n", u, c, got, expect); break; }
                }
                if (!ok) break;
            }
            if (!ok) break;
        }
        check(ok, "线段树合并后子树各颜色计数与暴力一致（200 组随机树）");
    }

    printf("=== 2. 李超线段树 ===\n");
    {
        bool ok = true;
        const int X0 = 0, X1 = 100;
        for (int iter = 0; iter < 200; iter++) {
            lchReset(1, X0, X1);
            int m = (int)(rng() % 12) + 1;
            vector<Line> lines;
            for (int i = 0; i < m; i++) {
                double k = (double)((int)(rng() % 21) - 10);
                double b = (double)((int)(rng() % 101) - 50);
                lines.push_back(Line(k, b));
                lchInsert(Line(k, b), 1, X0, X1);
            }
            for (int x = X0; x <= X1; x++) {
                double expect = -INFDBL;
                for (auto &ln : lines) expect = std::max(expect, ln.at(x));
                double got = lchQuery(1, X0, X1, x);
                if (fabs(got - expect) > 1e-6) {
                    ok = false;
                    printf("    x=%d got=%.3f expect=%.3f\n", x, got, expect);
                    break;
                }
            }
            if (!ok) break;
        }
        check(ok, "李超线段树单点最值与暴力一致（200 组随机直线）");
    }

    printf("=== 3. 可撤销并查集 ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 200; iter++) {
            int n = (int)(rng() % 12) + 2;
            initUF(n);
            vector<int> snapshots;          // 记录每个时间点的 histTop
            vector<vector<pair<int,int>>> edgeHistory;   // 每个时间点的边集（暴力用）

            vector<pair<int,int>> curEdges;
            for (int step = 0; step < 30; step++) {
                int op = (int)(rng() % 3);
                if (op == 0) {
                    int a = (int)(rng() % n) + 1, b = (int)(rng() % n) + 1;
                    uniteUF(a, b);
                    curEdges.push_back({a, b});
                } else if (op == 1 && !snapshots.empty()) {
                    // 撤销到某个历史点
                    int t = (int)(rng() % snapshots.size());
                    undoUF(snapshots[t]);
                    curEdges = edgeHistory[t];
                    snapshots.resize(t + 1);
                    edgeHistory.resize(t + 1);
                    continue;
                }
                snapshots.push_back(histTop);
                edgeHistory.push_back(curEdges);
                // 与暴力比较连通性
                for (int a = 1; a <= n && ok; a++)
                    for (int b = a + 1; b <= n; b++) {
                        vector<int> pa(n + 1);
                        for (int i = 1; i <= n; i++) pa[i] = i;
                        function<int(int)> f = [&](int x) { return pa[x] == x ? x : pa[x] = f(pa[x]); };
                        for (auto &e2 : curEdges) pa[f(e2.first)] = f(e2.second);
                        bool bruteConn = (f(a) == f(b));
                        bool gotConn = (findUF(a) == findUF(b));
                        if (bruteConn != gotConn) { ok = false; printf("    step=%d (%d,%d) 不一致\n", step, a, b); break; }
                    }
                if (!ok) break;
            }
            if (!ok) break;
        }
        check(ok, "可撤销并查集（合并+回滚）连通性与暴力一致（200 组随机操作序列）");
    }

    printf("=== 4. 吉司机线段树（区间 chmin + 求和） ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 300; iter++) {
            int n = (int)(rng() % 30) + 2;
            static int arr[GSMAXN], brute[GSMAXN];
            for (int i = 1; i <= n; i++) { arr[i] = (int)(rng() % 50); brute[i] = arr[i]; }
            gBuild(1, 1, n, arr);
            int q = (int)(rng() % 20) + 1;
            for (int t = 0; t < q; t++) {
                int l = (int)(rng() % n) + 1, r = (int)(rng() % n) + 1;
                if (l > r) std::swap(l, r);
                if (rng() % 2) {
                    int x = (int)(rng() % 50);
                    gChmin(1, 1, n, l, r, x);
                    for (int i = l; i <= r; i++) brute[i] = std::min(brute[i], x);
                } else {
                    ll expect = 0;
                    for (int i = l; i <= r; i++) expect += brute[i];
                    ll got = gQuerySum(1, 1, n, l, r);
                    if (got != expect) { ok = false; printf("    n=%d [%d,%d] got=%lld expect=%lld\n", n, l, r, got, expect); break; }
                }
            }
            if (!ok) break;
        }
        check(ok, "吉司机线段树区间 chmin + 求和与暴力一致（300 组随机测试）");
    }

    printf("=== 5. 左偏树（可并堆） ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 200 && ok; iter++) {
            int n = (int)(rng() % 20) + 2;
            initLT(n);
            for (int i = 1; i <= n; i++) ltVal[i] = (int)(rng() % 100) + 1;

            vector<int> dsu(n + 1);
            for (int i = 1; i <= n; i++) dsu[i] = i;
            function<int(int)> dfind = [&](int x) { return dsu[x] == x ? x : dsu[x] = dfind(dsu[x]); };
            vector<multiset<int>> ms(n + 1);
            vector<int> heapRoot(n + 1, 0);
            for (int i = 1; i <= n; i++) { ms[i].insert(ltVal[i]); heapRoot[i] = i; }

            for (int t = 0; t < 30 && ok; t++) {
                int a = (int)(rng() % n) + 1, b = (int)(rng() % n) + 1;
                int ha = dfind(a), hb = dfind(b);
                if (ha == hb) continue;
                int newRoot = mergeLT(heapRoot[ha], heapRoot[hb]);
                dsu[hb] = ha;
                for (int x : ms[hb]) ms[ha].insert(x);
                ms[hb].clear();
                heapRoot[ha] = newRoot;
                if (ltVal[newRoot] != *ms[ha].begin()) {
                    ok = false;
                    printf("    合并后堆顶错误 newRoot=%d val=%d expect=%d\n",
                           newRoot, ltVal[newRoot], *ms[ha].begin());
                    break;
                }
                if (rng() % 2 && !ms[ha].empty()) {          // 随机弹出堆顶
                    int expectedMin = *ms[ha].begin();
                    if (ltVal[newRoot] != expectedMin) { ok = false; break; }
                    heapRoot[ha] = popLT(newRoot);
                    ms[ha].erase(ms[ha].begin());
                    if (heapRoot[ha] != 0 && !ms[ha].empty() &&
                        ltVal[heapRoot[ha]] != *ms[ha].begin()) {
                        ok = false; printf("    弹出后堆顶错误\n"); break;
                    }
                }
            }
        }
        check(ok, "左偏树合并/弹出后堆顶与 multiset 一致（200 组随机操作）");
    }

    printf("=== 6. 笛卡尔树 ===\n");
    {
        bool ok = true, okSum = true;
        for (int iter = 0; iter < 300; iter++) {
            int n = (int)(rng() % 30) + 1;
            for (int i = 1; i <= n; i++) ctA[i] = (int)(rng() % 50) + 1;
            int root = buildCartesian(n);
            // 验证中序遍历 = 原序列
            vector<int> inord;
            function<void(int)> inorder = [&](int u) {
                if (!u) return;
                inorder(ctL[u]); inord.push_back(u); inorder(ctR[u]);
            };
            inorder(root);
            if ((int)inord.size() != n) { ok = false; printf("    中序长度错\n"); break; }
            for (int i = 0; i < n; i++) if (inord[i] != i + 1) { ok = false; printf("    中序不是原序列\n"); break; }
            // 验证堆性质（小根）
            function<bool(int)> heapOk = [&](int u) {
                if (!u) return true;
                if (ctL[u] && ctA[ctL[u]] < ctA[u]) return false;
                if (ctR[u] && ctA[ctR[u]] < ctA[u]) return false;
                return heapOk(ctL[u]) && heapOk(ctR[u]);
            };
            if (!heapOk(root)) { ok = false; printf("    堆性质不满足\n"); break; }
            // 用笛卡尔树求"所有子区间最小值之和"并与暴力比较
            function<pair<ll,int>(int)> solveSub = [&](int u) -> pair<ll,int> {
                if (!u) return {0, 0};
                auto [sumL, szL] = solveSub(ctL[u]);
                auto [sumR, szR] = solveSub(ctR[u]);
                ll contrib = (ll)ctA[u] * (szL + 1) * (szR + 1);
                return {sumL + sumR + contrib, szL + szR + 1};
            };
            ll got = solveSub(root).first;
            ll expect = 0;
            for (int l = 1; l <= n; l++) {
                int mn = INT_MAX;
                for (int r = l; r <= n; r++) { mn = std::min(mn, ctA[r]); expect += mn; }
            }
            if (got != expect) { okSum = false; printf("    n=%d got=%lld expect=%lld\n", n, got, expect); break; }
        }
        check(ok, "笛卡尔树的中序 = 原序列且满足小根堆性质");
        check(okSum, "用笛卡尔树求所有子区间最小值之和与暴力一致（300 组）");
    }

    printf("=== 7. 分块（区间加 + 区间求和） ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 300; iter++) {
            int n = (int)(rng() % 150) + 1;
            static ll brute[BKMAXN];
            for (int i = 1; i <= n; i++) { bkA[i] = (ll)(rng() % 100); brute[i] = bkA[i]; }
            buildBlock(n);
            int q = (int)(rng() % 40) + 1;
            for (int t = 0; t < q; t++) {
                int l = (int)(rng() % n) + 1, r = (int)(rng() % n) + 1;
                if (l > r) std::swap(l, r);
                if (rng() % 2) {
                    ll v = (ll)((int)(rng() % 21) - 10);
                    rangeAdd(l, r, v);
                    for (int i = l; i <= r; i++) brute[i] += v;
                } else {
                    ll expect = 0;
                    for (int i = l; i <= r; i++) expect += brute[i];
                    ll got = rangeSum(l, r);
                    if (got != expect) { ok = false; printf("    n=%d [%d,%d] got=%lld expect=%lld\n", n, l, r, got, expect); break; }
                }
            }
            if (!ok) break;
        }
        check(ok, "分块区间加/区间求和与暴力一致（300 组随机测试）");
    }

    printf("=== 8. 莫队（区间不同颜色数） ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 200; iter++) {
            int n = (int)(rng() % 150) + 1;
            for (int i = 1; i <= n; i++) moA[i] = (int)(rng() % 10) + 1;
            int q = (int)(rng() % 50) + 1;
            vector<MoQuery> qs(q);
            vector<int> expect(q);
            for (int i = 0; i < q; i++) {
                int l = (int)(rng() % n) + 1, r = (int)(rng() % n) + 1;
                if (l > r) std::swap(l, r);
                qs[i] = {l, r, i};
                vector<char> seen(20, 0);
                int c = 0;
                for (int p = l; p <= r; p++) if (!seen[moA[p]]) { seen[moA[p]] = 1; ++c; }
                expect[i] = c;
            }
            moBlockSize = std::max(1, (int)(n / std::sqrt((double)q)));
            vector<MoQuery> order = qs;
            std::sort(order.begin(), order.end(), cmpMo);
            fill(moCnt, moCnt + MOMAXN * 4, 0);
            moDistinct = 0;
            int curL = 1, curR = 0;
            vector<int> got(q);
            for (auto &qq : order) {
                while (curL > qq.l) moAdd(--curL);
                while (curR < qq.r) moAdd(++curR);
                while (curL < qq.l) moDel(curL++);
                while (curR > qq.r) moDel(curR--);
                got[qq.id] = moDistinct;
            }
            for (int i = 0; i < q; i++)
                if (got[i] != expect[i]) { ok = false; printf("    询问 %d got=%d expect=%d\n", i, got[i], expect[i]); break; }
            if (!ok) break;
        }
        check(ok, "莫队求区间不同颜色数与暴力一致（200 组随机询问）");
    }

    printf("=== 9. 三元环计数 ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 300; iter++) {
            int n = (int)(rng() % 14) + 3;
            vector<pair<int,int>> edges;
            for (int i = 1; i <= n; i++)
                for (int j = i + 1; j <= n; j++)
                    if (rng() % 3 == 0) edges.push_back({i, j});
            ll got = countTriangles(n, edges);
            // 暴力：用邻接矩阵枚举三元组
            vector<vector<char>> adj(n + 1, vector<char>(n + 1, 0));
            for (auto &e2 : edges) { adj[e2.first][e2.second] = 1; adj[e2.second][e2.first] = 1; }
            ll expect = 0;
            for (int i = 1; i <= n; i++)
                for (int j = i + 1; j <= n; j++)
                    for (int k = j + 1; k <= n; k++)
                        if (adj[i][j] && adj[j][k] && adj[i][k]) ++expect;
            if (got != expect) { ok = false; printf("    n=%d got=%lld expect=%lld\n", n, got, expect); break; }
        }
        check(ok, "三元环计数与暴力枚举一致（300 组随机图）");
    }

    printf("\n===== 通过 %d 项，失败 %d 项 =====\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
