// 第 03 章 模板验证：点分治 / DSU on tree / 长链剖分 / 虚树 / 树哈希
// 编译: g++ -std=c++17 -O2 -o verify_ds.exe tools/verify_ds.cpp
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

const int TN = 60;                 // 测试用的最大点数

// ==================================================================
// 教程模板 1：点分治（统计距离 <= K 的点对数）
// ==================================================================
struct PEdge { int to, w, nxt; } pe[TN * 2];
int phead[TN], pecnt = 1;
void addPEdge(int u, int v, int w) {
    pe[++pecnt] = {v, w, phead[u]}; phead[u] = pecnt;
    pe[++pecnt] = {u, w, phead[v]}; phead[v] = pecnt;
}
bool premoved[TN];
int psz[TN];
int pK;

int getSize(int u, int fa) {
    psz[u] = 1;
    for (int i = phead[u]; i; i = pe[i].nxt) {
        int v = pe[i].to;
        if (v == fa || premoved[v]) continue;
        psz[u] += getSize(v, u);
    }
    return psz[u];
}
int getCentroid(int u, int fa, int tot) {
    for (int i = phead[u]; i; i = pe[i].nxt) {
        int v = pe[i].to;
        if (v == fa || premoved[v]) continue;
        if (psz[v] * 2 > tot) return getCentroid(v, u, tot);
    }
    return u;
}
int pdis[TN], pdisCnt;
void collectDis(int u, int fa, int d) {
    pdis[++pdisCnt] = d;
    for (int i = phead[u]; i; i = pe[i].nxt) {
        int v = pe[i].to;
        if (v == fa || premoved[v]) continue;
        collectDis(v, u, d + pe[i].w);
    }
}
ll calcPD(int u, int d0) {
    pdisCnt = 0;
    collectDis(u, 0, d0);
    sort(pdis + 1, pdis + pdisCnt + 1);
    ll res = 0;
    for (int l = 1, r = pdisCnt; l < r; ) {
        if (pdis[l] + pdis[r] <= pK) { res += r - l; ++l; }
        else --r;
    }
    return res;
}
ll pdAns = 0;
void dividePD(int u) {
    premoved[u] = true;
    pdAns += calcPD(u, 0);
    for (int i = phead[u]; i; i = pe[i].nxt) {
        int v = pe[i].to;
        if (premoved[v]) continue;
        pdAns -= calcPD(v, pe[i].w);
        int tot = getSize(v, u);
        int c = getCentroid(v, u, tot);
        dividePD(c);
    }
}

// ==================================================================
// 教程模板 2：DSU on tree（子树众数颜色编号之和）
// ==================================================================
const int DN = 60;
vector<int> dg[DN];
int dcol[DN];
int dsz[DN], heavy[DN];
int cntCol[DN];
int maxCnt;
ll curSum;
ll ansDsu[DN];

void dfsSizeDsu(int u, int fa) {
    dsz[u] = 1; heavy[u] = 0;
    for (int v : dg[u]) {
        if (v == fa) continue;
        dfsSizeDsu(v, u);
        dsz[u] += dsz[v];
        if (heavy[u] == 0 || dsz[v] > dsz[heavy[u]]) heavy[u] = v;
    }
}
void addSubtree(int u, int fa, int delta) {
    int c = dcol[u];
    cntCol[c] += delta;
    if (cntCol[c] > maxCnt) { maxCnt = cntCol[c]; curSum = c; }
    else if (cntCol[c] == maxCnt) { curSum += c; }
    for (int v : dg[u]) if (v != fa) addSubtree(v, u, delta);
}
void dfsDsu(int u, int fa, bool keep) {
    for (int v : dg[u])
        if (v != fa && v != heavy[u]) dfsDsu(v, u, false);
    if (heavy[u]) dfsDsu(heavy[u], u, true);
    for (int v : dg[u])
        if (v != fa && v != heavy[u]) addSubtree(v, u, +1);
    cntCol[dcol[u]]++;
    if (cntCol[dcol[u]] > maxCnt) { maxCnt = cntCol[dcol[u]]; curSum = dcol[u]; }
    else if (cntCol[dcol[u]] == maxCnt) { curSum += dcol[u]; }
    ansDsu[u] = curSum;
    if (!keep) {
        addSubtree(u, fa, -1);
        maxCnt = 0; curSum = 0;
    }
}

// ==================================================================
// 教程模板 3：长链剖分（CF1009F Dominant Indices）
// ==================================================================
const int LN = 65;
vector<int> lg2[LN];
int ldep[LN], lson[LN];
int lbuf[LN * 2];
int *lf[LN];
int lans[LN];
int lbufPtr;

void dfsDeep(int u, int fa) {
    for (int v : lg2[u]) {
        if (v == fa) continue;
        dfsDeep(v, u);
        if (ldep[v] > ldep[lson[u]]) lson[u] = v;
    }
    ldep[u] = ldep[lson[u]] + 1;
}
void solveLong(int u, int fa) {
    lf[u][0] = 1;
    if (lson[u]) {
        lf[lson[u]] = lf[u] + 1;
        solveLong(lson[u], u);
        lans[u] = lans[lson[u]] + 1;
    }
    for (int v : lg2[u]) {
        if (v == fa || v == lson[u]) continue;
        lf[v] = lbuf + lbufPtr; lbufPtr += ldep[v];
        solveLong(v, u);
        for (int d = 1; d <= ldep[v]; d++) {
            lf[u][d] += lf[v][d - 1];
            if (lf[u][d] > lf[u][lans[u]] ||
                (lf[u][d] == lf[u][lans[u]] && d < lans[u]))
                lans[u] = d;
        }
    }
    if (lf[u][lans[u]] == 1) lans[u] = 0;
}

// ==================================================================
// 教程模板 4：虚树
// ==================================================================
const int VN = 60;
vector<int> vg[VN];
int vdfn[VN], vdep[VN], vsiz[VN], vfa[VN][8], vtimer;

void dfsInitV(int u, int fa) {
    vdfn[u] = ++vtimer;
    vsiz[u] = 1;
    vfa[u][0] = fa;
    vdep[u] = vdep[fa] + 1;
    for (int k = 1; k < 8; k++) vfa[u][k] = vfa[vfa[u][k - 1]][k - 1];
    for (int v : vg[u]) {
        if (v == fa) continue;
        dfsInitV(v, u);
        vsiz[u] += vsiz[v];
    }
}
int lcaV(int a, int b) {
    if (vdep[a] < vdep[b]) swap(a, b);
    for (int k = 7; k >= 0; k--) if (vdep[a] - (1 << k) >= vdep[b]) a = vfa[a][k];
    if (a == b) return a;
    for (int k = 7; k >= 0; k--) if (vfa[a][k] != vfa[b][k]) { a = vfa[a][k]; b = vfa[b][k]; }
    return vfa[a][0];
}
vector<pair<int,int>> vEdges;
int vstk[VN], vtop;

void buildVirtualTree(vector<int> key) {
    sort(key.begin(), key.end(), [](int a, int b) { return vdfn[a] < vdfn[b]; });
    int m = (int)key.size();
    for (int i = 1; i < m; i++) key.push_back(lcaV(key[i - 1], key[i]));
    sort(key.begin(), key.end(), [](int a, int b) { return vdfn[a] < vdfn[b]; });
    key.erase(unique(key.begin(), key.end()), key.end());

    vEdges.clear();
    vtop = 0;
    for (int u : key) {
        if (vtop == 0) { vstk[++vtop] = u; continue; }
        int l = lcaV(u, vstk[vtop]);
        while (vtop >= 2 && vdep[vstk[vtop - 1]] >= vdep[l]) {
            vEdges.push_back({vstk[vtop - 1], vstk[vtop]});
            --vtop;
        }
        if (vstk[vtop] != l) {
            vEdges.push_back({l, vstk[vtop]});
            vstk[vtop] = l;
        }
        vstk[++vtop] = u;
    }
    while (vtop >= 2) {
        vEdges.push_back({vstk[vtop - 1], vstk[vtop]});
        --vtop;
    }
}

// ==================================================================
// 教程模板 5：树哈希（有根）
// ==================================================================
const int HN = 60;
vector<int> hg[HN];

unsigned long long splitmix(unsigned long long x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
unsigned long long dfsHash(int u, int fa) {
    unsigned long long h = 1;
    for (int v : hg[u]) {
        if (v == fa) continue;
        h += splitmix(dfsHash(v, u));
    }
    return h;
}

// ==================================================================
// 暴力参照
// ==================================================================
int tn;                            // 点数
int distM[TN][TN];

// 全源最短路（Floyd），用于点分治对拍
void buildDist(int n, const vector<array<int,3>> &edges) {
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++) distM[i][j] = (i == j ? 0 : (int)1e9);
    for (auto &ed : edges) {
        int u = ed[0], v = ed[1], w = ed[2];
        distM[u][v] = min(distM[u][v], w);
        distM[v][u] = min(distM[v][u], w);
    }
    for (int k = 1; k <= n; k++)
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                if (distM[i][k] + distM[k][j] < distM[i][j])
                    distM[i][j] = distM[i][k] + distM[k][j];
}

// ==================================================================
int main() {
    setvbuf(stdout, NULL, _IONBF, 0);

    printf("=== 1. 点分治 ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 200; iter++) {
            int n = (int)(rng() % 20) + 2;
            vector<array<int,3>> edges;
            for (int i = 2; i <= n; i++) {
                int p = (int)(rng() % (i - 1)) + 1;
                int w = (int)(rng() % 5) + 1;
                edges.push_back({p, i, w});
            }
            // 建点分治图
            pecnt = 1;
            for (int i = 1; i <= n; i++) { phead[i] = 0; premoved[i] = false; }
            for (auto &ed : edges) addPEdge(ed[0], ed[1], ed[2]);
            pK = (int)(rng() % 12);

            pdAns = 0;
            getSize(1, 0);
            int c = getCentroid(1, 0, getSize(1, 0));
            dividePD(c);

            // 暴力
            buildDist(n, edges);
            ll expect = 0;
            for (int i = 1; i <= n; i++)
                for (int j = i + 1; j <= n; j++)
                    if (distM[i][j] <= pK) ++expect;
            if (pdAns != expect) {
                ok = false;
                printf("    n=%d K=%d 得到 %lld 期望 %lld\n", n, pK, pdAns, expect);
                break;
            }
        }
        check(ok, "点分治统计距离 ≤ K 的点对数与 Floyd 暴力一致（200 组随机树）");
    }

    printf("=== 2. DSU on tree ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 200; iter++) {
            int n = (int)(rng() % 20) + 2;
            vector<pair<int,int>> treeEdges;
            for (int i = 1; i <= n; i++) dg[i].clear();
            for (int i = 2; i <= n; i++) {
                int p = (int)(rng() % (i - 1)) + 1;
                dg[p].push_back(i); dg[i].push_back(p);
                treeEdges.push_back({p, i});
            }
            for (int i = 1; i <= n; i++) dcol[i] = (int)(rng() % 6) + 1;

            // 以 1 为根定父子关系（与 dsu on tree 的子树定义保持一致）
            vector<int> par(n + 1, 0);
            {
                function<void(int,int)> setPar = [&](int x, int f) {
                    par[x] = f;
                    for (int v : dg[x]) if (v != f) setPar(v, x);
                };
                setPar(1, 0);
            }

            // 我的 dsu on tree（注意：cntCol 要整体清空，颜色值可能大于 n）
            for (int i = 0; i < DN; i++) cntCol[i] = 0;
            for (int i = 1; i <= n; i++) ansDsu[i] = 0;
            maxCnt = 0; curSum = 0;
            dfsSizeDsu(1, 0);
            dfsDsu(1, 0, false);

            // 暴力：对每个 u 遍历它的子树统计（注意要排除父节点，只走子树）
            for (int u = 1; u <= n; u++) {
                vector<int> sub;
                function<void(int,int)> collect = [&](int x, int fa) {
                    sub.push_back(x);
                    for (int v : dg[x]) if (v != fa) collect(v, x);
                };
                collect(u, par[u]);
                vector<int> cnt(10, 0);            // 下标是颜色值（1..6），不能按 n 开
                int best = 0;
                for (int x : sub) best = max(best, ++cnt[dcol[x]]);
                ll sum = 0;
                vector<int> seen(10, 0);
                for (int x : sub) {
                    int c = dcol[x];
                    if (cnt[c] == best && !seen[c]) { seen[c] = 1; sum += c; }
                }
                if (ansDsu[u] != sum) {
                    ok = false;
                    printf("    n=%d u=%d 得到 %lld 期望 %lld\n", n, u, ansDsu[u], sum);
                    break;
                }
            }
            if (!ok) break;
        }
        check(ok, "DSU on tree 子树众数统计与暴力一致（200 组随机树）");
    }

    printf("=== 3. 长链剖分（CF1009F） ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 150; iter++) {
            int n = (int)(rng() % 18) + 2;
            for (int i = 1; i <= n; i++) lg2[i].clear();
            for (int i = 2; i <= n; i++) {
                int p = (int)(rng() % (i - 1)) + 1;
                lg2[p].push_back(i); lg2[i].push_back(p);
            }
            for (int i = 1; i <= n; i++) { lson[i] = 0; lans[i] = 0; ldep[i] = 0; }
            fill(lbuf, lbuf + LN * 2, 0);          // 公用缓冲必须清空（跨轮次会残留）
            dfsDeep(1, 0);                         // 先求长儿子与长链长度
            lf[1] = lbuf;                          // 根的段从 buf 开始
            lbufPtr = ldep[1];                     // 关键：分配指针必须从"根段之后"开始
            solveLong(1, 0);

            // 以 1 为根定父子关系
            vector<int> lpar(n + 1, 0);
            {
                function<void(int,int)> setPar = [&](int x, int f) {
                    lpar[x] = f;
                    for (int v : lg2[x]) if (v != f) setPar(v, x);
                };
                setPar(1, 0);
            }

            // 暴力：对每个 u 统计它子树内每层的点数（从 u 出发，只走子树）
            for (int u = 1; u <= n; u++) {
                vector<int> cnt(n + 1, 0);
                function<void(int,int,int)> dfsCnt = [&](int x, int fa, int d) {
                    cnt[d]++;
                    for (int v : lg2[x]) if (v != fa) dfsCnt(v, x, d + 1);
                };
                dfsCnt(u, lpar[u], 0);
                int bestD = 0;
                for (int d = 1; d <= n; d++) if (cnt[d] > cnt[bestD]) bestD = d;
                if (lans[u] != bestD) {
                    ok = false;
                    printf("    n=%d u=%d 得到 %d 期望 %d\n", n, u, lans[u], bestD);
                    break;
                }
            }
            if (!ok) break;
        }
        check(ok, "长链剖分求 Dominant Indices 与暴力一致（150 组随机树）");
    }

    printf("=== 4. 虚树 ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 300; iter++) {
            int n = (int)(rng() % 25) + 3;
            for (int i = 1; i <= n; i++) vg[i].clear();
            for (int i = 2; i <= n; i++) {
                int p = (int)(rng() % (i - 1)) + 1;
                vg[p].push_back(i); vg[i].push_back(p);
            }
            vtimer = 0;
            for (int i = 0; i <= n; i++) {
                vdep[i] = 0; vsiz[i] = 0;
                for (int k = 0; k < 8; k++) vfa[i][k] = 0;
            }
            dfsInitV(1, 0);

            int k = (int)(rng() % n) + 1;
            vector<int> key;
            vector<char> isKey(n + 1, 0);
            for (int i = 0; i < k; i++) {
                int x = (int)(rng() % n) + 1;
                if (!isKey[x]) { isKey[x] = 1; key.push_back(x); }
            }
            buildVirtualTree(key);

            // 期望的点集 = 关键点 + 两两 LCA
            vector<char> expect(n + 1, 0);
            for (int x : key) expect[x] = 1;
            for (size_t i = 0; i < key.size(); i++)
                for (size_t j = i + 1; j < key.size(); j++)
                    expect[lcaV(key[i], key[j])] = 1;

            vector<char> got(n + 1, 0);
            for (auto &e2 : vEdges) { got[e2.first] = 1; got[e2.second] = 1; }
            for (int x : key) {
                // 单个关键点（无 LCA 边）时也要在点集中
                if (vEdges.empty()) got[x] = 1;
            }
            // 校验点集
            for (int i = 1; i <= n; i++) {
                if (expect[i] != got[i]) {
                    ok = false;
                    printf("    n=%d 点集不一致 i=%d expect=%d got=%d\n", n, i, (int)expect[i], (int)got[i]);
                    break;
                }
            }
            if (!ok) break;
            // 校验边数 = 点数 - 1（虚树是树）且每条边是祖先关系
            int vc = 0;
            for (int i = 1; i <= n; i++) vc += got[i];
            if ((int)vEdges.size() != vc - 1) {
                ok = false;
                printf("    n=%d 边数=%d 点数=%d\n", n, (int)vEdges.size(), vc);
                break;
            }
            for (auto &e2 : vEdges) {
                int a = e2.first, b = e2.second;
                if (!(vdfn[a] <= vdfn[b] && vdfn[b] < vdfn[a] + vsiz[a])) {
                    ok = false;
                    printf("    n=%d 虚树边 (%d,%d) 不是祖先关系\n", n, a, b);
                    break;
                }
            }
            if (!ok) break;
        }
        check(ok, "虚树点集 = 关键点 ∪ LCA，且构成合法的祖先树（300 组随机实例）");
    }

    printf("=== 5. 树哈希与同构判定 ===\n");
    {
        // 同构：随机树 + 随机置换，哈希应相同
        bool okSame = true;
        for (int iter = 0; iter < 300; iter++) {
            int n = (int)(rng() % 15) + 2;
            vector<pair<int,int>> edges;
            for (int i = 2; i <= n; i++) edges.push_back({(int)(rng() % (i - 1)) + 1, i});

            // 原树
            for (int i = 1; i <= n; i++) hg[i].clear();
            for (auto &e2 : edges) { hg[e2.first].push_back(e2.second); hg[e2.second].push_back(e2.first); }
            unsigned long long h1 = dfsHash(1, 0);

            // 随机置换标号
            vector<int> perm(n + 1);
            for (int i = 1; i <= n; i++) perm[i] = i;
            shuffle(perm.begin() + 1, perm.end(), rng);
            for (int i = 1; i <= n; i++) hg[i].clear();
            for (auto &e2 : edges) {
                int u = perm[e2.first], v = perm[e2.second];
                hg[u].push_back(v); hg[v].push_back(u);
            }
            unsigned long long h2 = dfsHash(perm[1], 0);       // 以对应节点为根
            if (h1 != h2) {
                okSame = false;
                printf("    n=%d 同构树哈希不同\n", n);
                break;
            }
        }
        check(okSame, "同构的有根树（随机重标号）哈希相同（300 组）");

        // 严谨性验证：用 AHU 规范化括号串做"真实同构判定"（与哈希无关）
        {
            auto ahuStr = [&](auto &&self, int u, int fa) -> string {
                vector<string> subs;
                for (int v : hg[u]) if (v != fa) subs.push_back(self(self, v, u));
                sort(subs.begin(), subs.end());
                string s = "(";
                for (auto &t : subs) s += t;
                s += ")";
                return s;
            };
            vector<unsigned long long> hs;
            vector<string> ss;
            for (int iter = 0; iter < 300; iter++) {
                int n = (int)(rng() % 10) + 2;
                for (int i = 1; i <= n; i++) hg[i].clear();
                for (int i = 2; i <= n; i++) {
                    int p = (int)(rng() % (i - 1)) + 1;
                    hg[p].push_back(i); hg[i].push_back(p);
                }
                hs.push_back(dfsHash(1, 0));
                ss.push_back(ahuStr(ahuStr, 1, 0));
            }
            bool okImp = true;
            int collide = 0;
            for (size_t i = 0; i < hs.size(); i++)
                for (size_t j = i + 1; j < hs.size(); j++) {
                    bool sameShape = (ss[i] == ss[j]);
                    bool sameHash = (hs[i] == hs[j]);
                    if (sameShape && !sameHash) okImp = false;      // 同构必须同哈希
                    if (!sameShape && sameHash) ++collide;          // 不同构却同哈希 = 碰撞
                }
            check(okImp, "树哈希满足「同构 ⟹ 同哈希」（300 棵随机树两两比较）");
            long long pairs = (long long)hs.size() * (hs.size() - 1) / 2;
            printf("    （不同构但哈希相同的对数：%d / %lld）\n", collide, pairs);
        }
    }

    printf("\n===== 通过 %d 项，失败 %d 项 =====\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
