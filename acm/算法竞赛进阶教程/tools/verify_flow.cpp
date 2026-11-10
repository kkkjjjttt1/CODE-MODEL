// 第 02 章 模板验证：Dinic 最大流 / 最小割 / 费用流 / 二分图匹配 / KM / 上下界可行流
// 编译: g++ -std=c++17 -O2 -o verify_flow.exe tools/verify_flow.cpp
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

const int MAXN = 60, MAXM = 400;
const ll INF = (ll)4e18;

// ==================================================================
// 教程模板：Dinic
// ==================================================================
struct Edge { int to, nxt; ll cap; } e[MAXM * 2];
int head[MAXN], ecnt = 1;
int n, S, T;

void init_flow(int nn) {
    n = nn; ecnt = 1;
    for (int i = 1; i <= n + 2; i++) head[i] = 0;
}
void addEdge(int u, int v, ll cap) {
    e[++ecnt] = {v, head[u], cap}; head[u] = ecnt;
    e[++ecnt] = {u, head[v], 0};   head[v] = ecnt;
}
int level[MAXN], cur[MAXN];

bool bfs() {
    for (int i = 1; i <= n; i++) level[i] = -1;
    static int q[MAXN];
    int qh = 0, qt = 0;
    level[S] = 0; q[qt++] = S;
    while (qh < qt) {
        int u = q[qh++];
        for (int i = head[u]; i; i = e[i].nxt) {
            int v = e[i].to;
            if (e[i].cap > 0 && level[v] < 0) { level[v] = level[u] + 1; q[qt++] = v; }
        }
    }
    return level[T] >= 0;
}
ll dfs(int u, ll f) {
    if (u == T) return f;
    ll used = 0;
    for (int &i = cur[u]; i; i = e[i].nxt) {
        int v = e[i].to;
        if (e[i].cap > 0 && level[v] == level[u] + 1) {
            ll d = dfs(v, min(f - used, e[i].cap));
            if (d > 0) {
                e[i].cap -= d;
                e[i ^ 1].cap += d;
                used += d;
                if (used == f) break;
            }
        }
    }
    return used;
}
ll dinic() {
    ll flow = 0;
    while (bfs()) {
        for (int i = 1; i <= n; i++) cur[i] = head[i];
        ll f;
        while ((f = dfs(S, INF)) > 0) flow += f;
    }
    return flow;
}
vector<char> minCutSide() {
    vector<char> vis(n + 1, 0);
    static int q[MAXN];
    int qh = 0, qt = 0;
    vis[S] = 1; q[qt++] = S;
    while (qh < qt) {
        int u = q[qh++];
        for (int i = head[u]; i; i = e[i].nxt)
            if (e[i].cap > 0 && !vis[e[i].to]) { vis[e[i].to] = 1; q[qt++] = e[i].to; }
    }
    return vis;
}

// ==================================================================
// 教程模板：费用流（SPFA 版）
// ==================================================================
struct Edge2 { int to, nxt; ll cap, cost; } e2[MAXM * 2];
int head2[MAXN], ecnt2 = 1;
int n2, S2, T2;
const ll INF2 = (ll)4e18;

void init_mcmf(int nn) {
    n2 = nn; ecnt2 = 1;
    for (int i = 1; i <= nn; i++) head2[i] = 0;
}
void addEdge2(int u, int v, ll cap, ll cost) {
    e2[++ecnt2] = {v, head2[u], cap, cost};  head2[u] = ecnt2;
    e2[++ecnt2] = {u, head2[v], 0, -cost};   head2[v] = ecnt2;
}
ll dis2[MAXN];
int preE[MAXN];
bool inq[MAXN];

bool spfa_mcmf() {
    for (int i = 1; i <= n2; i++) { dis2[i] = INF2; inq[i] = false; preE[i] = 0; }
    queue<int> q;
    dis2[S2] = 0; q.push(S2); inq[S2] = true;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        inq[u] = false;
        for (int i = head2[u]; i; i = e2[i].nxt) {
            int v = e2[i].to;
            if (e2[i].cap > 0 && dis2[u] + e2[i].cost < dis2[v]) {
                dis2[v] = dis2[u] + e2[i].cost;
                preE[v] = i;
                if (!inq[v]) { inq[v] = true; q.push(v); }
            }
        }
    }
    return dis2[T2] < INF2;
}
pair<ll, ll> mcmf() {
    ll flow = 0, cost = 0;
    while (spfa_mcmf()) {
        ll f = INF2;
        for (int v = T2; v != S2; v = e2[preE[v] ^ 1].to) f = min(f, e2[preE[v]].cap);
        for (int v = T2; v != S2; v = e2[preE[v] ^ 1].to) {
            e2[preE[v]].cap -= f;
            e2[preE[v] ^ 1].cap += f;
        }
        flow += f;
        cost += f * dis2[T2];
    }
    return {flow, cost};
}

// ==================================================================
// 教程模板：匈牙利
// ==================================================================
const int MAXL = 25;
vector<int> adjL[MAXL];
int matchR[MAXL];

bool tryKuhn(int u, vector<int> &vis, int tag) {
    for (int v : adjL[u]) {
        if (vis[v] == tag) continue;
        vis[v] = tag;
        if (!matchR[v] || tryKuhn(matchR[v], vis, tag)) { matchR[v] = u; return true; }
    }
    return false;
}
int hungarian(int nL) {
    fill(matchR, matchR + MAXL, 0);
    vector<int> vis(MAXL, 0);
    int res = 0;
    for (int u = 1; u <= nL; u++) if (tryKuhn(u, vis, u)) ++res;
    return res;
}

// ==================================================================
// 教程模板：KM（最大权完美匹配）
// ==================================================================
const int MAXK = 12;
typedef ll llw;
const llw KINF = (llw)4e18;
int kn;
llw kw[MAXK][MAXK];
llw lx[MAXK], ly[MAXK];
int matchY[MAXK], matchU[MAXK];      // matchY[v]=右部 v 匹配的左部点；matchU[u]=左部 u 匹配的右部点

// BFS 版 KM（O(n^3)）：每次增广用 BFS 扩展交错树，避免递归版"顶标调整不收敛"的坑
llw km(int n_) {
    kn = n_;
    for (int i = 1; i <= kn; i++) {
        lx[i] = 0; ly[i] = 0; matchY[i] = 0; matchU[i] = 0;
        for (int j = 1; j <= kn; j++) lx[i] = max(lx[i], kw[i][j]);
    }
    for (int u = 1; u <= kn; u++) {
        vector<char> visX(kn + 1, 0), visY(kn + 1, 0);
        vector<llw> slack(kn + 1, KINF);
        vector<int> pre(kn + 1, 0);
        int cur = u;
        visX[cur] = 1;
        while (true) {
            for (int v = 1; v <= kn; v++) {
                if (visY[v]) continue;
                llw gap = lx[cur] + ly[v] - kw[cur][v];
                if (gap < slack[v]) { slack[v] = gap; pre[v] = cur; }
            }
            llw d = KINF; int nxt = -1;
            for (int v = 1; v <= kn; v++)
                if (!visY[v] && slack[v] < d) { d = slack[v]; nxt = v; }
            if (nxt < 0) return -KINF;                 // 不存在完美匹配
            for (int i = 1; i <= kn; i++) {
                if (visX[i]) lx[i] -= d;
                if (visY[i]) ly[i] += d;
                else slack[i] -= d;
            }
            if (matchY[nxt] == 0) {                    // 增广
                int v = nxt;
                while (v) {
                    int u2 = pre[v];
                    int nv = matchU[u2];
                    matchY[v] = u2;
                    matchU[u2] = v;
                    v = nv;
                }
                break;
            } else {
                visY[nxt] = 1;
                cur = matchY[nxt];
                visX[cur] = 1;
            }
        }
    }
    llw res = 0;
    for (int v = 1; v <= kn; v++) if (matchY[v]) res += kw[matchY[v]][v];
    return res;
}

// ==================================================================
// 暴力参照
// ==================================================================
// 暴力最大流/最小割：枚举所有 s-t 割（小图）
ll bruteMinCut(int nn, int s, int t, const vector<array<ll,3>> &edges) {
    ll best = INF;
    for (int mask = 0; mask < (1 << nn); mask++) {
        if (!((mask >> (s - 1)) & 1)) continue;       // s 必须在 S 侧
        if ((mask >> (t - 1)) & 1) continue;          // t 必须在 T 侧
        ll cut = 0;
        for (auto &ed : edges) {
            int u = (int)ed[0], v = (int)ed[1]; ll c = ed[2];
            bool su = (mask >> (u - 1)) & 1, sv = (mask >> (v - 1)) & 1;
            if (su && !sv) cut += c;
        }
        best = min(best, cut);
    }
    return best;
}

// 暴力最小费用最大流：枚举每条边的流量（小容量小图）
pair<ll,ll> bruteMCMF(int nn, int s, int t, const vector<array<ll,4>> &edges) {
    // edges: {u, v, cap, cost}
    int m = (int)edges.size();
    ll bestCost = INF, bestFlow = -1;
    vector<int> f(m, 0);
    function<void(int)> rec = [&](int idx) {
        if (idx == m) {
            // 检查流量守恒 + 计算流量与费用
            vector<ll> bal(nn + 1, 0);
            ll cost = 0;
            for (int i = 0; i < m; i++) {
                bal[(int)edges[i][0]] -= f[i];
                bal[(int)edges[i][1]] += f[i];
                cost += (ll)f[i] * edges[i][3];
            }
            for (int v = 1; v <= nn; v++) {
                if (v == s || v == t) continue;
                if (bal[v] != 0) return;
            }
            ll flow = bal[t];                       // t 的净流入 = 流量值
            if (flow != -bal[s]) return;
            if (flow < 0) return;
            if (flow > bestFlow) { bestFlow = flow; bestCost = cost; }
            else if (flow == bestFlow && cost < bestCost) bestCost = cost;
            return;
        }
        for (int val = 0; val <= (int)edges[idx][2]; val++) {
            f[idx] = val;
            rec(idx + 1);
        }
        f[idx] = 0;
    };
    rec(0);
    if (bestFlow < 0) return {0, 0};
    return {bestFlow, bestCost};
}

// 暴力二分图最大匹配：状压 DP（右部点用位掩码）
int bruteMatching(int nL, int nR, const vector<pair<int,int>> &edges) {
    vector<int> adjMask(nL + 1, 0);
    for (auto &pr : edges) adjMask[pr.first] |= (1 << (pr.second - 1));
    vector<int> dp(1 << nR, -1);
    dp[0] = 0;
    for (int u = 1; u <= nL; u++) {
        vector<int> ndp = dp;
        for (int mask = 0; mask < (1 << nR); mask++) {
            if (dp[mask] < 0) continue;
            for (int v = 1; v <= nR; v++) {
                if (!((adjMask[u] >> (v - 1)) & 1)) continue;
                if ((mask >> (v - 1)) & 1) continue;
                int nm = mask | (1 << (v - 1));
                ndp[nm] = max(ndp[nm], dp[mask] + 1);
            }
        }
        dp = ndp;
    }
    int best = 0;
    for (int mask = 0; mask < (1 << nR); mask++) best = max(best, dp[mask]);
    return best;
}

// 暴力 KM：枚举全排列
llw bruteKM(int nn, llw w[][MAXK]) {
    vector<int> perm(nn);
    for (int i = 0; i < nn; i++) perm[i] = i + 1;
    llw best = -KINF;
    do {
        llw s = 0;
        for (int i = 1; i <= nn; i++) s += w[i][perm[i - 1]];
        best = max(best, s);
    } while (next_permutation(perm.begin(), perm.end()));
    return best;
}

// ==================================================================
int main() {
    setvbuf(stdout, NULL, _IONBF, 0);      // 便于定位卡点
    printf("=== 1. Dinic 最大流 / 最小割 ===\n");
    {
        bool okFlow = true;
        for (int iter = 0; iter < 500; iter++) {
            int nn = (int)(rng() % 6) + 2;              // 点 1..nn，S=1，T=nn
            int m = (int)(rng() % 9);
            vector<array<ll,3>> edges;
            for (int i = 0; i < m; i++) {
                int u = (int)(rng() % nn) + 1, v = (int)(rng() % nn) + 1;
                if (u == v) continue;
                ll c = (ll)(rng() % 9) + 1;
                edges.push_back({(ll)u, (ll)v, c});
            }
            S = 1; T = nn;
            init_flow(nn);
            for (auto &ed : edges) addEdge((int)ed[0], (int)ed[1], ed[2]);
            ll flow = dinic();
            ll cut = bruteMinCut(nn, S, T, edges);
            if (flow != cut) {
                okFlow = false;
                printf("    nn=%d 最大流=%lld 最小割=%lld\n", nn, flow, cut);
                break;
            }
            // 最小割 S 侧：验证割容量 = flow，且 S 侧包含 s、不含 t
            vector<char> side = minCutSide();
            if (!side[S] || side[T]) { okFlow = false; printf("    S 侧划分错误\n"); break; }
            ll cutCap = 0;
            for (auto &ed : edges)
                if (side[(int)ed[0]] && !side[(int)ed[1]]) cutCap += ed[2];
            if (cutCap != flow) { okFlow = false; printf("    S 侧割容量=%lld 流量=%lld\n", cutCap, flow); break; }
        }
        check(okFlow, "Dinic 最大流 == 暴力枚举最小割（500 组随机图，含 S 侧划分验证）");
    }

    printf("=== 2. 费用流（SPFA 版） ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 300; iter++) {
            int nn = (int)(rng() % 4) + 2;              // 小图，最多 4 个点
            int m = (int)(rng() % 4) + 1;
            vector<array<ll,4>> edges;
            for (int i = 0; i < m; i++) {
                int u = (int)(rng() % nn) + 1, v = (int)(rng() % nn) + 1;
                if (u == v) continue;
                ll c = (ll)(rng() % 3);                 // 容量 0..2
                ll w = (ll)(rng() % 7);                 // 费用 0..6（非负）
                if (c == 0) continue;
                edges.push_back({(ll)u, (ll)v, c, w});
            }
            if (edges.empty()) continue;
            S2 = 1; T2 = nn;
            init_mcmf(nn);
            for (auto &ed : edges) addEdge2((int)ed[0], (int)ed[1], ed[2], ed[3]);
            pair<ll,ll> got = mcmf();
            pair<ll,ll> expect = bruteMCMF(nn, S2, T2, edges);
            if (got.first != expect.first || got.second != expect.second) {
                ok = false;
                printf("    nn=%d m=%d 得到(%lld,%lld) 暴力(%lld,%lld)\n",
                       nn, (int)edges.size(), got.first, got.second, expect.first, expect.second);
                break;
            }
        }
        check(ok, "SPFA 费用流与暴力枚举流量一致（最小费用最大流）");
    }

    printf("=== 3. 二分图最大匹配 ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 600; iter++) {
            int nL = (int)(rng() % 5) + 1, nR = (int)(rng() % 5) + 1;
            vector<pair<int,int>> edges;
            for (int u = 1; u <= nL; u++)
                for (int v = 1; v <= nR; v++)
                    if (rng() % 3 == 0) edges.push_back({u, v});
            for (int u = 1; u <= nL; u++) adjL[u].clear();
            for (auto &pr : edges) adjL[pr.first].push_back(pr.second);
            int got = hungarian(nL);
            int expect = bruteMatching(nL, nR, edges);
            if (got != expect) {
                ok = false;
                printf("    nL=%d nR=%d 得到 %d 期望 %d\n", nL, nR, got, expect);
                break;
            }
            // 同时用 Dinic 版验证
            S = nL + nR + 1; T = nL + nR + 2;
            init_flow(nL + nR + 2);
            for (int u = 1; u <= nL; u++) addEdge(S, u, 1);
            for (int v = 1; v <= nR; v++) addEdge(nL + v, T, 1);
            for (auto &pr : edges) addEdge(pr.first, nL + pr.second, 1);
            ll dinicMatch = dinic();
            if (dinicMatch != expect) {
                ok = false;
                printf("    Dinic 匹配 %lld 期望 %d\n", dinicMatch, expect);
                break;
            }
        }
        check(ok, "匈牙利 & Dinic 二分图匹配与状压 DP 暴力一致");
    }

    printf("=== 4. KM 带权完美匹配 ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 200; iter++) {
            int nn = (int)(rng() % 6) + 1;
            for (int i = 1; i <= nn; i++)
                for (int j = 1; j <= nn; j++)
                    kw[i][j] = (llw)(rng() % 21);       // 0..20
            llw got = km(nn);
            llw expect = bruteKM(nn, kw);
            if (got != expect) {
                ok = false;
                printf("    nn=%d 得到 %lld 期望 %lld\n", nn, got, expect);
                break;
            }
        }
        check(ok, "KM 最大权完美匹配与全排列暴力一致");
    }

    printf("=== 5. 上下界可行流（小规模暴力） ===\n");
    {
        // 暴力：枚举每条边的流量（lo..hi），检查流量守恒
        bool ok = true;
        for (int iter = 0; iter < 400; iter++) {
            int nn = (int)(rng() % 4) + 2;
            int m = (int)(rng() % 4) + 1;
            vector<array<ll,4>> bes;                    // {u, v, lo, hi}
            for (int i = 0; i < m; i++) {
                int u = (int)(rng() % nn) + 1, v = (int)(rng() % nn) + 1;
                if (u == v) continue;
                ll lo = (ll)(rng() % 2), hi = lo + (ll)(rng() % 2);
                bes.push_back({(ll)u, (ll)v, lo, hi});
            }
            if (bes.empty()) continue;
            // 暴力判定
            int mm = (int)bes.size();
            vector<int> f(mm, 0);
            bool bruteFeasible = false;
            function<void(int)> rec = [&](int idx) {
                if (bruteFeasible) return;
                if (idx == mm) {
                    vector<ll> bal(nn + 1, 0);
                    for (int i = 0; i < mm; i++) {
                        bal[(int)bes[i][0]] -= f[i];
                        bal[(int)bes[i][1]] += f[i];
                    }
                    for (int v = 1; v <= nn; v++) if (bal[v] != 0) return;
                    bruteFeasible = true;
                    return;
                }
                for (int val = (int)bes[idx][2]; val <= (int)bes[idx][3]; val++) {
                    f[idx] = val;
                    rec(idx + 1);
                }
            };
            rec(0);

            // 我的算法：无源汇上下界可行流
            int SS = nn + 1, TT = nn + 2;
            init_flow(nn + 2);
            vector<ll> d(nn + 2, 0);
            vector<int> pos(mm, 0);
            for (int i = 0; i < mm; i++) {
                pos[i] = ecnt + 1;
                addEdge((int)bes[i][0], (int)bes[i][1], bes[i][3] - bes[i][2]);
                d[(int)bes[i][0]] -= bes[i][2];
                d[(int)bes[i][1]] += bes[i][2];
            }
            ll need = 0;
            for (int v = 1; v <= nn; v++) {
                if (d[v] > 0) { addEdge(SS, v, d[v]); need += d[v]; }
                else if (d[v] < 0) addEdge(v, TT, -d[v]);
            }
            S = SS; T = TT;
            ll got = dinic();
            bool myFeasible = (got == need);

            if (myFeasible != bruteFeasible) {
                ok = false;
                printf("    nn=%d m=%d 算法=%d 暴力=%d\n", nn, mm, (int)myFeasible, (int)bruteFeasible);
                break;
            }
            // 若可行，验证还原出的流量满足上下界与守恒
            if (myFeasible) {
                vector<ll> flow(mm), bal(nn + 2, 0);
                for (int i = 0; i < mm; i++) {
                    flow[i] = bes[i][2] + (bes[i][3] - bes[i][2] - e[pos[i]].cap);
                    if (flow[i] < bes[i][2] || flow[i] > bes[i][3]) { ok = false; break; }
                    bal[(int)bes[i][0]] -= flow[i];
                    bal[(int)bes[i][1]] += flow[i];
                }
                for (int v = 1; v <= nn && ok; v++) if (bal[v] != 0) ok = false;
            }
            if (!ok) break;
        }
        check(ok, "无源汇上下界可行流：存在性与还原流量均与暴力一致");
    }

    printf("\n===== 通过 %d 项，失败 %d 项 =====\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
