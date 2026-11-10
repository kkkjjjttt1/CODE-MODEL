// 第 05/06/07 章 模板验证：CDQ / 整体二分 / FWT / 子集卷积 / 类欧 / min_25 / DLX / 双向BFS
// 编译: g++ -std=c++17 -O2 -o verify_adv.exe tools/verify_adv.cpp
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

// ==================================================================
// 1. CDQ 分治：三维偏序
// ==================================================================
const int CDQN = 200;
struct Node { int a, b, c, cnt, ans; } cdqP[CDQN];
int cdqValMax, cdqBit[CDQN];

void cdqBitAdd(int i, int v) { for (; i <= cdqValMax; i += i & -i) cdqBit[i] += v; }
int cdqBitSum(int i) { int s = 0; for (; i > 0; i -= i & -i) s += cdqBit[i]; return s; }

bool cmpCdqA(const Node &x, const Node &y) {
    if (x.a != y.a) return x.a < y.a;
    if (x.b != y.b) return x.b < y.b;
    return x.c < y.c;
}
bool cmpCdqB(const Node &x, const Node &y) {
    if (x.b != y.b) return x.b < y.b;
    return x.c < y.c;
}

void cdq(int l, int r) {
    if (l == r) return;
    int mid = (l + r) >> 1;
    cdq(l, mid);
    cdq(mid + 1, r);
    sort(cdqP + l, cdqP + mid + 1, cmpCdqB);
    sort(cdqP + mid + 1, cdqP + r + 1, cmpCdqB);
    int i = l;
    for (int j = mid + 1; j <= r; j++) {
        while (i <= mid && cdqP[i].b <= cdqP[j].b) { cdqBitAdd(cdqP[i].c, cdqP[i].cnt); ++i; }
        cdqP[j].ans += cdqBitSum(cdqP[j].c);
    }
    for (int t = l; t < i; t++) cdqBitAdd(cdqP[t].c, -cdqP[t].cnt);
    inplace_merge(cdqP + l, cdqP + mid + 1, cdqP + r + 1, cmpCdqB);
}

// ==================================================================
// 2. 整体二分：静态区间第 k 小
// ==================================================================
const int OAN = 300;
struct OQuery { int l, r, k, id; };
OQuery oQs[OAN], oQl[OAN], oQr[OAN];
int oA[OAN], oAns[OAN], oBit[OAN], oN, oM;

void oBitAdd(int i, int v) { for (; i <= oN; i += i & -i) oBit[i] += v; }
int oBitSum(int i) { int s = 0; for (; i > 0; i -= i & -i) s += oBit[i]; return s; }

void solveAll(vector<OQuery> qs, int vl, int vr) {
    if (qs.empty()) return;
    if (vl == vr) {
        for (const auto &q : qs) oAns[q.id] = vl;
        return;
    }
    int mid = (vl + vr) >> 1;
    // 关键：只加入「值落在当前值域左半 [vl, mid]」的元素
    for (int i = 1; i <= oN; i++) if (vl <= oA[i] && oA[i] <= mid) oBitAdd(i, 1);
    vector<OQuery> leftQs, rightQs;
    for (auto q : qs) {
        int c = oBitSum(q.r) - oBitSum(q.l - 1);
        if (c >= q.k) leftQs.push_back(q);
        else { q.k -= c; rightQs.push_back(q); }
    }
    for (int i = 1; i <= oN; i++) if (vl <= oA[i] && oA[i] <= mid) oBitAdd(i, -1);
    solveAll(leftQs, vl, mid);
    solveAll(rightQs, mid + 1, vr);
}

// ==================================================================
// 3. FWT
// ==================================================================
void fwtOr(vector<ll> &a, bool inv) {
    int n = (int)a.size();
    for (int len = 1; len < n; len <<= 1)
        for (int i = 0; i < n; i += (len << 1))
            for (int j = i; j < i + len; j++)
                a[j + len] += inv ? -a[j] : a[j];
}
void fwtAnd(vector<ll> &a, bool inv) {
    int n = (int)a.size();
    for (int len = 1; len < n; len <<= 1)
        for (int i = 0; i < n; i += (len << 1))
            for (int j = i; j < i + len; j++)
                a[j] += inv ? -a[j + len] : a[j + len];
}
void fwtXor(vector<ll> &a, bool inv) {
    int n = (int)a.size();
    for (int len = 1; len < n; len <<= 1)
        for (int i = 0; i < n; i += (len << 1))
            for (int j = i; j < i + len; j++) {
                ll u = a[j], v = a[j + len];
                a[j] = u + v; a[j + len] = u - v;
            }
    if (inv) for (auto &x : a) x /= n;
}
// type: 0=OR, 1=AND, 2=XOR
vector<ll> walshConv(vector<ll> a, vector<ll> b, int type) {
    int n = (int)a.size();
    if (type == 0) { fwtOr(a, false); fwtOr(b, false); }
    else if (type == 1) { fwtAnd(a, false); fwtAnd(b, false); }
    else { fwtXor(a, false); fwtXor(b, false); }
    for (int i = 0; i < n; i++) a[i] *= b[i];
    if (type == 0) fwtOr(a, true);
    else if (type == 1) fwtAnd(a, true);
    else fwtXor(a, true);
    return a;
}

// 子集卷积
vector<ll> subsetConvolution(const vector<ll> &a, const vector<ll> &b) {
    int n = (int)a.size();
    int k = 0;
    while ((1 << k) < n) ++k;
    vector<vector<ll>> A(k + 1, vector<ll>(n, 0)), B(k + 1, vector<ll>(n, 0));
    for (int S = 0; S < n; S++) {
        int pc = __builtin_popcount((unsigned)S);
        A[pc][S] = a[S];
        B[pc][S] = b[S];
    }
    for (int i = 0; i <= k; i++) { fwtOr(A[i], false); fwtOr(B[i], false); }
    vector<vector<ll>> C(k + 1, vector<ll>(n, 0));
    for (int i = 0; i <= k; i++)
        for (int j = 0; i + j <= k; j++)
            for (int S = 0; S < n; S++)
                C[i + j][S] += A[i][S] * B[j][S];
    for (int i = 0; i <= k; i++) fwtOr(C[i], true);
    vector<ll> c(n, 0);
    for (int S = 0; S < n; S++) c[S] = C[__builtin_popcount((unsigned)S)][S];
    return c;
}

// ==================================================================
// 4. 类欧几里得
// ==================================================================
ll floorSum(ll n, ll m, ll a, ll b) {
    ll ans = 0;
    while (true) {
        if (a >= m) { ans += (n - 1) * n / 2 * (a / m); a %= m; }
        if (b >= m) { ans += n * (b / m); b %= m; }
        ll yMax = a * n + b;
        if (yMax < m) break;
        n = yMax / m;
        b = yMax % m;
        swap(m, a);
    }
    return ans;
}

// ==================================================================
// 5. min_25 筛（f(p^k) = p^k(p^k-1)）
// ==================================================================
const int MN25 = 100005;
const ll MOD25 = 1000000007, INV2 = 500000004, INV6 = 166666668;
ll mn25N, mn25Sq;
ll mn25W[MN25 * 2];
int mn25Id1[MN25], mn25Id2[MN25], mn25M;
ll mn25G1[MN25 * 2], mn25G2[MN25 * 2];
ll mn25Pre1[MN25], mn25Pre2[MN25];
int mn25Primes[MN25], mn25Pcnt;
bool mn25Vis[MN25];

inline ll sumPow1(ll x) { return x % MOD25 * ((x + 1) % MOD25) % MOD25 * INV2 % MOD25; }
inline ll sumPow2(ll x) { return x % MOD25 * ((x + 1) % MOD25) % MOD25 * ((2 * x + 1) % MOD25) % MOD25 * INV6 % MOD25; }
inline ll fSingle(ll pk) { ll t = pk % MOD25; return t * ((t - 1 + MOD25) % MOD25) % MOD25; }
inline int getId25(ll x) { return x <= mn25Sq ? mn25Id1[x] : mn25Id2[mn25N / x]; }

void initMin25(ll n) {
    mn25N = n;
    mn25Sq = (ll)sqrt((double)n);
    while (mn25Sq * mn25Sq > n) --mn25Sq;
    while ((mn25Sq + 1) * (mn25Sq + 1) <= n) ++mn25Sq;

    mn25Pcnt = 0;
    mn25Pre1[0] = mn25Pre2[0] = 0;
    for (int i = 0; i <= mn25Sq; i++) mn25Vis[i] = false;
    for (int i = 2; i <= mn25Sq; i++) {
        if (!mn25Vis[i]) {
            mn25Primes[++mn25Pcnt] = i;
            mn25Pre1[mn25Pcnt] = (mn25Pre1[mn25Pcnt - 1] + i) % MOD25;
            mn25Pre2[mn25Pcnt] = (mn25Pre2[mn25Pcnt - 1] + (ll)i * i) % MOD25;
        }
        for (int j = 1; j <= mn25Pcnt && (ll)i * mn25Primes[j] <= mn25Sq; j++) {
            mn25Vis[i * mn25Primes[j]] = true;
            if (i % mn25Primes[j] == 0) break;
        }
    }

    mn25M = 0;
    for (ll l = 1, r; l <= n; l = r + 1) {
        r = n / (n / l);
        mn25W[++mn25M] = n / l;
        if (mn25W[mn25M] <= mn25Sq) mn25Id1[mn25W[mn25M]] = mn25M;
        else mn25Id2[n / mn25W[mn25M]] = mn25M;
        mn25G1[mn25M] = (sumPow1(mn25W[mn25M]) - 1 + MOD25) % MOD25;
        mn25G2[mn25M] = (sumPow2(mn25W[mn25M]) - 1 + MOD25) % MOD25;
    }
    for (int j = 1; j <= mn25Pcnt; j++) {
        ll p = mn25Primes[j];
        for (int i = 1; i <= mn25M && p * p <= mn25W[i]; i++) {
            int k = getId25(mn25W[i] / p);
            mn25G1[i] = (mn25G1[i] - p % MOD25 * ((mn25G1[k] - mn25Pre1[j - 1] + MOD25) % MOD25) % MOD25 + MOD25) % MOD25;
            mn25G2[i] = (mn25G2[i] - p % MOD25 * (p % MOD25) % MOD25 * ((mn25G2[k] - mn25Pre2[j - 1] + MOD25) % MOD25) % MOD25 + MOD25) % MOD25;
        }
    }
}

ll S25(ll x, int y) {
    if (mn25Primes[y] >= x) return 0;
    int k = getId25(x);
    ll ans = ((mn25G2[k] - mn25Pre2[y] + MOD25) % MOD25 - (mn25G1[k] - mn25Pre1[y] + MOD25) % MOD25 + MOD25) % MOD25;
    for (int i = y + 1; i <= mn25Pcnt && (ll)mn25Primes[i] * mn25Primes[i] <= x; i++) {
        ll pe = mn25Primes[i];
        for (int e = 1; pe * mn25Primes[i] <= x; e++, pe *= mn25Primes[i]) {
            ans = (ans + fSingle(pe) * S25(x / pe, i) + fSingle(pe * mn25Primes[i])) % MOD25;
        }
    }
    return ans;
}

// 暴力：线性筛求同一个 f 的前缀和
vector<ll> brutePrefixMin25(int n) {
    vector<ll> pkVal(n + 1, 1), fval(n + 1, 0), pre(n + 1, 0);
    vector<int> rest(n + 1, 1), primes;
    vector<bool> comp(n + 1, false);
    fval[1] = 1;
    for (int i = 2; i <= n; i++) {
        if (!comp[i]) { primes.push_back(i); pkVal[i] = i; rest[i] = 1; }
        for (int j = 0; j < (int)primes.size() && (ll)i * primes[j] <= n; j++) {
            int p = primes[j], v = i * p;
            comp[v] = true;
            if (i % p == 0) { pkVal[v] = pkVal[i] * p; rest[v] = rest[i]; break; }
            pkVal[v] = p; rest[v] = i;
        }
    }
    for (int i = 1; i <= n; i++) {
        if (i == 1) fval[i] = 1;
        else {
            ll t = pkVal[i] % MOD25;
            fval[i] = t * ((t - 1 + MOD25) % MOD25) % MOD25 * fval[rest[i]] % MOD25;
        }
        pre[i] = (pre[i - 1] + fval[i]) % MOD25;
    }
    return pre;
}

// ==================================================================
// 6. DLX（精确覆盖）
// ==================================================================
const int DLXNODE = 20005;
int dlxL[DLXNODE], dlxR[DLXNODE], dlxU[DLXNODE], dlxD[DLXNODE];
int dlxRow[DLXNODE], dlxCol[DLXNODE], dlxHead[200], dlxSize[200];
int dlxTot, dlxColCnt;
vector<int> dlxAns;
bool dlxFound;

void dlxInit(int colCnt) {
    dlxColCnt = colCnt;
    dlxTot = colCnt;
    for (int i = 0; i <= colCnt; i++) {
        dlxL[i] = i - 1; dlxR[i] = i + 1;
        dlxU[i] = dlxD[i] = i;
        dlxSize[i] = 0; dlxRow[i] = 0; dlxCol[i] = i;
    }
    dlxL[0] = colCnt; dlxR[colCnt] = 0;
    for (int i = 1; i < 200; i++) dlxHead[i] = 0;    // head 按「行」索引，必须清够
}
void dlxInsert(int r, int c) {
    ++dlxTot;
    dlxRow[dlxTot] = r; dlxCol[dlxTot] = c;
    dlxU[dlxTot] = dlxU[c]; dlxD[dlxTot] = c;
    dlxD[dlxU[c]] = dlxTot; dlxU[c] = dlxTot;
    if (!dlxHead[r]) { dlxHead[r] = dlxL[dlxTot] = dlxR[dlxTot] = dlxTot; }
    else {
        dlxL[dlxTot] = dlxL[dlxHead[r]];
        dlxR[dlxTot] = dlxHead[r];
        dlxR[dlxL[dlxHead[r]]] = dlxTot;
        dlxL[dlxHead[r]] = dlxTot;
    }
    ++dlxSize[c];
}
void dlxRemoveCol(int c) {
    dlxR[dlxL[c]] = dlxR[c];
    dlxL[dlxR[c]] = dlxL[c];
    for (int i = dlxD[c]; i != c; i = dlxD[i])
        for (int j = dlxR[i]; j != i; j = dlxR[j]) {
            dlxU[dlxD[j]] = dlxU[j];
            dlxD[dlxU[j]] = dlxD[j];
            --dlxSize[dlxCol[j]];
        }
}
void dlxResumeCol(int c) {
    for (int i = dlxU[c]; i != c; i = dlxU[i])
        for (int j = dlxL[i]; j != i; j = dlxL[j]) {
            ++dlxSize[dlxCol[j]];
            dlxU[dlxD[j]] = j;
            dlxD[dlxU[j]] = j;
        }
    dlxR[dlxL[c]] = c;
    dlxL[dlxR[c]] = c;
}
void dlxDance(int depth) {
    if (dlxFound) return;
    if (dlxR[0] == 0) { dlxFound = true; return; }
    int c = dlxR[0];
    for (int i = dlxR[0]; i != 0; i = dlxR[i]) if (dlxSize[i] < dlxSize[c]) c = i;
    dlxRemoveCol(c);
    for (int i = dlxD[c]; i != c; i = dlxD[i]) {
        dlxAns.push_back(dlxRow[i]);
        for (int j = dlxR[i]; j != i; j = dlxR[j]) dlxRemoveCol(dlxCol[j]);
        dlxDance(depth + 1);
        if (dlxFound) return;                 // 关键：找到解直接返回，不要再做任何恢复
        for (int j = dlxL[i]; j != i; j = dlxL[j]) dlxResumeCol(dlxCol[j]);
        dlxAns.pop_back();
    }
    dlxResumeCol(c);
}

// ==================================================================
// 7. 双向 BFS（数字变换：+1 / *2，范围 [0, L]）
// ==================================================================
int bidirBFS(int start, int goal, int L, int maxStep) {
    if (start == goal) return 0;
    vector<int> distS(L + 1, -1), distG(L + 1, -1);
    queue<int> qs, qg;
    distS[start] = 0; qs.push(start);
    distG[goal] = 0;  qg.push(goal);
    auto expand = [&](int x) {
        vector<int> r;
        if (x + 1 <= L) r.push_back(x + 1);
        if (x - 1 >= 0) r.push_back(x - 1);
        if (x * 2 <= L) r.push_back(x * 2);
        if (x % 2 == 0) r.push_back(x / 2);      // 必须有逆操作，双向 BFS 才成立
        return r;
    };
    while (!qs.empty() && !qg.empty()) {
        bool useS = (qs.size() <= qg.size());
        queue<int> &cur = useS ? qs : qg;
        vector<int> &distCur = useS ? distS : distG;
        vector<int> &distOther = useS ? distG : distS;
        int layer = (int)cur.size();
        for (int t = 0; t < layer; t++) {
            int s = cur.front(); cur.pop();
            int d = distCur[s];
            if (d >= maxStep) continue;
            for (int nx : expand(s)) {
                if (distCur[nx] != -1) continue;
                if (distOther[nx] != -1) return d + 1 + distOther[nx];
                distCur[nx] = d + 1;
                cur.push(nx);
            }
        }
    }
    return -1;
}

// ==================================================================
int main() {
    setvbuf(stdout, NULL, _IONBF, 0);

    printf("=== 1. CDQ 分治（三维偏序） ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 300 && ok; iter++) {
            int n = (int)(rng() % 40) + 1;
            cdqValMax = 8;
            vector<array<int,3>> raw(n);
            for (int i = 0; i < n; i++)
                for (int k = 0; k < 3; k++) raw[i][k] = (int)(rng() % cdqValMax) + 1;
            // 排序 + 去重合并
            sort(raw.begin(), raw.end());
            int m = 0;
            for (int i = 0; i < n; ) {
                int j = i;
                while (j < n && raw[j] == raw[i]) ++j;
                ++m;
                cdqP[m].a = raw[i][0]; cdqP[m].b = raw[i][1]; cdqP[m].c = raw[i][2];
                cdqP[m].cnt = j - i; cdqP[m].ans = 0;
                i = j;
            }
            fill(cdqBit, cdqBit + CDQN, 0);
            cdq(1, m);
            // 暴力统计：对每个元素，数出所有各维都 <= 它的元素的 cnt 之和
            // （CDQ 后数组顺序被重排，但元素集合与 ans 绑定不变，所以可以就地比较）
            for (int i = 1; i <= m && ok; i++) {
                ll brute = 0;
                for (int j = 1; j <= m; j++) {
                    if (j == i) continue;
                    if (cdqP[j].a <= cdqP[i].a && cdqP[j].b <= cdqP[i].b && cdqP[j].c <= cdqP[i].c)
                        brute += cdqP[j].cnt;
                }
                if (cdqP[i].ans != brute) {
                    ok = false;
                    printf("    m=%d i=%d got=%d expect=%lld\n", m, i, cdqP[i].ans, brute);
                    break;
                }
            }
        }
        check(ok, "CDQ 三维偏序的 ans 与暴力统计一致（300 组随机三元组）");
    }

    printf("=== 2. 整体二分（静态区间第 k 小） ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 300 && ok; iter++) {
            oN = (int)(rng() % 40) + 1;
            int valMax = 20;
            for (int i = 1; i <= oN; i++) oA[i] = (int)(rng() % valMax) + 1;
            oM = (int)(rng() % 20) + 1;
            vector<array<int,4>> saved(oM + 1);        // 保存原始询问（solveAll 会打乱 oQs）
            for (int i = 1; i <= oM; i++) {
                int l = (int)(rng() % oN) + 1, r = (int)(rng() % oN) + 1;
                if (l > r) swap(l, r);
                int k = (int)(rng() % (r - l + 1)) + 1;
                oQs[i] = {l, r, k, i};
                saved[i] = {l, r, k, i};
            }
            fill(oBit, oBit + OAN, 0);
            solveAll(vector<OQuery>(oQs + 1, oQs + oM + 1), 1, valMax);
            for (int i = 1; i <= oM && ok; i++) {
                vector<int> seg;
                for (int p = saved[i][0]; p <= saved[i][1]; p++) seg.push_back(oA[p]);
                sort(seg.begin(), seg.end());
                int expect = seg[saved[i][2] - 1];
                if (oAns[i] != expect) {
                    ok = false;
                    printf("    询问 %d [%d,%d] k=%d got=%d expect=%d\n",
                           i, saved[i][0], saved[i][1], saved[i][2], oAns[i], expect);
                }
            }
        }
        check(ok, "整体二分求区间第 k 小与排序暴力一致（300 组随机测试）");
    }

    printf("=== 3. FWT（或/与/异或卷积） ===\n");
    {
        bool ok = true;
        const int K = 4, N = 1 << K;
        for (int iter = 0; iter < 100 && ok; iter++) {
            vector<ll> a(N), b(N);
            for (int i = 0; i < N; i++) { a[i] = (ll)(rng() % 10); b[i] = (ll)(rng() % 10); }
            for (int type = 0; type < 3; type++) {
                vector<ll> c = walshConv(a, b, type);
                vector<ll> expect(N, 0);
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++) {
                        int k = (type == 0) ? (i | j) : (type == 1 ? (i & j) : (i ^ j));
                        expect[k] += a[i] * b[j];
                    }
                for (int i = 0; i < N; i++)
                    if (c[i] != expect[i]) { ok = false; printf("    type=%d i=%d got=%lld expect=%lld\n", type, i, c[i], expect[i]); break; }
                if (!ok) break;
            }
        }
        check(ok, "OR / AND / XOR 三种 FWT 卷积与暴力一致（100 组随机数组）");
    }

    printf("=== 4. 子集卷积 ===\n");
    {
        bool ok = true;
        const int K = 3, N = 1 << K;
        for (int iter = 0; iter < 200 && ok; iter++) {
            vector<ll> a(N), b(N);
            for (int i = 0; i < N; i++) { a[i] = (ll)(rng() % 8); b[i] = (ll)(rng() % 8); }
            vector<ll> c = subsetConvolution(a, b);
            for (int S = 0; S < N && ok; S++) {
                ll expect = 0;
                for (int T = 0; T < N; T++)
                    if ((T & S) == T) expect += a[T] * b[S ^ T];    // T ⊆ S
                if (c[S] != expect) { ok = false; printf("    S=%d got=%lld expect=%lld\n", S, c[S], expect); }
            }
        }
        check(ok, "子集卷积 c[S] = Σ_{T⊆S} a[T]b[S\\T] 与暴力一致");
    }

    printf("=== 5. 类欧几里得（floorSum） ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 500 && ok; iter++) {
            ll n = (ll)(rng() % 50) + 1;
            ll m = (ll)(rng() % 30) + 1;
            ll a = (ll)(rng() % 60);
            ll b = (ll)(rng() % 60);
            ll expect = 0;
            for (ll i = 0; i < n; i++) expect += (a * i + b) / m;
            ll got = floorSum(n, m, a, b);
            if (got != expect) { ok = false; printf("    n=%lld m=%lld a=%lld b=%lld got=%lld expect=%lld\n", n, m, a, b, got, expect); }
        }
        check(ok, "floorSum 与暴力求和一致（500 组随机参数）");
    }

    printf("=== 6. min_25 筛 ===\n");
    {
        bool ok = true;
        const int N = 100000;
        vector<ll> pre = brutePrefixMin25(N);
        for (int n : {1, 2, 10, 100, 999, 5000, 20000, 99999, 100000}) {
            initMin25(n);
            ll got = (S25(n, 0) + 1) % MOD25;      // 加 f(1) = 1
            if (got != pre[n]) {
                ok = false;
                printf("    n=%d got=%lld expect=%lld\n", n, got, pre[n]);
                break;
            }
        }
        check(ok, "min_25 筛结果与线性筛前缀和一致（f(p^k)=p^k(p^k-1)，n ≤ 1e5）");
    }

    printf("=== 7. DLX（精确覆盖） ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 300 && ok; iter++) {
            int rows = (int)(rng() % 6) + 1;
            int cols = (int)(rng() % 6) + 1;
            vector<vector<int>> mat(rows, vector<int>(cols, 0));
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++) mat[i][j] = (int)(rng() % 2);
            // 暴力：枚举行子集，检查每列是否恰好被覆盖一次
            bool bruteSol = false;
            for (int mask = 0; mask < (1 << rows); mask++) {
                vector<int> cnt(cols, 0);
                for (int i = 0; i < rows; i++)
                    if ((mask >> i) & 1)
                        for (int j = 0; j < cols; j++) cnt[j] += mat[i][j];
                bool good = true;
                for (int j = 0; j < cols; j++) if (cnt[j] != 1) { good = false; break; }
                if (good) { bruteSol = true; break; }
            }
            // DLX
            dlxInit(cols);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    if (mat[i][j]) dlxInsert(i + 1, j + 1);
            dlxAns.clear();
            dlxFound = false;
            dlxDance(0);
            if (dlxFound != bruteSol) {
                ok = false;
                printf("    rows=%d cols=%d dlx=%d brute=%d\n", rows, cols, (int)dlxFound, (int)bruteSol);
                break;
            }
            // 若有解，验证解的合法性
            if (dlxFound) {
                vector<int> cnt(cols, 0);
                for (int r : dlxAns)
                    for (int j = 0; j < cols; j++) cnt[j] += mat[r - 1][j];
                for (int j = 0; j < cols; j++) if (cnt[j] != 1) { ok = false; break; }
            }
        }
        check(ok, "DLX 精确覆盖的可行性判定与解均与暴力一致（300 组随机矩阵）");
    }

    printf("=== 8. 双向 BFS ===\n");
    {
        bool ok = true;
        for (int iter = 0; iter < 200 && ok; iter++) {
            int L = 200;
            int s = (int)(rng() % 50) + 1, t = (int)(rng() % 50) + 1;
            // 单向 BFS 求最短步数
            vector<int> dist(L + 1, -1);
            queue<int> q;
            dist[s] = 0; q.push(s);
            while (!q.empty()) {
                int u = q.front(); q.pop();
                vector<int> nx;
                if (u + 1 <= L) nx.push_back(u + 1);
                if (u - 1 >= 0) nx.push_back(u - 1);
                if (u * 2 <= L) nx.push_back(u * 2);
                if (u % 2 == 0) nx.push_back(u / 2);
                for (int v : nx) if (dist[v] == -1) { dist[v] = dist[u] + 1; q.push(v); }
            }
            int got = bidirBFS(s, t, L, 100);
            if (got != dist[t]) {
                ok = false;
                printf("    s=%d t=%d bidir=%d bfs=%d\n", s, t, got, dist[t]);
            }
        }
        check(ok, "双向 BFS 与单向 BFS 最短步数一致（200 组）");
    }

    printf("\n===== 通过 %d 项，失败 %d 项 =====\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
