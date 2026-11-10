// 多项式 / 线性代数 / 组合 / 几何 / 博弈 / 单纯形 模板验证
// 编译: g++ -std=c++17 -O2 -o verify_more.exe tools/verify_more.cpp
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <bitset>
#include <random>
#include <chrono>
#include <functional>
#include <tuple>

using namespace std;
using ll = long long;
using i128 = __int128;

static int g_fail = 0, g_pass = 0;
static void check(bool cond, const char *name) {
    if (cond) { ++g_pass; printf("  [ok]   %s\n", name); }
    else      { ++g_fail; printf("  [FAIL] %s\n", name); }
}
static mt19937 rng(20240607);

// 通用整数快速幂（模素数）
ll qpow(ll a, ll e, ll mod) {
    ll r = 1 % mod;
    a %= mod; if (a < 0) a += mod;
    for (; e; e >>= 1, a = (ll)((i128)a * a % mod))
        if (e & 1) r = (ll)((i128)r * a % mod);
    return r;
}
ll my_gcd(ll a, ll b) { return b ? my_gcd(b, a % b) : a; }
ll my_exgcd(ll a, ll b, ll &x, ll &y) {
    if (b == 0) { x = 1; y = 0; return a; }
    ll d = my_exgcd(b, a % b, y, x);
    y -= a / b * x;
    return d;
}

// ==================================================================
// 第 7 章：NTT 与多项式全家桶
// ==================================================================
namespace NTT {
    const int MOD = 998244353;
    const int G   = 3;

    ll qpow(ll a, ll e) {
        ll r = 1 % MOD;
        for (a %= MOD; e; e >>= 1, a = a * a % MOD) if (e & 1) r = r * a % MOD;
        return r;
    }

    void ntt(vector<int> &a, bool inv) {
        int n = (int)a.size();
        for (int i = 1, j = 0; i < n; i++) {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) swap(a[i], a[j]);
        }
        for (int len = 2; len <= n; len <<= 1) {
            int w = (int)qpow(G, (MOD - 1) / len);
            if (inv) w = (int)qpow(w, MOD - 2);
            for (int i = 0; i < n; i += len) {
                ll wn = 1;
                for (int j = 0; j < len / 2; j++) {
                    int u = a[i + j];
                    int v = (int)(wn * a[i + j + len / 2] % MOD);
                    int x = u + v; if (x >= MOD) x -= MOD;
                    int y = u - v; if (y < 0) y += MOD;
                    a[i + j] = x;
                    a[i + j + len / 2] = y;
                    wn = wn * w % MOD;
                }
            }
        }
        if (inv) {
            ll inv_n = qpow(n, MOD - 2);
            for (int &x : a) x = (int)(x * inv_n % MOD);
        }
    }

    vector<int> convolution(vector<int> a, vector<int> b) {
        int n = (int)a.size(), m = (int)b.size();
        if (!n || !m) return {};
        if (min(n, m) <= 60) {
            vector<int> c(n + m - 1, 0);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                    c[i + j] = (c[i + j] + (ll)a[i] * b[j]) % MOD;
            return c;
        }
        int sz = 1;
        while (sz < n + m - 1) sz <<= 1;
        a.resize(sz); b.resize(sz);
        ntt(a, false); ntt(b, false);
        for (int i = 0; i < sz; i++) a[i] = (int)((ll)a[i] * b[i] % MOD);
        ntt(a, true);
        a.resize(n + m - 1);
        return a;
    }
}

vector<int> invs;
void init_invs(int n) {
    invs.assign(n + 1, 0);
    if (n >= 1) invs[1] = 1;
    for (int i = 2; i <= n; i++)
        invs[i] = (ll)(NTT::MOD - NTT::MOD / i) * invs[NTT::MOD % i] % NTT::MOD;
}

vector<int> poly_derivative(const vector<int> &f) {
    int n = (int)f.size();
    vector<int> d(n > 0 ? n - 1 : 0, 0);
    for (int i = 1; i < n; i++) d[i - 1] = (ll)f[i] * i % NTT::MOD;
    return d;
}

vector<int> poly_integral(const vector<int> &f) {
    int n = (int)f.size();
    vector<int> g(n + 1, 0);
    for (int i = 0; i < n; i++) g[i + 1] = (ll)f[i] * invs[i + 1] % NTT::MOD;
    return g;
}

vector<int> poly_inv(const vector<int> &f, int n) {
    vector<int> g(1, (int)NTT::qpow(f[0], NTT::MOD - 2));
    int cur = 1;
    while (cur < n) {
        int nxt = min(cur << 1, n);
        int sz = 1;
        while (sz < 2 * nxt) sz <<= 1;
        vector<int> a(f.begin(), f.begin() + min((int)f.size(), nxt));
        a.resize(sz);
        vector<int> b(g);
        b.resize(sz);
        NTT::ntt(a, false); NTT::ntt(b, false);
        for (int i = 0; i < sz; i++) {
            ll fg = (ll)a[i] * b[i] % NTT::MOD;
            a[i] = (int)((2 - fg + NTT::MOD) % NTT::MOD * b[i] % NTT::MOD);
        }
        NTT::ntt(a, true);
        g.assign(a.begin(), a.begin() + nxt);
        cur = nxt;
    }
    return g;
}

vector<int> poly_ln(const vector<int> &f, int n) {
    vector<int> d = poly_derivative(f);
    vector<int> invf = poly_inv(f, n);
    vector<int> q = NTT::convolution(d, invf);
    q.resize(n > 0 ? n - 1 : 0);
    vector<int> res = poly_integral(q);
    res.resize(n);
    return res;
}

vector<int> poly_exp(const vector<int> &f, int n) {
    vector<int> g(1, 1);
    int cur = 1;
    while (cur < n) {
        int nxt = min(cur << 1, n);
        vector<int> lg = poly_ln(g, nxt);
        vector<int> t(nxt, 0);
        for (int i = 0; i < nxt; i++) {
            int fi = (i < (int)f.size()) ? f[i] : 0;
            t[i] = ((fi - lg[i]) % NTT::MOD + NTT::MOD) % NTT::MOD;
        }
        t[0] = (t[0] + 1) % NTT::MOD;
        g = NTT::convolution(g, t);
        g.resize(nxt);
        cur = nxt;
    }
    return g;
}

vector<int> poly_sqrt(const vector<int> &f, int n, int g0 = 1) {
    vector<int> g(1, g0);
    int cur = 1;
    while (cur < n) {
        int nxt = min(cur << 1, n);
        vector<int> invg = poly_inv(g, nxt);
        vector<int> t = NTT::convolution(f, invg);
        t.resize(nxt);
        g.resize(nxt);
        const int INV2 = (NTT::MOD + 1) / 2;
        for (int i = 0; i < nxt; i++)
            g[i] = (ll)(g[i] + t[i]) % NTT::MOD * INV2 % NTT::MOD;
        cur = nxt;
    }
    return g;
}

vector<int> poly_pow(const vector<int> &f, ll k, int n) {
    vector<int> lg = poly_ln(f, n);
    for (int &x : lg) x = (int)((ll)x * (k % NTT::MOD) % NTT::MOD);
    return poly_exp(lg, n);
}

// 朴素卷积
vector<int> naive_conv(const vector<int> &a, const vector<int> &b) {
    if (a.empty() || b.empty()) return {};
    vector<int> c(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); i++)
        for (size_t j = 0; j < b.size(); j++)
            c[i + j] = (c[i + j] + (ll)a[i] * b[j]) % NTT::MOD;
    return c;
}

bool eq_prefix(const vector<int> &a, const vector<int> &b, int n) {
    for (int i = 0; i < n; i++) {
        int x = (i < (int)a.size()) ? a[i] : 0;
        int y = (i < (int)b.size()) ? b[i] : 0;
        if (x != y) return false;
    }
    return true;
}

// ==================================================================
// 拉格朗日插值 / 线性递推
// ==================================================================
ll lagrange(const vector<ll> &x, const vector<ll> &y, ll k, ll MOD) {
    int n = (int)x.size();
    ll res = 0;
    for (int i = 0; i < n; i++) {
        ll num = 1, den = 1;
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            num = num * ((k - x[j]) % MOD + MOD) % MOD;
            den = den * ((x[i] - x[j]) % MOD + MOD) % MOD;
        }
        res = (res + y[i] % MOD * num % MOD * qpow(den, MOD - 2, MOD)) % MOD;
    }
    return res;
}

ll lagrange_consecutive(const vector<ll> &y, ll k, ll MOD) {
    int n = (int)y.size() - 1;
    if (k >= 0 && k <= n) return y[(int)k] % MOD;
    vector<ll> pre(n + 2, 1), suf(n + 2, 1), fac(n + 2, 1), ifac(n + 2, 1);
    for (int i = 0; i <= n; i++) pre[i + 1] = pre[i] * ((k - i) % MOD + MOD) % MOD;
    for (int i = n; i >= 0; i--) suf[i] = suf[i + 1] * ((k - i) % MOD + MOD) % MOD;
    for (int i = 1; i <= n + 1; i++) fac[i] = fac[i - 1] * i % MOD;
    ifac[n + 1] = qpow(fac[n + 1], MOD - 2, MOD);
    for (int i = n + 1; i >= 1; i--) ifac[i - 1] = ifac[i] * i % MOD;
    ll res = 0;
    for (int i = 0; i <= n; i++) {
        ll num = pre[i] * suf[i + 1] % MOD;
        ll den = ifac[i] * ifac[n - i] % MOD;
        ll term = y[i] % MOD * num % MOD * den % MOD;
        if ((n - i) & 1) res = (res - term + MOD) % MOD;
        else res = (res + term) % MOD;
    }
    return res;
}

ll kitamasa(vector<ll> a, vector<ll> c, ll n, ll MOD) {
    int k = (int)c.size();
    if (n < k) return a[(int)n] % MOD;
    for (ll &v : a) v = (v % MOD + MOD) % MOD;
    for (ll &v : c) v = (v % MOD + MOD) % MOD;

    auto mul = [&](const vector<ll> &u, const vector<ll> &v) {
        vector<ll> t(2 * k, 0);
        for (int i = 0; i < k; i++) {
            if (!u[i]) continue;
            for (int j = 0; j < k; j++)
                t[i + j] = (t[i + j] + u[i] * v[j]) % MOD;
        }
        for (int i = 2 * k - 2; i >= k; i--) {
            if (!t[i]) continue;
            for (int j = 1; j <= k; j++)
                t[i - j] = (t[i - j] + t[i] * c[j - 1]) % MOD;
        }
        t.resize(k);
        return t;
    };

    vector<ll> res(k, 0), base(k, 0);
    res[0] = 1;
    if (k == 1) base[0] = c[0] % MOD;
    else base[1] = 1;
    for (ll e = n; e; e >>= 1, base = mul(base, base))
        if (e & 1) res = mul(res, base);

    ll ans = 0;
    for (int i = 0; i < k; i++) ans = (ans + res[i] * a[i]) % MOD;
    return ans;
}

vector<ll> berlekamp_massey(const vector<ll> &s, ll MOD) {
    vector<ll> C(1, 1), B(1, 1);
    int L = 0, m = 1;
    ll b = 1;
    for (int n = 0; n < (int)s.size(); n++) {
        ll d = 0;
        for (int i = 0; i <= L && i < (int)C.size(); i++)
            d = (d + C[i] * s[n - i]) % MOD;
        d = (d % MOD + MOD) % MOD;
        if (d == 0) { ++m; continue; }
        vector<ll> T = C;
        ll coef = d * qpow(b, MOD - 2, MOD) % MOD;
        if ((int)C.size() < (int)B.size() + m) C.resize(B.size() + m, 0);
        for (int i = 0; i < (int)B.size(); i++)
            C[i + m] = ((C[i + m] - coef * B[i]) % MOD + MOD) % MOD;
        if (2 * L <= n) { L = n + 1 - L; B = T; b = d; m = 1; }
        else ++m;
    }
    vector<ll> res(L, 0);
    for (int i = 0; i < L; i++) res[i] = (MOD - C[i + 1] % MOD) % MOD;
    return res;
}

// ==================================================================
// 第 6 章：矩阵 / 高斯消元 / 行列式 / 线性基 / 矩阵树
// ==================================================================
const ll GMOD = 998244353;

// 模意义高斯消元
int gauss_mod(vector<vector<ll>> a, vector<ll> b, vector<ll> &x, ll P) {
    int n = (int)a.size(), m = (int)a[0].size();
    vector<int> where(m, -1);
    int row = 0;
    for (int col = 0; col < m && row < n; col++) {
        int piv = -1;
        for (int i = row; i < n; i++) if (a[i][col] % P) { piv = i; break; }
        if (piv < 0) continue;
        swap(a[piv], a[row]); swap(b[piv], b[row]);
        ll inv = qpow(a[row][col], P - 2, P);
        for (int j = col; j < m; j++) a[row][j] = a[row][j] % P * inv % P;
        b[row] = b[row] % P * inv % P;
        for (int i = 0; i < n; i++) {
            if (i == row || a[i][col] % P == 0) continue;
            ll f = a[i][col] % P;
            for (int j = col; j < m; j++)
                a[i][j] = ((a[i][j] - f * a[row][j]) % P + P) % P;
            b[i] = ((b[i] - f * b[row]) % P + P) % P;
        }
        where[col] = row++;
    }
    for (int i = row; i < n; i++) if (b[i] % P) return -1;
    x.assign(m, 0);
    for (int col = 0; col < m; col++)
        if (where[col] != -1) x[col] = b[where[col]] % P;
    return m - row;
}

// 异或方程组
int gauss_xor(vector<unsigned long long> coef, vector<unsigned long long> rhs, int m) {
    int n = (int)coef.size(), rank = 0;
    for (int col = 0; col < m && rank < n; col++) {
        int piv = -1;
        for (int i = rank; i < n; i++) if ((coef[i] >> col) & 1ULL) { piv = i; break; }
        if (piv < 0) continue;
        swap(coef[piv], coef[rank]); swap(rhs[piv], rhs[rank]);
        for (int i = 0; i < n; i++) {
            if (i == rank || !((coef[i] >> col) & 1ULL)) continue;
            coef[i] ^= coef[rank]; rhs[i] ^= rhs[rank];
        }
        rank++;
    }
    for (int i = rank; i < n; i++) if (rhs[i]) return -1;
    return m - rank;
}

// 模素数行列式
ll det_mod(vector<vector<ll>> a, ll P) {
    int n = (int)a.size();
    ll res = 1;
    for (int i = 0; i < n; i++) {
        int piv = -1;
        for (int j = i; j < n; j++) if (a[j][i] % P) { piv = j; break; }
        if (piv < 0) return 0;
        if (piv != i) { swap(a[piv], a[i]); res = (P - res) % P; }
        res = res * a[i][i] % P;
        ll inv = qpow(a[i][i], P - 2, P);
        for (int j = i + 1; j < n; j++) {
            if (a[j][i] % P == 0) continue;
            ll f = a[j][i] % P * inv % P;
            for (int k = i; k < n; k++)
                a[j][k] = ((a[j][k] - f * a[i][k]) % P + P) % P;
        }
    }
    return res;
}

// 矩阵树定理
ll count_spanning_trees(int n, const vector<pair<int,int>> &edges, ll P) {
    if (n == 1) return 1 % P;
    vector<vector<ll>> L(n, vector<ll>(n, 0));
    for (const auto &e : edges) {
        int u = e.first, v = e.second;
        if (u == v) continue;
        L[u][u]++; L[v][v]++;
        L[u][v] = (L[u][v] - 1 + P) % P;
        L[v][u] = (L[v][u] - 1 + P) % P;
    }
    vector<vector<ll>> a(n - 1, vector<ll>(n - 1));
    for (int i = 0; i + 1 < n; i++)
        for (int j = 0; j + 1 < n; j++) a[i][j] = L[i][j] % P;
    return det_mod(a, P);
}

// 线性基
struct LinearBasis {
    static const int LOG = 62;
    ll b[LOG + 1];
    bool hasZero;
    LinearBasis() { memset(b, 0, sizeof(b)); hasZero = false; }
    bool insert(ll x) {
        for (int i = LOG; i >= 0; i--) {
            if (!((x >> i) & 1)) continue;
            if (!b[i]) { b[i] = x; return true; }
            x ^= b[i];
        }
        hasZero = true;
        return false;
    }
    ll max_xor() const {
        ll res = 0;
        for (int i = LOG; i >= 0; i--) res = max(res, res ^ b[i]);
        return res;
    }
    bool can(ll x) const {
        for (int i = LOG; i >= 0; i--)
            if ((x >> i) & 1) { if (!b[i]) return false; x ^= b[i]; }
        return true;
    }
    void normalize() {
        for (int i = LOG; i >= 0; i--) {
            if (!b[i]) continue;
            for (int j = i - 1; j >= 0; j--)
                if ((b[i] >> j) & 1) b[i] ^= b[j];
        }
    }
    ll kth(ll k) const {
        vector<ll> v;
        for (int i = 0; i <= LOG; i++) if (b[i]) v.push_back(b[i]);
        int m = (int)v.size();
        if (m < 62 && k > (1LL << m)) return -1;
        --k;
        ll res = 0;
        for (int i = 0; i < m; i++) if ((k >> i) & 1) res ^= v[i];
        return res;
    }
};

// 矩阵快速幂（斐波那契）
template <int MOD>
struct Matrix {
    int n, m;
    vector<vector<ll>> a;
    Matrix(int n = 0, int m = 0) : n(n), m(m), a(n, vector<ll>(m, 0)) {}
    static Matrix identity(int n) {
        Matrix I(n, n);
        for (int i = 0; i < n; i++) I.a[i][i] = 1;
        return I;
    }
    Matrix operator*(const Matrix &o) const {
        Matrix r(n, o.m);
        for (int i = 0; i < n; i++)
            for (int k = 0; k < m; k++) {
                if (!a[i][k]) continue;
                ll v = a[i][k];
                for (int j = 0; j < o.m; j++)
                    r.a[i][j] = (r.a[i][j] + v * o.a[k][j]) % MOD;
            }
        return r;
    }
    Matrix pow(ll e) const {
        Matrix r = identity(n), b = *this;
        for (; e; e >>= 1, b = b * b) if (e & 1) r = r * b;
        return r;
    }
};

template <int MOD>
ll fib_matrix(ll n) {
    if (n <= 0) return 0;
    Matrix<MOD> T(2, 2);
    T.a[0][0] = 1; T.a[0][1] = 1;
    T.a[1][0] = 1; T.a[1][1] = 0;
    return T.pow(n).a[0][1];
}

// ==================================================================
int main() {
    printf("=== 1. NTT 卷积 ===\n");
    {
        bool ok = true;
        for (int it = 0; it < 40; it++) {
            int n = (int)(rng() % 80) + 1, m = (int)(rng() % 80) + 1;
            vector<int> a(n), b(m);
            for (int &x : a) x = (int)(rng() % NTT::MOD);
            for (int &x : b) x = (int)(rng() % NTT::MOD);
            vector<int> c1 = NTT::convolution(a, b);
            vector<int> c2 = naive_conv(a, b);
            if (c1 != c2) { ok = false; printf("    n=%d m=%d 不一致\n", n, m); break; }
        }
        check(ok, "NTT 卷积与朴素卷积一致（小规模，含暴力分支）");

        // 大规模（强制走 NTT 分支）
        bool ok2 = true;
        for (int it = 0; it < 5; it++) {
            int n = 300 + (int)(rng() % 200), m = 250 + (int)(rng() % 200);
            vector<int> a(n), b(m);
            for (int &x : a) x = (int)(rng() % NTT::MOD);
            for (int &x : b) x = (int)(rng() % NTT::MOD);
            vector<int> c1 = NTT::convolution(a, b);
            vector<int> c2 = naive_conv(a, b);
            if (c1 != c2) { ok2 = false; break; }
        }
        check(ok2, "NTT 卷积与朴素卷积一致（大规模，走 NTT 分支）");
    }

    printf("=== 2. 多项式求逆 / ln / exp / sqrt / pow ===\n");
    {
        init_invs(1 << 14);
        const int n = 64;

        // 求逆：f * inv(f) ≡ 1 (mod x^n)
        bool ok = true;
        for (int it = 0; it < 10; it++) {
            vector<int> f(n);
            f[0] = (int)(rng() % (NTT::MOD - 1)) + 1;
            for (int i = 1; i < n; i++) f[i] = (int)(rng() % NTT::MOD);
            vector<int> g = poly_inv(f, n);
            vector<int> prod = NTT::convolution(f, g);
            prod.resize(n + 1);
            bool good = (prod[0] == 1);
            for (int i = 1; i < n; i++) if (prod[i] != 0) good = false;
            if (!good) { ok = false; break; }
        }
        check(ok, "poly_inv: f * f^{-1} ≡ 1 (mod x^64)");

        // ln / exp 互逆
        bool ok2 = true;
        for (int it = 0; it < 8; it++) {
            vector<int> f(n, 0);
            f[0] = 1;
            for (int i = 1; i < n; i++) f[i] = (int)(rng() % NTT::MOD);
            vector<int> lg = poly_ln(f, n);
            vector<int> ex = poly_exp(lg, n);
            if (!eq_prefix(ex, f, n)) { ok2 = false; break; }
        }
        check(ok2, "poly_exp(poly_ln(f)) == f (mod x^64)");

        // exp 的正确性：对 f[0]=0，验证 exp(f)' = exp(f) * f'
        bool ok3 = true;
        for (int it = 0; it < 8; it++) {
            vector<int> f(n, 0);
            for (int i = 1; i < n; i++) f[i] = (int)(rng() % NTT::MOD);
            vector<int> e = poly_exp(f, n);
            vector<int> de = poly_derivative(e);
            vector<int> df = poly_derivative(f);
            vector<int> rhs = NTT::convolution(e, df);
            rhs.resize(n - 1);
            if (!eq_prefix(de, rhs, n - 1)) { ok3 = false; break; }
        }
        check(ok3, "poly_exp 满足 (e^f)' = e^f · f'");

        // sqrt：g^2 = f
        bool ok4 = true;
        for (int it = 0; it < 8; it++) {
            vector<int> g(n);
            g[0] = 1;
            for (int i = 1; i < n; i++) g[i] = (int)(rng() % NTT::MOD);
            vector<int> f = NTT::convolution(g, g);
            f.resize(n);
            vector<int> s = poly_sqrt(f, n);
            vector<int> ss = NTT::convolution(s, s);
            if (!eq_prefix(ss, f, n)) { ok4 = false; break; }
        }
        check(ok4, "poly_sqrt: (sqrt(f))^2 == f (mod x^64)");

        // pow：f^3 = f*f*f
        bool ok5 = true;
        for (int it = 0; it < 6; it++) {
            vector<int> f(n, 0);
            f[0] = 1;
            for (int i = 1; i < n; i++) f[i] = (int)(rng() % NTT::MOD);
            vector<int> p = poly_pow(f, 3, n);
            vector<int> t1 = NTT::convolution(f, f); t1.resize(n);
            vector<int> t2 = NTT::convolution(t1, f); t2.resize(n);
            if (!eq_prefix(p, t2, n)) { ok5 = false; break; }
        }
        check(ok5, "poly_pow: f^3 == f·f·f (mod x^64)");
    }

    printf("=== 3. 拉格朗日插值 ===\n");
    {
        const ll MOD = 1000000007;
        bool ok = true;
        for (int it = 0; it < 20; it++) {
            int deg = (int)(rng() % 5);                  // 次数 0..4
            vector<ll> c(deg + 1);
            for (ll &v : c) v = (ll)(rng() % MOD);
            auto eval = [&](ll x) {
                ll r = 0, p = 1;
                for (int i = 0; i <= deg; i++) { r = (r + c[i] * p) % MOD; p = p * x % MOD; }
                return r;
            };
            vector<ll> xs, ys;
            for (int i = 0; i <= deg; i++) { xs.push_back(i); ys.push_back(eval(i)); }
            for (int q = 0; q < 3; q++) {
                ll k = (ll)(rng() % 100) + 50;
                if (lagrange(xs, ys, k, MOD) != eval(k)) { ok = false; break; }
                if (lagrange_consecutive(ys, k, MOD) != eval(k)) { ok = false; break; }
            }
            if (!ok) break;
        }
        check(ok, "lagrange / lagrange_consecutive 与直接求值一致");

        // 大 k（不落在插值点上）
        {
            int deg = 4;
            vector<ll> c(deg + 1, 0);
            c[0] = 7; c[1] = 3; c[2] = 5; c[3] = 1; c[4] = 2;
            auto eval = [&](ll x) {
                ll r = 0, p = 1;
                for (int i = 0; i <= deg; i++) { r = (r + c[i] * p) % MOD; p = p * x % MOD; }
                return r;
            };
            vector<ll> ys;
            for (int i = 0; i <= deg; i++) ys.push_back(eval(i));
            ll k = 123456789012345LL % MOD;
            check(lagrange_consecutive(ys, k, MOD) == eval(k), "连续点插值在巨大 k 处正确");
        }
    }

    printf("=== 4. Kitamasa 与 Berlekamp-Massey ===\n");
    {
        const ll MOD = 1000000007;
        bool ok = true;
        for (int it = 0; it < 10; it++) {
            int k = (int)(rng() % 5) + 1;
            vector<ll> c(k), a(k);
            for (ll &v : c) v = (ll)(rng() % MOD);
            for (ll &v : a) v = (ll)(rng() % MOD);
            vector<ll> s = a;
            for (int i = k; i <= 60; i++) {
                ll v = 0;
                for (int j = 0; j < k; j++) v = (v + c[j] * s[i - 1 - j]) % MOD;
                s.push_back(v);
            }
            for (int n = 0; n <= 60; n++)
                if (kitamasa(a, c, n, MOD) != s[n]) { ok = false; printf("    k=%d n=%d\n", k, n); break; }
            if (!ok) break;
        }
        check(ok, "kitamasa 与朴素递推一致（n ≤ 60）");

        // BM：从序列反推递推，再用 kitamasa 预测后续项
        bool ok2 = true;
        for (int it = 0; it < 10; it++) {
            int k = (int)(rng() % 4) + 2;
            vector<ll> c(k), a(k);
            for (ll &v : c) v = (ll)(rng() % MOD);
            for (ll &v : a) v = (ll)(rng() % MOD);
            vector<ll> s = a;
            for (int i = k; i <= 120; i++) {
                ll v = 0;
                for (int j = 0; j < k; j++) v = (v + c[j] * s[i - 1 - j]) % MOD;
                s.push_back(v);
            }
            vector<ll> rc = berlekamp_massey(vector<ll>(s.begin(), s.begin() + 100), MOD);
            if (rc.empty() || (int)rc.size() > k + 2) { ok2 = false; break; }
            // 用反推出的递推预测后续项
            vector<ll> ra(rc.size());
            for (size_t i = 0; i < rc.size(); i++) ra[i] = s[i];
            for (int n = 100; n <= 120; n++)
                if (kitamasa(ra, rc, n, MOD) != s[n]) { ok2 = false; printf("    k=%d n=%d\n", k, n); break; }
            if (!ok2) break;
        }
        check(ok2, "BM 反推的递推能正确预测后续项");

        // 大 n：与矩阵快速幂对比斐波那契
        {
            bool ok3 = true;
            vector<ll> c = {1, 1}, a = {0, 1};
            for (ll n : {100LL, 1000LL, 12345LL, 1000000LL}) {
                ll k1 = kitamasa(a, c, n, 1000000007);
                ll k2 = fib_matrix<1000000007>(n);
                if (k1 != k2) { ok3 = false; }
            }
            check(ok3, "kitamasa 与矩阵快速幂求斐波那契一致");
        }
    }

    printf("=== 5. 矩阵快速幂 ===\n");
    {
        bool ok = true;
        vector<ll> f(200);
        f[0] = 0; f[1] = 1;
        for (int i = 2; i < 200; i++) f[i] = (f[i - 1] + f[i - 2]) % 1000000007;
        for (int i = 0; i < 200; i++)
            if (fib_matrix<1000000007>(i) != f[i]) { ok = false; break; }
        check(ok, "矩阵快速幂求斐波那契与递推一致");
    }

    printf("=== 6. 高斯消元（模意义） ===\n");
    {
        const ll P = 1000000007;
        bool ok = true;
        for (int it = 0; it < 30; it++) {
            int n = (int)(rng() % 8) + 1;
            vector<vector<ll>> A(n, vector<ll>(n));
            vector<ll> x0(n), b(n, 0);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) A[i][j] = (ll)(rng() % P);
                x0[i] = (ll)(rng() % P);
            }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) b[i] = (b[i] + A[i][j] * x0[j]) % P;
            vector<ll> x;
            int freev = gauss_mod(A, b, x, P);
            if (freev != 0) continue;                 // 跳过奇异矩阵
            if (x != x0) { ok = false; printf("    n=%d 解不一致\n", n); break; }
        }
        check(ok, "模意义高斯消元解出构造解");

        // 无解检测
        {
            vector<vector<ll>> A = {{1, 1}, {2, 2}};
            vector<ll> b = {1, 3}, x;
            check(gauss_mod(A, b, x, P) == -1, "模意义高斯消元能检出无解");
        }
    }

    printf("=== 7. 异或方程组 ===\n");
    {
        bool ok = true;
        for (int it = 0; it < 40; it++) {
            int m = (int)(rng() % 8) + 1, n = (int)(rng() % 8) + 1;
            vector<unsigned long long> coef(n), rhs(n, 0);
            unsigned long long x0 = rng() & ((m == 64) ? ~0ULL : ((1ULL << m) - 1));
            for (int i = 0; i < n; i++) {
                coef[i] = rng() & ((1ULL << m) - 1);
                rhs[i] = (unsigned long long)__builtin_popcountll(coef[i] & x0) & 1ULL;
            }
            int fv = gauss_xor(coef, rhs, m);
            if (fv < 0) { ok = false; break; }
        }
        check(ok, "异或方程组对构造解均可解");
    }

    printf("=== 8. 行列式 ===\n");
    {
        const ll P = 1000000007;
        // 与 O(n!) 暴力对比（n <= 6）
        auto brute_det = [](vector<vector<ll>> a, ll P) {
            int n = (int)a.size();
            vector<int> p(n);
            for (int i = 0; i < n; i++) p[i] = i;
            ll res = 0;
            do {
                ll prod = 1;
                int inv = 0;
                for (int i = 0; i < n; i++)
                    for (int j = i + 1; j < n; j++) if (p[i] > p[j]) ++inv;
                for (int i = 0; i < n; i++) prod = prod * a[i][p[i]] % P;
                if (inv & 1) res = (res - prod + P) % P;
                else res = (res + prod) % P;
            } while (next_permutation(p.begin(), p.end()));
            return res;
        };
        bool ok = true;
        for (int it = 0; it < 30; it++) {
            int n = (int)(rng() % 6) + 1;
            vector<vector<ll>> a(n, vector<ll>(n));
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) a[i][j] = (ll)(rng() % 1000);
            if (det_mod(a, P) != brute_det(a, P)) { ok = false; printf("    n=%d\n", n); break; }
        }
        check(ok, "det_mod 与 n! 暴力展开一致");

        // 奇异矩阵行列式为 0
        {
            vector<vector<ll>> a = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
            check(det_mod(a, P) == 0, "奇异矩阵行列式为 0");
        }
    }

    printf("=== 9. 矩阵树定理 ===\n");
    {
        const ll P = 1000000007;
        // 与暴力枚举边子集对比（小图）
        auto brute_trees = [](int n, const vector<pair<int,int>> &edges) {
            int m = (int)edges.size();
            ll cnt = 0;
            for (int mask = 0; mask < (1 << m); mask++) {
                if (__builtin_popcount(mask) != n - 1) continue;
                vector<int> fa(n);
                for (int i = 0; i < n; i++) fa[i] = i;
                function<int(int)> find = [&](int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); };
                bool ok = true;
                for (int i = 0; i < m && ok; i++) {
                    if (!((mask >> i) & 1)) continue;
                    int u = find(edges[i].first), v = find(edges[i].second);
                    if (u == v) ok = false;
                    else fa[u] = v;
                }
                if (ok) ++cnt;
            }
            return cnt;
        };
        bool ok = true;
        for (int it = 0; it < 30; it++) {
            int n = (int)(rng() % 4) + 2;
            vector<pair<int,int>> edges;
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    if (rng() % 2) edges.push_back({i, j});
            if (edges.empty()) continue;
            ll a = count_spanning_trees(n, edges, P);
            ll b = brute_trees(n, edges);
            if (a != b % P) { ok = false; printf("    n=%d 得到 %lld 期望 %lld\n", n, a, b); break; }
        }
        check(ok, "矩阵树定理与暴力枚举生成树一致");

        // 完全图 K_n 的生成树个数 = n^{n-2}
        {
            bool ok2 = true;
            for (int n = 2; n <= 6; n++) {
                vector<pair<int,int>> edges;
                for (int i = 0; i < n; i++)
                    for (int j = i + 1; j < n; j++) edges.push_back({i, j});
                ll expect = 1;
                for (int i = 0; i < n - 2; i++) expect = expect * n % 1000000007;
                if (count_spanning_trees(n, edges, 1000000007) != expect) ok2 = false;
            }
            check(ok2, "K_n 生成树数 = n^{n-2}（Cayley 公式）");
        }
    }

    printf("=== 10. 线性基 ===\n");
    {
        bool ok = true, ok2 = true, ok3 = true;
        for (int it = 0; it < 200; it++) {
            int n = (int)(rng() % 10) + 1;
            vector<ll> v(n);
            for (ll &x : v) x = (ll)(rng() % 1024);
            LinearBasis lb;
            for (ll x : v) lb.insert(x);

            // 最大异或和 = 暴力枚举子集
            ll best = 0;
            for (int mask = 0; mask < (1 << n); mask++) {
                ll x = 0;
                for (int i = 0; i < n; i++) if ((mask >> i) & 1) x ^= v[i];
                best = max(best, x);
            }
            if (lb.max_xor() != best) { ok = false; break; }

            // can()：暴力判断可达性
            vector<bool> reach(1024, false);
            for (int mask = 0; mask < (1 << n); mask++) {
                ll x = 0;
                for (int i = 0; i < n; i++) if ((mask >> i) & 1) x ^= v[i];
                reach[x] = true;
            }
            vector<ll> all;
            for (int x = 0; x < 1024; x++) if (reach[x]) all.push_back(x);
            sort(all.begin(), all.end());
            for (int x = 0; x < 1024; x++)
                if (lb.can(x) != reach[x]) { ok2 = false; break; }

            // kth
            LinearBasis lb2 = lb;
            lb2.normalize();
            for (size_t i = 0; i < all.size(); i++) {
                if (lb2.kth((ll)i + 1) != all[i]) { ok3 = false; printf("    第 %zu 小错误\n", i + 1); break; }
            }
            if (!ok2 || !ok3) break;
        }
        check(ok, "线性基最大异或和与暴力一致");
        check(ok2, "线性基 can() 与暴力可达性一致");
        check(ok3, "线性基第 k 小与暴力排序一致");
    }

    printf("\n===== 通过 %d 项，失败 %d 项 =====\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
