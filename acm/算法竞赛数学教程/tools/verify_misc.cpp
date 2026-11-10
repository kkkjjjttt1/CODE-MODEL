// 组合 / 几何 / 博弈 / 单纯形 / Burnside / 数值方法 模板验证
// 编译: g++ -std=c++17 -O2 -o verify_misc.exe tools/verify_misc.cpp
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <unordered_map>
#include <random>
#include <chrono>
#include <functional>

using namespace std;
using ll = long long;
using i128 = __int128;

static int g_fail = 0, g_pass = 0;
static void check(bool cond, const char *name) {
    if (cond) { ++g_pass; printf("  [ok]   %s\n", name); }
    else      { ++g_fail; printf("  [FAIL] %s\n", name); }
}
static mt19937 rng(20240607);

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
ll inv_exgcd(ll a, ll mod) {
    ll x, y;
    ll d = my_exgcd(a, mod, x, y);
    if (d != 1) return -1;
    return (x % mod + mod) % mod;
}

// ==================================================================
// 第 4 章：组合数
// ==================================================================
const int MAXF = 2000005;
ll fac[MAXF], ifac[MAXF];
ll CMOD;

void init_comb(int n, ll mod) {
    CMOD = mod;
    fac[0] = 1;
    for (int i = 1; i <= n; i++) fac[i] = fac[i - 1] * i % mod;
    ifac[n] = qpow(fac[n], mod - 2, mod);
    for (int i = n; i >= 1; i--) ifac[i - 1] = ifac[i] * i % mod;
}
ll Csmall(int n, int k) {
    if (k < 0 || k > n) return 0;
    return fac[n] * ifac[k] % CMOD * ifac[n - k] % CMOD;
}

ll lucas(ll n, ll k, ll p) {
    if (k == 0) return 1 % p;
    return Csmall((int)(n % p), (int)(k % p)) * lucas(n / p, k / p, p) % p;
}

// exLucas
ll fac_mod(ll n, ll p, ll pk) {
    if (n == 0) return 1 % pk;
    ll res = 1;
    for (ll i = 2; i <= pk; i++)
        if (i % p) res = (ll)((i128)res * i % pk);
    res = qpow(res, n / pk, pk);
    for (ll i = 2; i <= n % pk; i++)
        if (i % p) res = (ll)((i128)res * i % pk);
    return (ll)((i128)res * fac_mod(n / p, p, pk) % pk);
}
ll C_mod_pk(ll n, ll m, ll p, ll pk) {
    if (m < 0 || m > n) return 0;
    ll a = fac_mod(n, p, pk), b = fac_mod(m, p, pk), c = fac_mod(n - m, p, pk);
    ll k = 0;
    for (ll i = n; i; i /= p) k += i / p;
    for (ll i = m; i; i /= p) k -= i / p;
    for (ll i = n - m; i; i /= p) k -= i / p;
    ll res = a * inv_exgcd(b % pk, pk) % pk * inv_exgcd(c % pk, pk) % pk;
    res = (ll)((i128)res * qpow(p, k, pk) % pk);
    return res;
}
ll exlucas(ll n, ll m, ll mod) {
    if (m < 0 || m > n) return 0;
    ll res = 0, M = 1, tmp = mod;
    for (ll p = 2; p * p <= tmp; p++) {
        if (tmp % p) continue;
        ll pk = 1;
        while (tmp % p == 0) { tmp /= p; pk *= p; }
        ll r = C_mod_pk(n, m, p, pk);
        ll inv = inv_exgcd(M % pk, pk);
        ll t = (ll)((i128)((r - res) % pk + pk) % pk * inv % pk);
        res = (ll)(((i128)res + (i128)t * M) % (M * pk));
        M *= pk;
    }
    if (tmp > 1) {
        ll p = tmp, pk = tmp;
        ll r = C_mod_pk(n, m, p, pk);
        ll inv = inv_exgcd(M % pk, pk);
        ll t = (ll)((i128)((r - res) % pk + pk) % pk * inv % pk);
        res = (ll)(((i128)res + (i128)t * M) % (M * pk));
        M *= pk;
    }
    return (res % M + M) % M;
}

// 大整数组合数（i128，仅用于 n <= 60 的精确对照）
i128 bigC(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    i128 r = 1;
    for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
    return r;
}

// 用质因数分解独立计算 C(n,k) mod p（不依赖阶乘逆元）
ll brute_comb_mod(int n, int k, ll p) {
    if (k < 0 || k > n) return 0;
    vector<int> expo(n + 1, 0);
    for (int i = n - k + 1; i <= n; i++) {
        int t = i;
        for (int q = 2; q * q <= t; q++) while (t % q == 0) { expo[q]++; t /= q; }
        if (t > 1) expo[t]++;
    }
    for (int i = 1; i <= k; i++) {
        int t = i;
        for (int q = 2; q * q <= t; q++) while (t % q == 0) { expo[q]--; t /= q; }
        if (t > 1) expo[t]--;
    }
    ll r = 1 % p;
    for (int q = 2; q <= n; q++) {
        if (!expo[q]) continue;
        if (q % p == 0) return 0;               // 结果是 p 的倍数
        r = r * qpow(q % p, expo[q], p) % p;
    }
    return r;
}

// 卡特兰
vector<ll> catalan(int n, ll MOD) {
    vector<ll> C(n + 1), f(2 * n + 2), g(2 * n + 2);
    f[0] = 1;
    for (int i = 1; i <= 2 * n + 1; i++) f[i] = f[i - 1] * i % MOD;
    g[2 * n + 1] = qpow(f[2 * n + 1], MOD - 2, MOD);
    for (int i = 2 * n + 1; i >= 1; i--) g[i - 1] = g[i] * i % MOD;
    auto comb = [&](int a, int b) -> ll {
        if (b < 0 || b > a) return 0;
        return f[a] * g[b] % MOD * g[a - b] % MOD;
    };
    for (int i = 0; i <= n; i++) C[i] = comb(2 * i, i) * qpow(i + 1, MOD - 2, MOD) % MOD;
    return C;
}

// 第二类斯特林
vector<vector<ll>> stirling2(int n, ll MOD) {
    vector<vector<ll>> S(n + 1, vector<ll>(n + 1, 0));
    S[0][0] = 1 % MOD;
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= i; j++)
            S[i][j] = (S[i - 1][j - 1] + j % MOD * S[i - 1][j]) % MOD;
    return S;
}

// 整数拆分（五边形数定理）
vector<ll> partition_numbers(int n, ll MOD) {
    vector<ll> p(n + 1, 0);
    p[0] = 1 % MOD;
    for (int i = 1; i <= n; i++) {
        ll s = 0;
        for (int k = 1; ; k++) {
            ll g1 = (ll)k * (3 * k - 1) / 2;
            ll g2 = (ll)k * (3 * k + 1) / 2;
            if (g1 > i && g2 > i) break;
            ll sign = (k & 1) ? 1 : -1;
            if (g1 <= i) s += sign * p[i - g1];
            if (g2 <= i) s += sign * p[i - g2];
        }
        p[i] = ((s % MOD) + MOD) % MOD;
    }
    return p;
}

// 错排
vector<ll> derangement(int n, ll MOD) {
    vector<ll> D(n + 1);
    D[0] = 1 % MOD;
    if (n >= 1) D[1] = 0;
    for (int i = 2; i <= n; i++) D[i] = (i - 1) % MOD * ((D[i - 1] + D[i - 2]) % MOD) % MOD;
    return D;
}

// ==================================================================
// 第 9 章：计算几何
// ==================================================================
const double EPS = 1e-9;
int sgn(double x) { return x > EPS ? 1 : (x < -EPS ? -1 : 0); }

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    Point operator+(const Point &o) const { return Point(x + o.x, y + o.y); }
    Point operator-(const Point &o) const { return Point(x - o.x, y - o.y); }
    Point operator*(double k) const { return Point(x * k, y * k); }
    Point operator/(double k) const { return Point(x / k, y / k); }
    bool operator<(const Point &o) const {
        return sgn(x - o.x) != 0 ? x < o.x : y < o.y;
    }
    bool operator==(const Point &o) const { return sgn(x - o.x) == 0 && sgn(y - o.y) == 0; }
};
double dot(const Point &a, const Point &b) { return a.x * b.x + a.y * b.y; }
double cross(const Point &a, const Point &b) { return a.x * b.y - a.y * b.x; }
double cross(const Point &a, const Point &b, const Point &c) { return cross(b - a, c - a); }
double len(const Point &a) { return sqrt(dot(a, a)); }
double dist(const Point &a, const Point &b) { return len(a - b); }

bool on_segment(const Point &a, const Point &b, const Point &p) {
    return sgn(cross(a, b, p)) == 0 && sgn(dot(a - p, b - p)) <= 0;
}
bool segment_intersect_strict(const Point &a, const Point &b, const Point &c, const Point &d) {
    int c1 = sgn(cross(a, b, c)), c2 = sgn(cross(a, b, d));
    int c3 = sgn(cross(c, d, a)), c4 = sgn(cross(c, d, b));
    return c1 * c2 < 0 && c3 * c4 < 0;
}
double polygon_area2(const vector<Point> &p) {
    double s = 0;
    int n = (int)p.size();
    for (int i = 0; i < n; i++) s += cross(p[i], p[(i + 1) % n]);
    return s;
}
int point_in_polygon(const vector<Point> &poly, const Point &q) {
    int n = (int)poly.size(), cnt = 0;
    for (int i = 0; i < n; i++) {
        const Point &a = poly[i], &b = poly[(i + 1) % n];
        if (on_segment(a, b, q)) return 2;
        int k = sgn(cross(a, b, q));
        int d1 = sgn(a.y - q.y), d2 = sgn(b.y - q.y);
        if (k > 0 && d1 <= 0 && d2 > 0) ++cnt;
        if (k < 0 && d2 <= 0 && d1 > 0) --cnt;
    }
    return cnt != 0;
}

vector<Point> convex_hull(vector<Point> p) {
    sort(p.begin(), p.end());
    p.erase(unique(p.begin(), p.end()), p.end());
    int n = (int)p.size();
    if (n <= 2) return p;
    vector<Point> h(2 * n);
    int k = 0;
    for (int i = 0; i < n; i++) {
        while (k >= 2 && sgn(cross(h[k - 2], h[k - 1], p[i])) <= 0) k--;
        h[k++] = p[i];
    }
    for (int i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && sgn(cross(h[k - 2], h[k - 1], p[i])) <= 0) k--;
        h[k++] = p[i];
    }
    h.resize(k - 1);
    return h;
}

bool in_convex(const vector<Point> &h, const Point &q) {
    int n = (int)h.size();
    if (n == 1) return dist(h[0], q) < EPS;
    if (n == 2) return on_segment(h[0], h[1], q);
    if (sgn(cross(h[0], h[1], q)) < 0) return false;
    if (sgn(cross(h[0], h[n - 1], q)) > 0) return false;
    int l = 1, r = n - 1;
    while (r - l > 1) {
        int mid = (l + r) >> 1;
        if (sgn(cross(h[0], h[mid], q)) >= 0) l = mid; else r = mid;
    }
    return sgn(cross(h[l], h[l + 1 == n ? 0 : l + 1], q)) >= 0;
}

double convex_diameter(const vector<Point> &h) {
    int n = (int)h.size();
    if (n == 1) return 0;
    if (n == 2) return dist(h[0], h[1]);
    double res = 0;
    int j = 1;
    for (int i = 0; i < n; i++) {
        Point a = h[i], b = h[(i + 1) % n];
        while (fabs(cross(b - a, h[(j + 1) % n] - a)) > fabs(cross(b - a, h[j] - a)))
            j = (j + 1) % n;
        res = max(res, max(dist(a, h[j]), dist(b, h[j])));
    }
    return res;
}

// Pick 定理
pair<ll, ll> pick_theorem(const vector<pair<ll,ll>> &v) {
    int n = (int)v.size();
    ll area2 = 0, B = 0;
    for (int i = 0; i < n; i++) {
        ll x1 = v[i].first, y1 = v[i].second;
        ll x2 = v[(i + 1) % n].first, y2 = v[(i + 1) % n].second;
        area2 += x1 * y2 - x2 * y1;
        B += std::gcd(llabs(x2 - x1), llabs(y2 - y1));
    }
    ll S = llabs(area2) / 2;
    ll I = S - B / 2 + 1;
    return {I, B};
}

// ==================================================================
// 第 8 章：博弈
// ==================================================================
bool sg_take_win(int n, const vector<int> &take) {
    vector<int> g(n + 1, 0);
    for (int i = 1; i <= n; i++) {
        vector<char> mark(n + 2, 0);
        for (int t : take) if (i - t >= 0) mark[g[i - t]] = 1;
        int m = 0;
        while (mark[m]) ++m;
        g[i] = m;
    }
    return g[n] != 0;
}

bool anti_nim(const vector<ll> &a) {
    ll x = 0, cnt1 = 0, big = 0;
    for (ll v : a) {
        x ^= v;
        if (v == 1) ++cnt1;
        if (v > 1) ++big;
    }
    if (big == 0) return (cnt1 % 2 == 0);
    return x != 0;
}

bool staircase_nim(const vector<ll> &stair) {
    ll x = 0;
    for (size_t i = 1; i <= stair.size(); i += 2) x ^= stair[i - 1];
    return x != 0;
}

// ==================================================================
// 第 10 章：单纯形 / 三分 / Burnside
// ==================================================================
struct Simplex {
    static constexpr double SEPS = 1e-9;
    int m, n, tot;                      // m 个约束、n 个决策变量、tot = n + m 个变量
    vector<vector<double>> a;           // (m+1) x (tot+1)，最后一列是 RHS
    vector<int> basic;                  // 每行的基变量编号

    Simplex(const vector<vector<double>> &A, const vector<double> &b, const vector<double> &c) {
        m = (int)b.size();
        n = (int)c.size();
        tot = n + m;
        a.assign(m + 1, vector<double>(tot + 1, 0.0));
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) a[i][j] = A[i][j];
            a[i][n + i] = 1.0;                       // 松弛变量列（必须显式保存！）
            a[i][tot] = b[i];
        }
        for (int j = 0; j < n; j++) a[m][j] = -c[j];
        basic.resize(m);
        for (int i = 0; i < m; i++) basic[i] = n + i;
    }
    void pivot(int r, int c) {
        double k = a[r][c];
        for (int j = 0; j <= tot; j++) a[r][j] /= k;
        for (int i = 0; i <= m; i++) {
            if (i == r || fabs(a[i][c]) < SEPS) continue;
            double f = a[i][c];
            for (int j = 0; j <= tot; j++) a[i][j] -= f * a[r][j];
        }
        basic[r] = c;
    }
    bool solve(double &best) {
        while (true) {
            int c = -1;
            for (int j = 0; j < tot; j++) if (a[m][j] < -SEPS) { c = j; break; }
            if (c < 0) break;
            int r = -1;
            for (int i = 0; i < m; i++) {
                if (a[i][c] <= SEPS) continue;
                if (r < 0 || a[i][tot] / a[i][c] < a[r][tot] / a[r][c] - SEPS) r = i;
            }
            if (r < 0) return false;
            pivot(r, c);
        }
        best = a[m][tot];
        return true;
    }
    vector<double> solution() const {
        vector<double> x(n, 0.0);
        for (int i = 0; i < m; i++) if (basic[i] < n) x[basic[i]] = a[i][tot];
        return x;
    }
};

double ternary_search(double l, double r, const function<double(double)> &f) {
    for (int it = 0; it < 200; it++) {
        double m1 = l + (r - l) / 3.0;
        double m2 = r - (r - l) / 3.0;
        if (f(m1) < f(m2)) l = m1; else r = m2;
    }
    return (l + r) / 2;
}

ll necklace_count(ll n, ll k, ll MOD, bool rotationOnly = false) {
    ll sum = 0;
    for (ll i = 0; i < n; i++)
        sum = (sum + qpow(k % MOD, std::gcd(i, n), MOD)) % MOD;
    if (rotationOnly) return sum % MOD * qpow(n % MOD, MOD - 2, MOD) % MOD;
    if (n & 1) {
        sum = (sum + (n % MOD) * qpow(k % MOD, (n + 1) / 2, MOD)) % MOD;
    } else {
        ll half = n / 2 % MOD;
        sum = (sum + half * qpow(k % MOD, n / 2 + 1, MOD)) % MOD;
        sum = (sum + half * qpow(k % MOD, n / 2, MOD)) % MOD;
    }
    return sum % MOD * qpow(2 * n % MOD, MOD - 2, MOD) % MOD;
}

// ==================================================================
int main() {
    printf("=== 1. 组合数 / Lucas / exLucas ===\n");
    {
        const ll P = 1000000007;

        // 参照 1：i128 精确组合数（n <= 60，不会溢出）
        init_comb(60, P);
        bool ok0 = true;
        for (int it = 0; it < 200; it++) {
            int n = (int)(rng() % 61), k = (int)(rng() % 61);
            if (Csmall(n, k) != (ll)(bigC(n, k) % P)) { ok0 = false; printf("    C(%d,%d)\n", n, k); break; }
        }
        check(ok0, "组合数与 i128 精确值一致（n <= 60）");

        // 参照 2：帕斯卡三角 mod P
        const int PN = 1000;
        vector<vector<ll>> pas(PN + 1, vector<ll>(PN + 1, 0));
        for (int i = 0; i <= PN; i++) {
            pas[i][0] = 1;
            for (int j = 1; j <= i; j++)
                pas[i][j] = (pas[i - 1][j - 1] + pas[i - 1][j]) % P;
        }
        init_comb(PN, P);
        bool ok = true;
        for (int it = 0; it < 500; it++) {
            int n = (int)(rng() % (PN + 1)), k = (int)(rng() % (PN + 1));
            if (Csmall(n, k) != pas[n][k]) { ok = false; printf("    C(%d,%d)\n", n, k); break; }
        }
        check(ok, "阶乘+逆元组合数与帕斯卡三角一致（n <= 1000）");

        // Lucas：与质因数分解法算出的 C(n,k) mod 997 对比
        const ll LP = 997;
        init_comb((int)LP - 1, LP);
        bool ok2 = true;
        for (int it = 0; it < 300; it++) {
            int n = (int)(rng() % 3000), k = (int)(rng() % 3000);
            if (k > n) swap(n, k);
            ll expect = brute_comb_mod(n, k, LP);
            ll got = lucas(n, k, LP);
            if (got != expect) { ok2 = false; printf("    C(%d,%d) 得到 %lld 期望 %lld\n", n, k, got, expect); break; }
        }
        check(ok2, "Lucas 与质因数分解法一致（n <= 3000, p = 997）");

        // exLucas：模数任意
        bool ok3 = true;
        ll mods[] = {12, 100, 1000, 999983, 1000000000LL, 97 * 89, 2 * 3 * 5 * 7};
        for (ll mod : mods) {
            for (int it = 0; it < 60; it++) {
                int n = (int)(rng() % 60), k = (int)(rng() % 60);
                ll expect = (ll)(bigC(n, k) % mod);
                ll got = exlucas(n, k, mod);
                if (got != expect) { ok3 = false; printf("    C(%d,%d) mod %lld 得到 %lld 期望 %lld\n", n, k, mod, got, expect); break; }
            }
            if (!ok3) break;
        }
        check(ok3, "exLucas 与精确组合数取模一致（含合数模）");
    }

    printf("=== 2. 卡特兰 / 斯特林 / 拆分数 / 错排 ===\n");
    {
        const ll MOD = 1000000007;
        // 卡特兰：与递推 C_n = Σ C_i C_{n-1-i} 对比
        {
            vector<ll> C = catalan(60, MOD);
            vector<ll> D(61, 0);
            D[0] = 1;
            for (int n = 1; n <= 60; n++)
                for (int i = 0; i < n; i++)
                    D[n] = (D[n] + D[i] * D[n - 1 - i]) % MOD;
            check(C == D, "卡特兰闭式与递推一致（n ≤ 60）");
        }
        // 第二类斯特林：与暴力集合划分对比
        {
            vector<vector<ll>> S = stirling2(7, MOD);
            auto brute = [&](int n, int k) -> ll {
                // 枚举每个元素放到哪个盒子，最后除以 k!
                ll cnt = 0;
                vector<int> assign(n, 0);
                ll total = 1;
                for (int i = 0; i < n; i++) total *= k;
                for (ll mask = 0; mask < total; mask++) {
                    ll t = mask;
                    vector<int> used(k, 0);
                    for (int i = 0; i < n; i++) { assign[i] = (int)(t % k); t /= k; used[assign[i]] = 1; }
                    bool full = true;
                    for (int i = 0; i < k; i++) if (!used[i]) full = false;
                    if (full) ++cnt;
                }
                ll f = 1;
                for (int i = 1; i <= k; i++) f *= i;
                return cnt / f;
            };
            bool ok = true;
            for (int n = 1; n <= 6; n++)
                for (int k = 1; k <= n; k++)
                    if (S[n][k] != brute(n, k) % MOD) { ok = false; }
            check(ok, "第二类斯特林数与暴力划分一致（n ≤ 6）");
        }
        // 整数拆分：与暴力 DP 对比
        {
            const int N = 200;
            vector<ll> p = partition_numbers(N, MOD);
            vector<vector<ll>> dp(N + 1, vector<ll>(N + 1, 0));
            dp[0][0] = 1;
            for (int i = 1; i <= N; i++)
                for (int j = 1; j <= i; j++)
                    dp[i][j] = (dp[i - 1][j - 1] + dp[i - j][j]) % MOD;
            bool ok = true;
            for (int i = 0; i <= N; i++) {
                ll s = 0;
                for (int j = 0; j <= i; j++) s = (s + dp[i][j]) % MOD;
                if (p[i] != s) { ok = false; printf("    p[%d]=%lld 期望 %lld\n", i, p[i], s); break; }
            }
            check(ok, "五边形数定理求拆分数与 DP 一致（n ≤ 200）");
        }
        // 错排：与暴力枚举排列对比
        {
            vector<ll> D = derangement(9, MOD);
            bool ok = true;
            for (int n = 1; n <= 8; n++) {
                vector<int> perm(n);
                for (int i = 0; i < n; i++) perm[i] = i;
                ll cnt = 0;
                do {
                    bool fixed = false;
                    for (int i = 0; i < n; i++) if (perm[i] == i) fixed = true;
                    if (!fixed) ++cnt;
                } while (next_permutation(perm.begin(), perm.end()));
                if (D[n] != cnt % MOD) { ok = false; printf("    D[%d]=%lld 期望 %lld\n", n, D[n], cnt); break; }
            }
            check(ok, "错排递推与暴力枚举排列一致（n ≤ 8）");
        }
    }

    printf("=== 3. 二项式反演 ===\n");
    {
        const ll MOD = 1000000007;
        init_comb(200, MOD);
        bool ok = true;
        for (int it = 0; it < 30; it++) {
            int n = (int)(rng() % 30) + 1;
            vector<ll> g(n + 1);
            for (ll &v : g) v = (ll)(rng() % MOD);
            vector<ll> f(n + 1, 0);
            for (int k = 0; k <= n; k++)
                for (int i = k; i <= n; i++)
                    f[k] = (f[k] + Csmall(i, k) * g[i]) % MOD;
            for (int k = 0; k <= n; k++) {
                ll s = 0;
                for (int i = k; i <= n; i++) {
                    ll term = Csmall(i, k) * f[i] % MOD;
                    if ((i - k) & 1) s = (s - term + MOD) % MOD;
                    else s = (s + term) % MOD;
                }
                if (s != g[k]) { ok = false; break; }
            }
            if (!ok) break;
        }
        check(ok, "二项式反演 g(k) = Σ(-1)^{i-k}C(i,k)f(i) 能还原 g");
    }

    printf("=== 4. 凸包与旋转卡壳 ===\n");
    {
        bool ok = true, ok2 = true, ok3 = true;
        for (int it = 0; it < 60; it++) {
            int n = (int)(rng() % 40) + 3;
            vector<Point> p(n);
            for (auto &q : p) q = Point((double)(rng() % 100), (double)(rng() % 100));
            vector<Point> h = convex_hull(p);
            if ((int)h.size() < 3) continue;
            // 1) 所有点都在凸包内
            for (auto &q : p) if (!in_convex(h, q)) { ok = false; break; }
            // 2) 凸包顶点都是原集合中的点
            for (auto &q : h) {
                bool found = false;
                for (auto &r : p) if (sgn(dist(q, r)) == 0) found = true;
                if (!found) ok2 = false;
            }
            // 3) 直径与暴力最远点对一致
            double bd = 0;
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++) bd = max(bd, dist(p[i], p[j]));
            if (fabs(convex_diameter(h) - bd) > 1e-6) { ok3 = false; }
            if (!ok || !ok2 || !ok3) break;
        }
        check(ok, "凸包包含所有输入点");
        check(ok2, "凸包顶点均为输入点");
        check(ok3, "旋转卡壳求直径 = 暴力最远点对");
    }

    printf("=== 5. 点在多边形内 / Pick 定理 ===\n");
    {
        // 正方形 (0,0)-(10,0)-(10,10)-(0,10)
        vector<Point> sq = {Point(0,0), Point(10,0), Point(10,10), Point(0,10)};
        bool ok = true;
        if (point_in_polygon(sq, Point(5,5)) != 1) ok = false;
        if (point_in_polygon(sq, Point(15,5)) != 0) ok = false;
        if (point_in_polygon(sq, Point(0,5)) != 2) ok = false;
        if (point_in_polygon(sq, Point(10,10)) != 2) ok = false;
        // 三角形
        vector<Point> tri = {Point(0,0), Point(10,0), Point(0,10)};
        if (point_in_polygon(tri, Point(1,1)) != 1) ok = false;
        if (point_in_polygon(tri, Point(6,6)) != 0) ok = false;
        if (point_in_polygon(tri, Point(5,0)) != 2) ok = false;
        check(ok, "射线法判断点在多边形内/边界");

        // Pick 定理：与暴力数格点对比
        bool ok2 = true;
        for (int it = 0; it < 20; it++) {
            // 随机凸多边形（保证简单）
            vector<pair<ll,ll>> pts;
            int n = (int)(rng() % 5) + 3;
            for (int i = 0; i < n; i++) pts.push_back({(ll)(rng() % 12), (ll)(rng() % 12)});
            // 用凸包排序成简单多边形
            sort(pts.begin(), pts.end());
            pts.erase(unique(pts.begin(), pts.end()), pts.end());
            if ((int)pts.size() < 3) continue;
            // 构造凸包（整点）
            vector<pair<ll,ll>> hull;
            int m = (int)pts.size();
            auto crs = [](pair<ll,ll> O, pair<ll,ll> A, pair<ll,ll> B) {
                return (A.first - O.first) * (B.second - O.second) - (A.second - O.second) * (B.first - O.first);
            };
            vector<pair<ll,ll>> hp(2 * m);
            int k = 0;
            for (int i = 0; i < m; i++) {
                while (k >= 2 && crs(hp[k-2], hp[k-1], pts[i]) <= 0) k--;
                hp[k++] = pts[i];
            }
            for (int i = m - 2, t = k + 1; i >= 0; i--) {
                while (k >= t && crs(hp[k-2], hp[k-1], pts[i]) <= 0) k--;
                hp[k++] = pts[i];
            }
            hp.resize(k - 1);
            if ((int)hp.size() < 3) continue;
            auto pr = pick_theorem(hp);
            // 暴力：枚举包围盒内的整点，用射线法判断
            ll minx = 100, maxx = -100, miny = 100, maxy = -100;
            for (auto &q : hp) {
                minx = min(minx, q.first); maxx = max(maxx, q.first);
                miny = min(miny, q.second); maxy = max(maxy, q.second);
            }
            vector<Point> poly;
            for (auto &q : hp) poly.push_back(Point((double)q.first, (double)q.second));
            ll inner = 0, bound = 0;
            for (ll x = minx; x <= maxx; x++)
                for (ll y = miny; y <= maxy; y++) {
                    int r = point_in_polygon(poly, Point((double)x, (double)y));
                    if (r == 1) ++inner;
                    else if (r == 2) ++bound;
                }
            if (pr.first != inner || pr.second != bound) {
                ok2 = false;
                printf("    得到 (I=%lld,B=%lld) 暴力 (I=%lld,B=%lld)\n", pr.first, pr.second, inner, bound);
                break;
            }
        }
        check(ok2, "Pick 定理与暴力数格点一致（S = I + B/2 - 1）");
    }

    printf("=== 6. 博弈：SG / 反 Nim / 阶梯 Nim ===\n");
    {
        // sg_take：与直接 DP 胜负对比（两者互相独立地实现）
        bool ok = true;
        for (int it = 0; it < 60; it++) {
            int k = (int)(rng() % 4) + 1;
            vector<int> take(k);
            for (int &t : take) t = (int)(rng() % 5) + 1;
            int n = (int)(rng() % 60);
            // 直接 DP：win[i] = 存在 t 使 !win[i-t]
            vector<bool> win(n + 1, false);
            for (int i = 1; i <= n; i++) {
                bool w = false;
                for (int t : take) if (i - t >= 0 && !win[i - t]) w = true;
                win[i] = w;
            }
            if (sg_take_win(n, take) != win[n]) { ok = false; break; }
        }
        check(ok, "SG 判定与直接胜负 DP 一致（取石子游戏）");

        // 反 Nim：与暴力博弈树对比（堆数 ≤ 3，堆大小 ≤ 5）
        bool ok2 = true;
        function<bool(vector<int>&)> antiWin = [&](vector<int> &st) -> bool {
            // 返回 true 表示当前局面（轮到当前玩家）必胜
            bool anyMove = false;
            int m = (int)st.size();
            for (int i = 0; i < m; i++) {
                for (int take = 1; take <= st[i]; take++) {
                    anyMove = true;
                    st[i] -= take;
                    bool opp = antiWin(st);
                    st[i] += take;
                    if (!opp) return true;      // 走到对手必败态
                }
            }
            if (!anyMove) return true;           // 无法行动 => 取最后一个者输 => 当前玩家赢
            return false;
        };
        for (int it = 0; it < 120; it++) {
            int m = (int)(rng() % 3) + 1;
            vector<int> st(m);
            for (int &v : st) v = (int)(rng() % 6);
            vector<ll> a(st.begin(), st.end());
            vector<int> st2 = st;
            bool brute = antiWin(st2);
            if (anti_nim(a) != brute) {
                ok2 = false;
                printf("    局面 (");
                for (int v : st) printf("%d,", v);
                printf(") 公式=%d 暴力=%d\n", (int)anti_nim(a), (int)brute);
                break;
            }
        }
        check(ok2, "反 Nim 公式与暴力博弈树一致");

        // 阶梯 Nim：与暴力博弈树对比（级数 ≤ 3，每级 ≤ 3）
        bool ok3 = true;
        function<bool(vector<int>&)> stairWin = [&](vector<int> &st) -> bool {
            bool anyMove = false;
            int levels = (int)st.size();
            for (int i = 1; i < levels; i++) {
                for (int take = 1; take <= st[i]; take++) {
                    anyMove = true;
                    st[i] -= take;
                    st[i - 1] += take;
                    bool opp = stairWin(st);
                    st[i] += take;
                    st[i - 1] -= take;
                    if (!opp) return true;
                }
            }
            if (!anyMove) return false;          // 无石子可动 => 输
            return false;
        };
        for (int it = 0; it < 120; it++) {
            int levels = (int)(rng() % 3) + 2;
            vector<int> st(levels, 0);
            for (int i = 1; i < levels; i++) st[i] = (int)(rng() % 4);
            vector<ll> stair(st.begin() + 1, st.end());   // 第 1 级对应 stair[0]
            vector<int> st2 = st;
            bool brute = stairWin(st2);
            if (staircase_nim(stair) != brute) {
                ok3 = false;
                printf("    阶梯 (");
                for (int v : st) printf("%d,", v);
                printf(") 公式=%d 暴力=%d\n", (int)staircase_nim(stair), (int)brute);
                break;
            }
        }
        check(ok3, "阶梯 Nim 公式与暴力博弈树一致");
    }

    printf("=== 7. 单纯形法 ===\n");
    {
        // 与"枚举约束交点+可行性检查"的暴力对比
        auto brute_lp = [](const vector<vector<double>> &A0, const vector<double> &b0,
                           const vector<double> &c) -> double {
            int n = (int)c.size();
            // 把 x >= 0 也作为约束加入（-x_j <= 0），否则会漏掉位于坐标平面上的顶点
            vector<vector<double>> A = A0;
            vector<double> b = b0;
            for (int j = 0; j < n; j++) {
                vector<double> row(n, 0.0);
                row[j] = -1.0;
                A.push_back(row);
                b.push_back(0.0);
            }
            int m = (int)b.size();
            double best = -1e18;
            // 枚举 n 个约束的组合，解线性方程组（用 Cramer 法则）
            vector<int> idx(m);
            for (int i = 0; i < m; i++) idx[i] = i;
            vector<int> sel(n);
            function<void(int,int)> rec = [&](int start, int cnt) {
                if (cnt == n) {
                    // 解 A_sel x = b_sel
                    vector<vector<double>> M(n, vector<double>(n));
                    vector<double> rhs(n);
                    for (int i = 0; i < n; i++) {
                        M[i] = A[sel[i]];
                        rhs[i] = b[sel[i]];
                    }
                    double det = 0;
                    // 用高斯消元求 det 与解
                    vector<vector<double>> T = M;
                    vector<double> r = rhs;
                    double d = 1;
                    bool ok = true;
                    for (int i = 0; i < n; i++) {
                        int piv = i;
                        for (int j = i; j < n; j++) if (fabs(T[j][i]) > fabs(T[piv][i])) piv = j;
                        if (fabs(T[piv][i]) < 1e-9) { ok = false; break; }
                        swap(T[piv], T[i]); swap(r[piv], r[i]);
                        if (piv != i) d = -d;
                        for (int j = i + 1; j < n; j++) {
                            double f = T[j][i] / T[i][i];
                            for (int k = i; k < n; k++) T[j][k] -= f * T[i][k];
                            r[j] -= f * r[i];
                        }
                    }
                    det = d;
                    for (int i = 0; i < n; i++) det *= T[i][i];
                    if (!ok || fabs(det) < 1e-9) return;
                    // 回代
                    vector<double> x(n);
                    for (int i = n - 1; i >= 0; i--) {
                        double s = r[i];
                        for (int j = i + 1; j < n; j++) s -= T[i][j] * x[j];
                        x[i] = s / T[i][i];
                    }
                    for (int i = 0; i < n; i++) if (x[i] < -1e-7) return;
                    for (int i = 0; i < m; i++) {
                        double s = 0;
                        for (int j = 0; j < n; j++) s += A[i][j] * x[j];
                        if (s > b[i] + 1e-7) return;
                    }
                    double val = 0;
                    for (int j = 0; j < n; j++) val += c[j] * x[j];
                    best = max(best, val);
                    return;
                }
                for (int i = start; i < m; i++) {
                    sel[cnt] = i;
                    rec(i + 1, cnt + 1);
                }
            };
            rec(0, 0);
            return best;
        };

        bool ok = true;
        int tested = 0;
        for (int it = 0; it < 400 && tested < 40; it++) {
            int n = (int)(rng() % 3) + 1;
            int m = n + 1 + (int)(rng() % 2);
            vector<vector<double>> A(m, vector<double>(n));
            vector<double> b(m);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) A[i][j] = (double)((int)(rng() % 7) - 3);
                b[i] = (double)(rng() % 20) + 1;
            }
            vector<double> c(n);
            for (int j = 0; j < n; j++) c[j] = (double)((int)(rng() % 7) - 3);

            double expect = brute_lp(A, b, c);
            if (expect < -1e17) continue;          // 无可行顶点，跳过

            Simplex sp(A, b, c);
            double got;
            bool bounded = sp.solve(got);
            if (!bounded) continue;
            // 目标有界且暴力找到了最优顶点
            tested++;
            if (fabs(got - expect) > 1e-6) {
                ok = false;
                printf("    n=%d m=%d 单纯形 %.6f 暴力 %.6f\n", n, m, got, expect);
                break;
            }
        }
        check(ok && tested > 0, "单纯形法最优值与枚举顶点暴力一致");
        printf("    （对比了 %d 个随机线性规划）\n", tested);
    }

    printf("=== 8. 三分法 ===\n");
    {
        bool ok = true;
        for (int it = 0; it < 20; it++) {
            double a = (double)((int)(rng() % 21) - 10);
            double b = (double)((int)(rng() % 11) - 5);
            double c = (double)(rng() % 10);
            // 单峰（开口向下）函数：极大值点在 x = a + b/2
            auto f = [&](double x) { return -(x - a) * (x - a) + b * (x - a) + c; };
            double x = ternary_search(-50, 50, f);
            double opt = a + b / 2;
            if (fabs(x - opt) > 1e-4) { ok = false; printf("    得到 %.6f 期望 %.6f\n", x, opt); break; }
        }
        check(ok, "三分法求得单峰函数极大值点");
    }

    printf("=== 9. Burnside / Polya 项链计数 ===\n");
    {
        const ll MOD = 1000000007;
        // 暴力枚举所有 k^n 染色，用旋转/翻转生成等价类
        auto brute_necklace = [](int n, int k) -> ll {
            int total = 1;
            for (int i = 0; i < n; i++) total *= k;
            vector<int> fa(total);
            for (int i = 0; i < total; i++) fa[i] = i;
            function<int(int)> find = [&](int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); };
            auto unite = [&](int x, int y) { fa[find(x)] = find(y); };
            auto decode = [&](int code, vector<int> &v) {
                for (int i = 0; i < n; i++) { v[i] = code % k; code /= k; }
            };
            auto encode = [&](const vector<int> &v) {
                int code = 0, p = 1;
                for (int i = 0; i < n; i++) { code += v[i] * p; p *= k; }
                return code;
            };
            vector<int> v(n), w(n);
            for (int code = 0; code < total; code++) {
                decode(code, v);
                for (int s = 0; s < n; s++) {          // 旋转
                    for (int i = 0; i < n; i++) w[i] = v[(i + s) % n];
                    unite(code, encode(w));
                }
                for (int s = 0; s < n; s++) {          // 翻转
                    for (int i = 0; i < n; i++) w[i] = v[(s - i % n + n) % n];
                    unite(code, encode(w));
                }
            }
            ll cnt = 0;
            for (int i = 0; i < total; i++) if (find(i) == i) ++cnt;
            return cnt;
        };
        bool ok = true;
        for (int n = 1; n <= 6; n++) {
            for (int k = 1; k <= 3; k++) {
                if (n * k > 20) continue;
                ll b = brute_necklace(n, k);
                ll f = necklace_count(n, k, MOD, false);
                if (b % MOD != f) { ok = false; printf("    n=%d k=%d 公式=%lld 暴力=%lld\n", n, k, f, b); break; }
            }
            if (!ok) break;
        }
        check(ok, "项链计数（旋转+翻转）与暴力枚举轨道一致");

        bool ok2 = true;
        for (int n = 1; n <= 6; n++) {
            for (int k = 1; k <= 3; k++) {
                if (n * k > 20) continue;
                ll f = necklace_count(n, k, MOD, true);   // 只考虑旋转
                // 暴力只允许旋转
                int total = 1;
                for (int i = 0; i < n; i++) total *= k;
                vector<int> fa(total);
                for (int i = 0; i < total; i++) fa[i] = i;
                function<int(int)> find = [&](int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); };
                vector<int> v(n), w(n);
                for (int code = 0; code < total; code++) {
                    int t = code;
                    for (int i = 0; i < n; i++) { v[i] = t % k; t /= k; }
                    for (int s = 0; s < n; s++) {
                        for (int i = 0; i < n; i++) w[i] = v[(i + s) % n];
                        int c2 = 0, p = 1;
                        for (int i = 0; i < n; i++) { c2 += w[i] * p; p *= k; }
                        fa[find(code)] = find(c2);
                    }
                }
                ll cnt = 0;
                for (int i = 0; i < total; i++) if (find(i) == i) ++cnt;
                if (cnt % MOD != f) { ok2 = false; printf("    n=%d k=%d 公式=%lld 暴力=%lld\n", n, k, f, cnt); break; }
            }
            if (!ok2) break;
        }
        check(ok2, "项链计数（仅旋转）与暴力枚举轨道一致");
    }

    printf("\n===== 通过 %d 项，失败 %d 项 =====\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
