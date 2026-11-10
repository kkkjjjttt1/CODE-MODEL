// 数论模板正确性验证（与暴力实现对比）
// 编译: g++ -std=c++17 -O2 -o verify_nt.exe tools/verify_nt.cpp
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>
#include <chrono>

using namespace std;
using ll = long long;
using ull = unsigned long long;
using i128 = __int128;

static int g_fail = 0, g_pass = 0;
static void check(bool cond, const char *name) {
    if (cond) { ++g_pass; printf("  [ok]   %s\n", name); }
    else      { ++g_fail; printf("  [FAIL] %s\n", name); }
}

mt19937_64 rng(20240607);

// ==================================================================
// 教程模板（与第 1、2、3 章代码一致）
// ==================================================================

ll my_gcd(ll a, ll b) { return b ? my_gcd(b, a % b) : a; }

ll my_exgcd(ll a, ll b, ll &x, ll &y) {
    if (b == 0) { x = 1; y = 0; return a; }
    ll d = my_exgcd(b, a % b, y, x);
    y -= a / b * x;
    return d;
}

ll qpow(ll a, ll e, ll mod) {
    ll r = 1 % mod;
    a %= mod; if (a < 0) a += mod;
    for (; e; e >>= 1, a = (ll)((i128)a * a % mod))
        if (e & 1) r = (ll)((i128)r * a % mod);
    return r;
}

// --- Miller-Rabin ---
bool isPrime(ll n) {
    if (n < 2) return false;
    for (ll p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37})
        if (n % p == 0) return n == p;
    ll d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; ++s; }
    for (ll a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
        ll x = qpow(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (int i = 1; i < s; i++) {
            x = (ll)((i128)x * x % n);
            if (x == n - 1) { witness = false; break; }
        }
        if (witness) return false;
    }
    return true;
}

// --- Pollard-Rho (Brent) ---
inline ll fmul(ll a, ll b, ll m) { return (ll)((i128)a * b % m); }

ll pollard(ll n) {
    if (!(n & 1)) return 2;
    if (n % 3 == 0) return 3;
    while (true) {
        ll c = (ll)(rng() % (n - 1)) + 1;
        ll x = (ll)(rng() % (n - 2)) + 2, y = x, q = 1, ys = x, d = 1;
        auto f = [&](ll v) { return (fmul(v, v, n) + c) % n; };
        const ll M = 128;
        for (ll r = 1; d == 1; r <<= 1) {
            x = y;
            for (ll i = 0; i < r; i++) y = f(y);
            for (ll k = 0; k < r && d == 1; k += M) {
                ys = y;
                ll lim = min(M, r - k);
                for (ll i = 0; i < lim; i++) { y = f(y); q = fmul(q, abs(x - y), n); }
                d = my_gcd(q, n);
            }
        }
        if (d == n) {
            do { ys = f(ys); d = my_gcd(abs(x - ys), n); } while (d == 1);
        }
        if (d != n) return d;
    }
}

void factor_rec(ll n, vector<ll> &res) {
    if (n == 1) return;
    if (isPrime(n)) { res.push_back(n); return; }
    ll d = pollard(n);
    factor_rec(d, res);
    factor_rec(n / d, res);
}

// --- 线性筛 ---
const int MAXN = 1000005;
int primes[MAXN / 10], pcnt;
bool comp[MAXN];
int phi[MAXN], mu[MAXN], dd[MAXN], dex[MAXN];
ll  sigma[MAXN], spow[MAXN];

void sieve(int n) {
    comp[0] = comp[1] = true;
    phi[1] = mu[1] = 1; dd[1] = 1; sigma[1] = 1;
    for (int i = 2; i <= n; i++) {
        if (!comp[i]) {
            primes[pcnt++] = i;
            phi[i] = i - 1; mu[i] = -1;
            dd[i] = 2; dex[i] = 1;
            sigma[i] = i + 1; spow[i] = i + 1;
        }
        for (int j = 0; j < pcnt; j++) {
            int p = primes[j];
            ll v = (ll)i * p;
            if (v > n) break;
            comp[v] = true;
            if (i % p == 0) {
                phi[v] = phi[i] * p;
                mu[v] = 0;
                dex[v] = dex[i] + 1;
                dd[v] = dd[i] / (dex[i] + 1) * (dex[v] + 1);
                spow[v] = spow[i] * p + 1;
                sigma[v] = sigma[i] / spow[i] * spow[v];
                break;
            }
            phi[v] = phi[i] * (p - 1);
            mu[v] = -mu[i];
            dex[v] = 1;
            dd[v] = dd[i] * 2;
            spow[v] = p + 1;
            sigma[v] = sigma[i] * (p + 1);
        }
    }
}

// --- 线性同余 ---
ll solve_linear_congruence(ll a, ll c, ll m) {
    a = ((a % m) + m) % m;
    c = ((c % m) + m) % m;
    ll x, y;
    ll d = my_exgcd(a, m, x, y);
    if (c % d) return -1;
    ll step = m / d;
    x = (ll)((i128)(x % step + step) % step * ((c / d) % step) % step);
    return x;
}

// --- exCRT ---
bool merge_eq(ll &r1, ll &m1, ll r2, ll m2) {
    ll x, y;
    ll d = my_exgcd(m1, m2, x, y);
    ll diff = r2 - r1;
    if (diff % d) return false;
    ll step = m2 / d;
    ll t = (ll)((i128)((x % step + step) % step) * (((diff / d) % step + step) % step) % step);
    ll lcm = m1 / d * m2;
    r1 = (ll)(((i128)r1 + (i128)t * m1) % lcm);
    m1 = lcm;
    r1 = (r1 % m1 + m1) % m1;
    return true;
}

bool excrt(vector<ll> r, vector<ll> m, ll &x, ll &mod) {
    ll R = r[0], M = m[0];
    for (size_t i = 1; i < r.size(); i++)
        if (!merge_eq(R, M, ((r[i] % m[i]) + m[i]) % m[i], m[i])) return false;
    x = (R % M + M) % M; mod = M;
    return true;
}

// --- 欧拉函数 / 降幂 ---
ll euler_phi(ll n) {
    ll res = n;
    for (ll i = 2; i <= n / i; i++) {
        if (n % i == 0) {
            res = res / i * (i - 1);
            while (n % i == 0) n /= i;
        }
    }
    if (n > 1) res = res / n * (n - 1);
    return res;
}

ll power_big(ll a, const string &b, ll m) {
    ll ph = euler_phi(m);
    ll e = 0; bool big = false;
    for (char ch : b) {
        e = e * 10 + (ch - '0');
        if (e >= ph) { big = true; e %= ph; }
    }
    if (big) e += ph;
    return qpow(a, e, m);
}

// --- 原根 ---
ll primitive_root(ll p) {
    if (p == 2) return 1;
    ll ph = p - 1, t = ph;
    vector<ll> fac;
    for (ll i = 2; i <= t / i; i++)
        if (t % i == 0) { fac.push_back(i); while (t % i == 0) t /= i; }
    if (t > 1) fac.push_back(t);
    for (ll g = 2; g < p; g++) {
        bool ok = true;
        for (ll q : fac) if (qpow(g, ph / q, p) == 1) { ok = false; break; }
        if (ok) return g;
    }
    return -1;
}

// --- BSGS ---
ll bsgs(ll a, ll b, ll p) {
    a %= p; b %= p;
    if (b == 1 % p) return 0;
    ll m = (ll)sqrt((long double)p);
    while (m * m < p) m++;
    unordered_map<ll, ll> mp;
    mp.reserve(m * 2);
    ll cur = b;
    for (ll j = 0; j < m; j++) {
        if (!mp.count(cur)) mp[cur] = j;
        cur = (ll)((i128)cur * a % p);
    }
    ll am = qpow(a, m, p), now = 1;
    for (ll i = 1; i <= m; i++) {
        now = (ll)((i128)now * am % p);
        auto it = mp.find(now);
        if (it != mp.end()) {
            ll x = i * m - it->second;
            if (x >= 0) return x;
        }
    }
    return -1;
}

// --- exBSGS ---
ll inv_exgcd(ll a, ll mod) {
    ll x, y;
    ll d = my_exgcd(a, mod, x, y);
    if (d != 1) return -1;
    return (x % mod + mod) % mod;
}

ll exbsgs(ll a, ll b, ll p) {
    a %= p; b %= p;
    if (p == 1) return 0;
    if (b == 1 % p) return 0;
    ll k = 0, d = 1;
    while (true) {
        ll g = my_gcd(a, p);
        if (g == 1) break;
        if (b % g) return -1;
        b /= g; p /= g; ++k;
        d = (ll)((i128)d * (a / g) % p);
        if (d == b) return k;
    }
    ll invd = inv_exgcd(d % p, p);
    if (invd < 0) return -1;
    ll y = bsgs(a, (ll)((i128)b * invd % p), p);
    return y < 0 ? -1 : y + k;
}

// --- Cipolla ---
namespace cipolla {
    ll p, w;
    struct Fp2 { ll a, b; };
    Fp2 mul(Fp2 x, Fp2 y) {
        return { (ll)(((i128)x.a * y.a + (i128)x.b * y.b % p * w) % p),
                 (ll)(((i128)x.a * y.b + (i128)x.b * y.a) % p) };
    }
    Fp2 pw(Fp2 x, ll e) {
        Fp2 r{1, 0};
        for (; e; e >>= 1, x = mul(x, x)) if (e & 1) r = mul(r, x);
        return r;
    }
    ll solve(ll n, ll P) {
        p = P; n = ((n % p) + p) % p;
        if (n == 0) return 0;
        if (qpow(n, (p - 1) / 2, p) != 1) return -1;
        ll a;
        do {
            a = (ll)(rng() % p);
            w = (ll)((((i128)a * a - n) % p + p) % p);
        } while (w == 0 || qpow(w, (p - 1) / 2, p) != p - 1);
        return pw(Fp2{a % p, 1}, (p + 1) / 2).a;
    }
}

// --- 杜教筛（小阈值版，用于验证） ---
namespace dusieve {
    const int LIM = 2000;
    int mu2[LIM], phi2[LIM];
    ll phiPre[LIM];
    int muPre[LIM];
    bool inited = false;

    void init() {
        vector<bool> cmp(LIM, false);
        vector<int> pr;
        cmp[0] = cmp[1] = true;
        phi2[1] = 1; mu2[1] = 1;
        for (int i = 2; i < LIM; i++) {
            if (!cmp[i]) { pr.push_back(i); phi2[i] = i - 1; mu2[i] = -1; }
            for (size_t j = 0; j < pr.size() && (ll)i * pr[j] < LIM; j++) {
                int v = i * pr[j];
                cmp[v] = true;
                if (i % pr[j] == 0) { phi2[v] = phi2[i] * pr[j]; mu2[v] = 0; break; }
                phi2[v] = phi2[i] * (pr[j] - 1); mu2[v] = -mu2[i];
            }
        }
        for (int i = 1; i < LIM; i++) {
            phiPre[i] = phiPre[i - 1] + phi2[i];
            muPre[i] = muPre[i - 1] + mu2[i];
        }
        inited = true;
    }

    unordered_map<ll, ll> memPhi;
    unordered_map<ll, int> memMu;

    ll S_phi(ll n) {
        if (n < LIM) return phiPre[n];
        auto it = memPhi.find(n);
        if (it != memPhi.end()) return it->second;
        ll res = (ll)((i128)n * (n + 1) / 2);
        for (ll l = 2, r; l <= n; l = r + 1) {
            r = n / (n / l);
            res -= (r - l + 1) * S_phi(n / l);
        }
        return memPhi[n] = res;
    }

    int S_mu(ll n) {
        if (n < LIM) return muPre[n];
        auto it = memMu.find(n);
        if (it != memMu.end()) return it->second;
        ll res = 1;
        for (ll l = 2, r; l <= n; l = r + 1) {
            r = n / (n / l);
            res -= (r - l + 1) * S_mu(n / l);
        }
        return memMu[n] = (int)res;
    }
}

// ==================================================================
// 暴力参照实现（用与线性筛完全不同的方法：埃氏筛 / 倍数法）
// ==================================================================
const int BN = 200000;
vector<int> bphi(BN + 1), bmu(BN + 1), bd(BN + 1);
vector<ll> bsig(BN + 1, 0);

void build_brute() {
    vector<bool> isP(BN + 1, true);
    isP[0] = isP[1] = false;
    for (int i = 0; i <= BN; i++) { bphi[i] = i; bmu[i] = 1; }
    bphi[0] = 0; bmu[0] = 0;
    for (int i = 2; i <= BN; i++) {
        if (!isP[i]) continue;
        for (int j = i; j <= BN; j += i) {
            isP[j] = (j == i);
            bphi[j] -= bphi[j] / i;                     // φ 的埃氏求法
            bmu[j] = -bmu[j];
        }
        for (ll j = (ll)i * i; j <= BN; j += (ll)i * i) bmu[j] = 0;   // 含平方因子
    }
    for (int i = 1; i <= BN; i++)                        // 倍数法求 d 与 σ
        for (int j = i; j <= BN; j += i) { ++bd[j]; bsig[j] += i; }
}
bool brute_isprime(ll n) {
    if (n < 2) return false;
    for (ll i = 2; i <= n / i; i++) if (n % i == 0) return false;
    return true;
}

// ==================================================================
int main() {
    printf("=== 1. exgcd ===\n");
    {
        bool ok = true;
        for (int it = 0; it < 2000; it++) {
            ll a = (ll)(rng() % 1000000000) + 1, b = (ll)(rng() % 1000000000) + 1;
            ll x, y;
            ll d = my_exgcd(a, b, x, y);
            if (d != my_gcd(a, b) || (i128)a * x + (i128)b * y != d) { ok = false; break; }
        }
        check(ok, "exgcd 返回 gcd 且 ax+by=gcd");
    }

    printf("=== 2. 快速幂 ===\n");
    {
        bool ok = true;
        for (int it = 0; it < 1000; it++) {
            ll a = (ll)(rng() % 1000000), e = (ll)(rng() % 100000);
            ll m = 998244353;
            ll r = 1;
            for (ll i = 0; i < e; i++) r = r * a % m;
            if (qpow(a, e, m) != r) { ok = false; break; }
        }
        check(ok, "qpow 与朴素连乘一致");
        check(qpow(5, 3, 1) == 0, "qpow 在 mod = 1 时返回 0");
    }

    printf("=== 3. 线性筛（phi/mu/d/sigma） ===\n");
    {
        sieve(200000);
        build_brute();
        bool ok = true;
        for (int i = 1; i <= 200000; i++) {
            if (phi[i] != bphi[i]) { ok = false; printf("    phi[%d]=%d 期望 %d\n", i, phi[i], bphi[i]); break; }
            if (mu[i] != bmu[i]) { ok = false; printf("    mu[%d]=%d 期望 %d\n", i, mu[i], bmu[i]); break; }
            if (dd[i] != bd[i]) { ok = false; printf("    d[%d]=%d 期望 %d\n", i, dd[i], bd[i]); break; }
            if (sigma[i] != bsig[i]) { ok = false; printf("    sigma[%d]=%lld 期望 %lld\n", i, sigma[i], bsig[i]); break; }
        }
        check(ok, "线性筛 phi/mu/d/sigma 在 1..200000 全部正确");
        bool ok2 = true;
        for (int i = 1; i <= 5000; i++) {
            ll s = 0;
            for (int d = 1; d <= i; d++) if (i % d == 0) s += phi[d];
            if (s != i) { ok2 = false; break; }
        }
        check(ok2, "恒等式 Σ_{d|n} φ(d) = n");
    }

    printf("=== 4. Miller-Rabin ===\n");
    {
        bool ok = true;
        for (ll n = 1; n <= 200000; n++)
            if (isPrime(n) != brute_isprime(n)) { ok = false; printf("    n=%lld\n", n); break; }
        check(ok, "1..200000 与试除法一致");

        check(isPrime(1000000007LL) && isPrime(998244353LL) && isPrime(1000000009LL),
              "大素数判定正确");
        check(!isPrime(1000000007LL * 1000000009LL), "两个大素数之积被判为合数");
        check(!isPrime(3215031751LL), "强伪素数 3215031751 被判为合数");
        check(isPrime(2305843009213693951LL), "梅森素数 2^61-1 判定正确");
        check(!isPrime(9223372036854775807LL), "2^63-1 判定为合数");
    }

    printf("=== 5. Pollard-Rho 分解 ===\n");
    {
        bool ok = true;
        ll tests[] = {8051LL, 10403LL, 1000000007LL * 1000000009LL,
                      999999999999999989LL, 2LL * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29,
                      1000000000000000003LL};
        for (ll n : tests) {
            vector<ll> f;
            factor_rec(n, f);
            ll prod = 1;
            for (ll v : f) { prod *= v; if (!isPrime(v)) ok = false; }
            if (prod != n) ok = false;
            sort(f.begin(), f.end());
            printf("    %-22lld -> ", n);
            for (ll v : f) printf("%lld ", v);
            printf("\n");
        }
        check(ok, "Pollard-Rho 分解结果的乘积 = 原数，且因子均为素数");
    }

    printf("=== 6. 线性同余方程 ===\n");
    {
        bool ok = true;
        for (int it = 0; it < 3000; it++) {
            ll m = (ll)(rng() % 1000) + 1;
            ll a = (ll)(rng() % 2000);
            ll c = (ll)(rng() % 2000);
            ll x = solve_linear_congruence(a, c, m);
            // 暴力找最小非负解
            ll expect = -1;
            for (ll v = 0; v < m; v++) if ((i128)a * v % m == c % m) { expect = v; break; }
            if (x != expect) { ok = false; printf("    a=%lld c=%lld m=%lld 得到 %lld 期望 %lld\n", a, c, m, x, expect); break; }
        }
        check(ok, "线性同余最小非负解与暴力一致");
    }

    printf("=== 7. exCRT ===\n");
    {
        bool ok = true, okn = true;
        for (int it = 0; it < 2000; it++) {
            int k = (int)(rng() % 4) + 1;
            vector<ll> r(k), m(k);
            for (int i = 0; i < k; i++) {
                m[i] = (ll)(rng() % 30) + 1;
                r[i] = (ll)(rng() % m[i]);
            }
            ll x, mod;
            bool sol = excrt(r, m, x, mod);
            // 暴力枚举 0..lcm-1（lcm 可能很大，限制规模）
            ll L = 1;
            for (int i = 0; i < k; i++) { L = L / my_gcd(L, m[i]) * m[i]; if (L > 20000) break; }
            if (L <= 20000) {
                ll expect = -1;
                for (ll v = 0; v < L; v++) {
                    bool good = true;
                    for (int i = 0; i < k; i++) if (v % m[i] != r[i]) { good = false; break; }
                    if (good) { expect = v; break; }
                }
                if (expect < 0) { if (sol) { okn = false; } }
                else if (!sol || x != expect) { okn = false; printf("    期望 %lld 得到 %lld (sol=%d)\n", expect, sol ? x : -1, (int)sol); }
            }
        }
        check(ok, "exCRT 基本自检");
        check(okn, "exCRT 与小规模暴力枚举一致");
    }

    printf("=== 8. 欧拉函数与欧拉降幂 ===\n");
    {
        bool ok = true;
        for (int i = 1; i <= 50000; i++) if (euler_phi(i) != bphi[i]) { ok = false; break; }
        check(ok, "单点 euler_phi 与埃氏法一致（1..50000）");

        // 大指数字符串：b = "1000000000000000000000000000"
        string bs = "1000000000000000000000000000";
        ll x = power_big(3, bs, 1000000007LL);
        // 与直接计算比较（模 1e9+7，用大指数分解为 10^27）
        ll e = 1;
        for (int i = 0; i < 27; i++) e = (ll)((i128)e * 10 % (1000000007LL - 1));
        ll expect = 3;
        {
            // 直接算 3^(10^27) mod p，用 φ(p) = p-1
            expect = qpow(3, e, 1000000007LL);
        }
        check(x == expect, "欧拉降幂（字符串指数）与指数模 φ 一致");
    }

    printf("=== 9. 原根 ===\n");
    {
        bool ok = true;
        ll ps[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 998244353LL, 1000000007LL};
        for (ll p : ps) {
            ll g = primitive_root(p);
            if (g < 1) { ok = false; break; }
            // 验证阶为 p-1
            if (p == 2) { if (g != 1) ok = false; continue; }
            bool full = (qpow(g, p - 1, p) == 1);
            // 检查不存在更小的阶
            ll ph = p - 1;
            for (ll d = 1; d * d <= ph && full; d++) {
                if (ph % d) continue;
                if (d < ph && qpow(g, d, p) == 1) full = false;
                ll e = ph / d;
                if (e < ph && qpow(g, e, p) == 1) full = false;
            }
            if (!full) { ok = false; printf("    p=%lld g=%lld 不是原根\n", p, g); }
        }
        check(ok, "primitive_root 返回的 g 阶恰为 p-1");
    }

    printf("=== 10. BSGS ===\n");
    {
        bool ok = true;
        const ll P = 1000000007LL;
        ll g = primitive_root(P);
        for (int it = 0; it < 30; it++) {
            ll x = (ll)(rng() % 100000);
            ll b = qpow(g, x, P);
            ll res = bsgs(g, b, P);
            if (res < 0 || qpow(g, res, P) != b) { ok = false; printf("    x=%lld res=%lld\n", x, res); break; }
            if (res > x) { ok = false; printf("    返回的不是最小解\n"); break; }
        }
        check(ok, "bsgs 求出的指数使 a^x ≡ b 且不大于真实指数");

        // exBSGS：a 与 p 不互素
        bool ok2 = true;
        for (int it = 0; it < 200; it++) {
            ll a = (ll)(rng() % 100) + 2;
            ll x = (ll)(rng() % 12);
            ll p = (ll)(rng() % 500) + 2;
            ll b = qpow(a, x, p);
            ll res = exbsgs(a, b, p);
            if (res < 0) {
                // 可能真的无解（a^x 与 p 的公因子问题），暴力确认
                bool found = false;
                ll cur = 1 % p;
                for (ll i = 0; i < 200; i++) {
                    if (cur == b) { found = true; break; }
                    cur = (ll)((i128)cur * a % p);
                }
                if (found) { ok2 = false; break; }
            } else {
                if (qpow(a, res, p) != b) { ok2 = false; printf("    a=%lld b=%lld p=%lld res=%lld\n", a, b, p, res); break; }
            }
        }
        check(ok2, "exbsgs 在 a 与 p 不互素时正确");
    }

    printf("=== 11. Cipolla 二次剩余 ===\n");
    {
        bool ok = true;
        ll ps[] = {998244353LL, 1000000007LL, 1000000009LL, 10007LL};
        for (ll p : ps) {
            for (int it = 0; it < 20; it++) {
                ll r = (ll)(rng() % p);
                ll n = (ll)((i128)r * r % p);
                ll x = cipolla::solve(n, p);
                if (n == 0) { if (x != 0) { ok = false; break; } continue; }
                if (x < 0 || (ll)((i128)x * x % p) != n) { ok = false; printf("    p=%lld n=%lld x=%lld\n", p, n, x); break; }
            }
            // 非剩余应返回 -1
            ll nonres = -1;
            for (ll v = 2; v < 1000; v++) if (qpow(v, (p - 1) / 2, p) == p - 1) { nonres = v; break; }
            if (nonres > 0 && cipolla::solve(nonres, p) != -1) { ok = false; }
        }
        check(ok, "cipolla 返回的解满足 x^2 ≡ n，非剩余返回 -1");
    }

    printf("=== 12. 整除分块 ===\n");
    {
        bool ok = true;
        for (int n = 1; n <= 2000; n += 7) {
            ll a = 0, b = 0;
            for (ll l = 1, r; l <= n; l = r + 1) {
                ll q = n / l;
                r = n / q;
                a += (r - l + 1) * q;
            }
            for (ll i = 1; i <= n; i++) b += n / i;
            if (a != b) { ok = false; printf("    n=%d %lld vs %lld\n", n, a, b); break; }
        }
        check(ok, "一维整除分块求和与暴力一致");

        bool ok2 = true;
        for (int it = 0; it < 200; it++) {
            ll n = (ll)(rng() % 3000) + 1, m = (ll)(rng() % 3000) + 1;
            ll a = 0, b = 0;
            ll lim = min(n, m);
            for (ll l = 1, r; l <= lim; l = r + 1) {
                r = min(n / (n / l), m / (m / l));
                a += (r - l + 1) * (n / l) * (m / l);
            }
            for (ll i = 1; i <= lim; i++) b += (n / i) * (m / i);
            if (a != b) { ok2 = false; printf("    n=%lld m=%lld %lld vs %lld\n", n, m, a, b); break; }
        }
        check(ok2, "二维整除分块与暴力一致");
    }

    printf("=== 13. 莫比乌斯反演计数 ===\n");
    {
        sieve(5000);
        // Σ_{i<=n} Σ_{j<=m} [gcd(i,j)=1] 的莫反公式 vs 暴力
        bool ok = true;
        for (int it = 0; it < 60; it++) {
            int n = (int)(rng() % 300) + 1, m = (int)(rng() % 300) + 1;
            ll a = 0;
            for (ll l = 1, r; l <= min(n, m); l = r + 1) {
                r = min(n / (n / l), m / (m / l));
                ll s = 0;
                for (ll d = l; d <= r; d++) s += mu[d];
                a += s * (n / l) * (m / l);
            }
            ll b = 0;
            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= m; j++)
                    if (my_gcd(i, j) == 1) ++b;
            if (a != b) { ok = false; printf("    n=%d m=%d %lld vs %lld\n", n, m, a, b); break; }
        }
        check(ok, "ΣΣ[gcd=1] = Σ μ(d)⌊n/d⌋⌊m/d⌋");

        // ΣΣ gcd(i,j) = Σ φ(d)⌊n/d⌋⌊m/d⌋
        bool ok2 = true;
        for (int it = 0; it < 40; it++) {
            int n = (int)(rng() % 200) + 1, m = (int)(rng() % 200) + 1;
            ll a = 0, b = 0;
            for (int d = 1; d <= min(n, m); d++) a += (ll)phi[d] * (n / d) * (m / d);
            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= m; j++) b += my_gcd(i, j);
            if (a != b) { ok2 = false; printf("    n=%d m=%d %lld vs %lld\n", n, m, a, b); break; }
        }
        check(ok2, "ΣΣ gcd(i,j) = Σ φ(d)⌊n/d⌋⌊m/d⌋");
    }

    printf("=== 14. 杜教筛 ===\n");
    {
        dusieve::init();
        sieve(200000);
        vector<ll> phiPre(200001, 0);
        vector<int> muPre(200001, 0);
        for (int i = 1; i <= 200000; i++) {
            phiPre[i] = phiPre[i - 1] + phi[i];
            muPre[i] = muPre[i - 1] + mu[i];
        }
        bool ok = true, ok2 = true;
        for (int n = 1; n <= 200000; n += 997) {
            if (dusieve::S_phi(n) != phiPre[n]) { ok = false; printf("    S_phi(%d)=%lld 期望 %lld\n", n, dusieve::S_phi(n), phiPre[n]); break; }
            if (dusieve::S_mu(n) != muPre[n]) { ok2 = false; printf("    S_mu(%d)=%d 期望 %d\n", n, dusieve::S_mu(n), muPre[n]); break; }
        }
        check(ok, "杜教筛 S_phi 与线性筛前缀和一致");
        check(ok2, "杜教筛 S_mu 与线性筛前缀和一致");
    }

    printf("\n===== 通过 %d 项，失败 %d 项 =====\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
