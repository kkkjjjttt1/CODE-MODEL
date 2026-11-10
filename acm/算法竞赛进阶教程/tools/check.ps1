# 提取 Markdown 中所有 ```cpp 代码块并逐个做语法/语义检查
# 用法: powershell -NoProfile -ExecutionPolicy Bypass -File check.ps1 [-Root <教程目录>]
param(
    [string]$Root = "C:\Users\33317\Desktop\算法竞赛进阶教程"
)

$ErrorActionPreference = "Continue"
$utf8 = New-Object System.Text.UTF8Encoding($false)
$tmp = Join-Path $env:TEMP "acm_math_check"
if (Test-Path $tmp) { Remove-Item -Recurse -Force $tmp }
New-Item -ItemType Directory -Force -Path $tmp | Out-Null

# 公共前缀：只放「声明」，避免与代码块内的同名定义冲突，同时消除未声明标识符造成的误报。
# 注意：本地 MinGW GCC 8.1 的 <bits/stdc++.h> 自身有 bug，因此这里用具体头文件代替。
$header = @'
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
#include <unordered_set>
#include <queue>
#include <stack>
#include <bitset>
#include <random>
#include <chrono>
#include <functional>
#include <tuple>
#include <iostream>
using namespace std;
using ll = long long;
using ull = unsigned long long;
using i128 = __int128;
ll gcd(ll a, ll b);
ll lcm(ll a, ll b);
ll exgcd(ll a, ll b, ll &x, ll &y);
ll qpow(ll a, ll e, ll mod);
ll mul(ll a, ll b, ll mod);
bool isPrime(ll n);
ll pollard(ll n);
void factor(ll n, vector<ll> &res);
void sieve(int n);
ll euler_phi(ll n);
ll inv_exgcd(ll a, ll mod);
ll inv_fermat(ll a, ll mod);
ll primitive_root(ll p);
ll bsgs(ll a, ll b, ll p);
ll exbsgs(ll a, ll b, ll p);
ll C(int n, int k);
ll C_small(ll n, ll k, ll p);
ll lucas(ll n, ll k, ll p);
'@

$total = 0; $failed = 0; $counter = 0
$report = New-Object System.Collections.ArrayList

foreach ($md in (Get-ChildItem -Path $Root -Filter "*.md" | Sort-Object Name)) {
    $text = [System.IO.File]::ReadAllText($md.FullName, [System.Text.Encoding]::UTF8)
    $blocks = [regex]::Matches($text, '(?s)```cpp[ \t]*\r?\n(.*?)```')
    Write-Host ("=== {0}: {1} block(s)" -f $md.Name, $blocks.Count)
    $idx = 0
    $accum = ''''
    $accum = ""
    foreach ($b in $blocks) {
        $idx++; $total++; $counter++
        $code = $b.Groups[1].Value
        # 去掉 bits/stdc++.h（本地环境编不过），公共头已提供等价内容
        $code = [regex]::Replace($code, '(?m)^[ \t]*#include[ \t]*<bits/stdc\+\+\.h>[ \t]*\r?\n?', '')

        $ok = $false
        $out = @()
        $tag = ""
        $accumBody = $code
        if ($accum.Length -lt 40000) { $accumBody = $accum + "`n" + $code }
        $variants = @(
            @{ tag = "raw";     body = $code },
            @{ tag = "accum";   body = $accumBody },
            
            @{ tag = "wrapped"; body = "void __wrap__() {`n$code`n}" }
        )
        foreach ($v in $variants) {
            $tag = $v.tag
            $file = Join-Path $tmp ("check_{0}.cpp" -f $counter)
            [System.IO.File]::WriteAllText($file, ($header + "`n" + $v.body), $utf8)
            $out = & g++ -std=c++17 -fsyntax-only -w $file 2>&1
            if ($LASTEXITCODE -eq 0) { $ok = $true; break }
        }
        $accum += "`n" + $code
        if (-not $ok) {
            $failed++
            $msg = "FAIL {0} block #{1}  ({2})" -f $md.Name, $idx, $tag
            Write-Host $msg -ForegroundColor Red
            $out | Select-Object -First 8 | Write-Host
            [void]$report.Add(($msg + "`n" + ($out -join "`n")))
        }
    }
}

Write-Host ""
Write-Host ("总计 {0} 个代码块，失败 {1} 个。" -f $total, $failed)
if ($failed -gt 0) {
    $rp = Join-Path $env:TEMP "acm_math_check_report.txt"
    [System.IO.File]::WriteAllText($rp, ($report -join "`n`n"), $utf8)
    Write-Host ("详细报告: {0}" -f $rp)
}
