// ============================================================================
// bench_interval_tests2.cpp — 区间-区间运算法与特殊值（次正规/巨值/0）对照测试
// ============================================================================
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>

#include "nurbs/interval_arith.hpp"

using nurbs::Interval;

static bool ph(const std::string& s, double* out) {
    char* end = nullptr;
    *out = std::strtod(s.c_str(), &end);
    return end && *end == '\0';
}

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: %s file\n", argv[0]); return 1; }
    std::ifstream f(argv[1]);
    long n = 0;
    f >> n;
    long fails = 0, checked = 0;
    long by_op[4] = { 0, 0, 0, 0 };
    for (long i = 0; i < n; ++i) {
        std::string op, sa1, sa2, sb1, sb2, slo, shi;
        if (!(f >> op >> sa1 >> sa2 >> sb1 >> sb2 >> slo >> shi)) break;
        double a1, a2, b1, b2, lo, hi;
        ph(sa1, &a1); ph(sa2, &a2); ph(sb1, &b1); ph(sb2, &b2); ph(slo, &lo); ph(shi, &hi);
        Interval A(a1, a2), B(b1, b2), R;
        int idx = 0;
        if (op == "add") { R = nurbs::iadd(A, B); idx = 0; }
        else if (op == "sub") { R = nurbs::isub(A, B); idx = 1; }
        else if (op == "mul") { R = nurbs::imul(A, B); idx = 2; }
        else if (op == "div") { R = nurbs::idiv(A, B); idx = 3; }
        else { std::printf("bad op %s\n", op.c_str()); return 1; }
        ++checked;
        // 溢出时允许区间含 inf：只要求覆盖期望（lo<=lo, hi>=hi）；NaN 不接受
        bool ok = !std::isnan(R.lo) && !std::isnan(R.hi) && R.lo <= lo && R.hi >= hi;
        if (!ok) {
            ++fails; ++by_op[idx];
            if (fails <= 10)
                std::printf("FAIL %s A=[%a,%a] B=[%a,%a] expected [%a,%a] got [%a,%a]\n",
                            op.c_str(), a1, a2, b1, b2, lo, hi, R.lo, R.hi);
        }
    }
    std::printf("interval-interval tests: checked=%ld failures=%ld (add=%ld sub=%ld mul=%ld div=%ld)\n",
                checked, fails, by_op[0], by_op[1], by_op[2], by_op[3]);
    return fails == 0 ? 0 : 1;
}
