// ============================================================================
// bench_interval_tests.cpp — 区间算术核心的精确对照测试
//   读取 interval_tests.txt（Python Fraction 生成的期望向外区间），
//   校验 nurbs::iadd/isub/imul/idiv 的结果包含期望区间。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>

#include "nurbs/interval_arith.hpp"

using nurbs::Interval;

static bool parse_hex(const std::string& s, double* out) {
    char* end = nullptr;
    *out = std::strtod(s.c_str(), &end);
    return end && *end == '\0';
}

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: %s <interval_tests.txt>\n", argv[0]); return 1; }
    std::ifstream f(argv[1]);
    if (!f) { std::printf("cannot open %s\n", argv[1]); return 1; }
    long n = 0;
    f >> n;
    long fails = 0, checked = 0;
    long per_op_fail[4] = { 0, 0, 0, 0 };
    for (long i = 0; i < n; ++i) {
        std::string op;
        if (!(f >> op)) break;
        if (op == "div") {
            std::string sa, sb1, sb2, slo, shi;
            if (!(f >> sa >> sb1 >> sb2 >> slo >> shi)) break;
            double a, blo, bhi, lo, hi;
            parse_hex(sa, &a); parse_hex(sb1, &blo); parse_hex(sb2, &bhi);
            parse_hex(slo, &lo); parse_hex(shi, &hi);
            Interval r = nurbs::idiv(Interval(a, a), Interval(blo, bhi));
            ++checked;
            if (r.isEmpty() || !(r.lo <= lo && r.hi >= hi)) {
                ++fails; ++per_op_fail[3];
                if (fails <= 10) std::printf("FAIL div: a=%.17g div=[%.17g,%.17g] expected [%a,%a] got [%a,%a]\n",
                                             a, blo, bhi, lo, hi, r.lo, r.hi);
            }
        } else {
            std::string sa, sb, slo, shi;
            if (!(f >> sa >> sb >> slo >> shi)) break;
            double a, b, lo, hi;
            parse_hex(sa, &a); parse_hex(sb, &b); parse_hex(slo, &lo); parse_hex(shi, &hi);
            Interval r;
            int idx = 0;
            if (op == "add") { r = nurbs::iadd(Interval(a, a), Interval(b, b)); idx = 0; }
            else if (op == "sub") { r = nurbs::isub(Interval(a, a), Interval(b, b)); idx = 1; }
            else if (op == "mul") { r = nurbs::imul(Interval(a, a), Interval(b, b)); idx = 2; }
            else { std::printf("unknown op %s\n", op.c_str()); return 1; }
            ++checked;
            if (r.isEmpty() || !(r.lo <= lo && r.hi >= hi)) {
                ++fails; ++per_op_fail[idx];
                if (fails <= 10) std::printf("FAIL %s: a=%.17g b=%.17g expected [%a,%a] got [%a,%a]\n",
                                             op.c_str(), a, b, lo, hi, r.lo, r.hi);
            }
        }
    }
    std::printf("interval tests: checked=%ld, failures=%ld (add=%ld sub=%ld mul=%ld div=%ld)\n",
                checked, fails, per_op_fail[0], per_op_fail[1], per_op_fail[2], per_op_fail[3]);
    return fails == 0 ? 0 : 1;
}
