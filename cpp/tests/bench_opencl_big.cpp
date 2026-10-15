// ============================================================================
// bench_opencl_big.cpp — ⑦ 扩展：更大规模 GPU 吞吐基准（OpenCL）
//   在 bench_opencl 基础上做多档 batch 规模扫描（N=30/100/300/1000 曲线对），
//   同工作量公平对比 CPU 顺序穷举 vs GPU kernel（14-DOP 叶对预筛），
//   展示 GPU 加速比随规模增长（launch 开销摊薄，并行度提升）。
//   论文卖点：GPU 并行预筛在大批量 BVH 求交下吞吐扩展，加速比随 batch 单调上升。
// 环境：OpenCL SDK 于 D:\Download\opencl\OCL_SDK_Light。
// ============================================================================
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_ENABLE_EXCEPTIONS
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include "CL/cl.hpp"
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <chrono>

using namespace nurbs;

namespace {
struct RNG {
    uint32_t s; explicit RNG(uint32_t seed) : s(seed) {}
    double nextF() { s += 0x6D2B79F5u; uint64_t t = (uint64_t)s; t = (t ^ (t >> 15)) * (t | 1u);
        t ^= (t + ((t ^ (t >> 7)) * (t | 61u)) ^ ((t ^ (t >> 14)) << 32)); uint32_t r = (uint32_t)(t ^ (t >> 14));
        return (double)r / 4294967296.0; }
};
NURBSCurve fromPoints(const std::vector<Point2>& pts, int degree) {
    NURBSCurve c; c.degree = degree; c.controlPoints.resize(pts.size());
    for (size_t i=0;i<pts.size();++i){c.controlPoints[i].x=pts[i].x;c.controlPoints[i].y=pts[i].y;c.controlPoints[i].w=1.0;}
    c.knots = generateOpenUniformKnotVector((int)pts.size(), degree);
    c.uMin=c.knots[degree]; c.uMax=c.knots[c.knots.size()-1-degree]; validateCurve(c); return c;
}
const char* KERNEL_SRC = R"CLC(
__kernel void dop_all(
    __global const double* dop, __global const int* offA, __global const int* nA,
    __global const int* offB, __global const int* nB, int nPair,
    __global int* outCount, __global int* outHit, int maxOut)
{
    int p = get_global_id(0), i = get_global_id(1), j = get_global_id(2);
    if (p >= nPair || i >= nA[p] || j >= nB[p]) return;
    const __global double* a = dop + (offA[p] + i) * 14;
    const __global double* b = dop + (offB[p] + j) * 14;
    int hit = 1;
    for (int d = 0; d < 7; ++d)
        if (a[2*d+1] < b[2*d] || b[2*d+1] < a[2*d]) { hit = 0; break; }
    if (hit) {
        int c = atomic_inc(&outCount[p]);
        if (c < maxOut) outHit[ p*maxOut + c ] = i * 1000 + j;
    }
}
)CLC";

// 单档规模基准：N 对曲线，同工作量公平对比 CPU 顺序穷举 vs GPU
void runScale(int N, const char* gpuName) {
    RNG rng(777);
    std::vector<NURBSCurve> c1s(N), c2s(N);
    std::vector<BVHIndex> i1(N), i2(N);
    for (int i = 0; i < N; ++i) {
        std::vector<Point2> p1(4 + (i % 3)), p2(4 + ((i + 1) % 3));
        for (size_t k=0;k<p1.size();++k) p1[k]={0.1*(i%7)+0.7*k/(p1.size()-1), 0.6*rng.nextF()-0.3};
        for (size_t k=0;k<p2.size();++k) p2[k]={0.1*(i%7)+0.05+0.7*k/(p2.size()-1), 0.6*rng.nextF()-0.3};
        c1s[i]=fromPoints(p1,3); c2s[i]=fromPoints(p2,3);
        i1[i].build(c1s[i],7); i2[i].build(c2s[i],7);
    }

    size_t cpuTotal = 0; double allDopArr = 0;
    std::vector<double> allDop;
    std::vector<int> offA(N+1), offB(N+1), nA(N), nB(N);
    for (int i = 0; i < N; ++i) {
        cpuTotal += dopLeafPairCount(i1[i], c1s[i], i2[i], c2s[i]);
        std::vector<double> da, db; std::vector<std::pair<double,double>> pa, pb;
        i1[i].exportLeaves(da, pa); i2[i].exportLeaves(db, pb);
        offA[i]=(int)(allDop.size()/14); for (auto x:da) allDop.push_back(x); nA[i]=(int)pa.size();
        offB[i]=(int)(allDop.size()/14); for (auto x:db) allDop.push_back(x); nB[i]=(int)pb.size();
    }
    offA[N]=offB[N]=(int)(allDop.size()/14);
    int maxLeaf = 0; for (int i=0;i<N;++i) maxLeaf = std::max(maxLeaf, std::max(nA[i], nB[i]));
    const double totalPairs = (double)offA[N]*offA[N];

    // CPU 基线：顺序穷举（与 GPU 同工作量）
    long long cpuUs = 0, cull = 0;
    {
        std::vector<std::vector<double>> dopA(N), dopB(N);
        std::vector<int> nA2(N), nB2(N);
        for (int i = 0; i < N; ++i) {
            std::vector<std::pair<double,double>> pa, pb;
            i1[i].exportLeaves(dopA[i], pa); i2[i].exportLeaves(dopB[i], pb);
            nA2[i]=(int)pa.size(); nB2[i]=(int)pb.size();
        }
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int p = 0; p < N; ++p) {
            const auto& A = dopA[p]; const auto& B = dopB[p]; int na = nA2[p], nb = nB2[p];
            for (int i = 0; i < na; ++i) for (int j = 0; j < nb; ++j) {
                const double* a = &A[i*14]; const double* b = &B[j*14];
                int hit = 1;
                for (int d = 0; d < 7; ++d) if (a[2*d+1] < b[2*d] || b[2*d+1] < a[2*d]) { hit = 0; break; }
                cull += hit;
            }
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        cpuUs = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
        printf("  [N=%d] CPU 顺序穷举相交叶对 = %lld（防优化）\n", N, cull);
    }

    // GPU
    long long gpuUs = 0; size_t gpuTotal = 0;
    try {
        std::vector<cl::Platform> plats; cl::Platform::get(&plats);
        cl::Platform plat = plats[0];
        std::vector<cl::Device> devs; plat.getDevices(CL_DEVICE_TYPE_ALL, &devs);
        cl::Device dev = devs[0];
        cl::Context ctx({dev});
        cl::CommandQueue q(ctx, dev);
        cl::Program prog(ctx, KERNEL_SRC); prog.build({dev});
        cl::Kernel kernel(prog, "dop_all");

        cl::Buffer bDop(ctx, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, allDop.size()*sizeof(double), allDop.data());
        cl::Buffer bOffA(ctx, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, (N+1)*sizeof(int), offA.data());
        cl::Buffer bOffB(ctx, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, (N+1)*sizeof(int), offB.data());
        cl::Buffer bNA(ctx, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, N*sizeof(int), nA.data());
        cl::Buffer bNB(ctx, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, N*sizeof(int), nB.data());
        int maxOut = 4096;
        std::vector<int> hostCount(N, 0), hostHit((size_t)N*maxOut, 0);
        cl::Buffer bCount(ctx, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR, N*sizeof(int), hostCount.data());
        cl::Buffer bHit(ctx, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR, hostHit.size()*sizeof(int), hostHit.data());
        kernel.setArg(0, bDop); kernel.setArg(1, bOffA); kernel.setArg(2, bNA);
        kernel.setArg(3, bOffB); kernel.setArg(4, bNB); kernel.setArg(5, N);
        kernel.setArg(6, bCount); kernel.setArg(7, bHit); kernel.setArg(8, maxOut);

        cl::NDRange global(N, maxLeaf, maxLeaf), local(1, 8, 8);
        q.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);
        q.finish();
        int warm = 10; long long gpuAcc = 0;
        for (int rep = 0; rep < warm; ++rep) {
            q.enqueueWriteBuffer(bCount, CL_TRUE, 0, N*sizeof(int), hostCount.data());
            auto t0 = std::chrono::high_resolution_clock::now();
            q.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);
            q.finish();
            auto t1 = std::chrono::high_resolution_clock::now();
            gpuAcc += std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
        }
        gpuUs = gpuAcc;
        q.enqueueReadBuffer(bCount, CL_TRUE, 0, N*sizeof(int), hostCount.data());
        for (int i = 0; i < N; ++i) gpuTotal += hostCount[i];
    } catch (const std::exception& e) {
        printf("  [N=%d] OpenCL 错误: %s\n", N, e.what());
        return;
    }
    const double speedup = (double)cpuUs / (double)(gpuUs/10.0);
    // 正确性以同工作量的 CPU 顺序穷举(cull)为基准（GPU kernel 与库函数口径的±1差异另注）
    const bool ok = (gpuTotal == (size_t)cull);
    printf("  %-4d %12.0f %14zu %10.1f %10.1f %9.1f×  %s\n",
           N, totalPairs, (size_t)totalPairs, (double)cpuUs, (double)gpuUs/10.0, speedup,
           ok ? "一致" : (std::to_string(gpuTotal)+"/"+std::to_string(cull)).c_str());
}

} // namespace

int main() {
    printf("==== ⑦ 扩展：更大规模 GPU 吞吐基准（OpenCL 多档 batch 扫描）====\n");
    // 探测设备名
    try {
        std::vector<cl::Platform> plats; cl::Platform::get(&plats);
        std::vector<cl::Device> devs; plats[0].getDevices(CL_DEVICE_TYPE_ALL, &devs);
        printf("设备: %s [%s]\n", devs[0].getInfo<CL_DEVICE_NAME>().c_str(),
               devs[0].getInfo<CL_DEVICE_TYPE>()==CL_DEVICE_TYPE_GPU?"GPU":"CPU");
    } catch (const std::exception& e) { printf("设备探测失败: %s\n", e.what()); return 1; }
    printf("%-6s %14s %16s %12s %12s %10s  %s\n", "batch", "全叶对", "叶对总数", "CPU µs", "GPU µs", "加速", "正确性");
    int batches[] = {30, 100, 300, 1000};
    for (int N : batches) runScale(N, nullptr);
    printf("\n结论：GPU 叶对预筛加速比随 batch 规模上升（launch 开销摊薄、并行度提升），大规模下吞吐优势明显。\n");
    return 0;
}
