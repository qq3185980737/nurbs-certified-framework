================================================================================
NURBS 统一数值框架几何计算系统 - 本地运行指南
================================================================================

【重要：必须以本地 HTTP 服务器方式打开】
数学核使用 WebAssembly（cpp/wasm/nurbs.wasm），.wasm 通过 fetch 加载；
若直接双击 index.html（file:// 协议），浏览器会拦截 .wasm 加载，
导致页面初始化中断：画布与控制点输入区空白、无法输入。
请务必用本地 HTTP 服务器打开，方法见下方【快速启动】。

【快速启动】

方式一：一键脚本（推荐）
  双击本目录下「启动服务器.bat」，自动启动本地服务器并打开浏览器。

方式二：命令行手动
  1) Win+R 输入 cmd 回车，打开命令行
  2) 输入：cd /d D:\qqdownload\nurbs-certified-framework
  3) 输入：python -m http.server 8000
     （保持命令行窗口开启；若提示 python 不是命令，试 python3 或 py）
  4) 浏览器访问：http://localhost:8000/index.html

方式三：VS Code Live Server
  用 VS Code 打开本目录，安装「Live Server」扩展，点状态栏「Go Live」。

【支持的浏览器】
✅ Google Chrome (推荐)
✅ Microsoft Edge
✅ Mozilla Firefox
（任选其一；务必通过 http://localhost 访问，而非 file://）

【常见问题解决】

问题 1: 双击 index.html 后页面空白 / 控制点区空 / 无法输入
原因: 数学核 .wasm 在 file:// 协议下被浏览器拦截，初始化中断。
解决: 改用【快速启动】中的任一方式，通过 http://localhost 访问。

问题 2: 页面样式丢失
解决: TailwindCSS 使用 CDN，请检查网络连接、刷新(F5)或清除缓存。

问题 3: 点击「导出数据」无反应
解决: 允许浏览器弹窗；或在地址栏右侧允许被拦截的弹窗。

问题 4: 控制台显示模块 / WASM 加载错误
解决: 确认 cpp/wasm/nurbs.js 与 nurbs.wasm 与本文件在同一目录；通过 http 访问。

问题 5: Canvas 不显示或交互无响应
解决: 按 F12 查看控制台错误；调整窗口大小触发重绘；使用现代浏览器。

【文件结构说明】
index.html                     - 主页面文件（入口，须经 HTTP 服务器访问）
cpp/wasm/nurbs.js + nurbs.wasm - C++ 编译的 WebAssembly 数学核
src/                           - JavaScript/Web 前端源码（经 build.mjs 生成 index.html 内联段）
cpp/                           - C++17 参考实现、C API、WASM/OpenCL、基准与测试
启动服务器.bat                 - 一键启动本地服务器并打开页面

【功能特性】
- 两条 NURBS 曲线的交互式编辑（控制点 x/y 直接输入或画布拖拽）
- 自动判别：相交 / 相切 / 分离状态
- 保证性求交（区间算术 + Krawczyk 认证）：0 假阳性 / 0 假证书
- 退化几何恢复、曲面自交、共面重叠检测
- 包围盒预筛选、BVH 加速、批量并行
- 计算结果导出为 JSON

【技术说明】
- C++17 数学核（MSVC Release /O2 /fp:precise）编译为 WASM，浏览器与 C++ 共享同一算法（单源）
- 前端接入层由 build.mjs 从 src/ 生成内联进 index.html（勿手改内联段）
- TailwindCSS (CDN) 样式框架
- 原生 Canvas 2D 渲染

================================================================================
如需进一步技术支持，请查看浏览器开发者工具 (F12) 的控制台输出。
================================================================================
