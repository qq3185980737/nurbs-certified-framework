@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
cd /d "D:\qqdownload\nurbs-certified-framework\cpp"
cl /nologo /O2 /fp:precise /std:c++17 /utf-8 /EHsc /Iinclude tests\bench_surface_overlap.cpp /Fe:build\bench_surface_overlap.exe > build\_ov.log 2>&1
