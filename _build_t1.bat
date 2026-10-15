@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
cl /nologo /O2 /fp:precise /std:c++17 /EHsc "D:\qqdownload\NURBS 统一数值框架几何计算系统\cpp\tests\bench_true_selfintersect.cpp" /Fe:"D:\qqdownload\NURBS 统一数值框架几何计算系统\cpp\build\bench_true_selfintersect.exe"
