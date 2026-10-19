@echo off
cd /d "%~dp0"
echo ============================================================
echo  NURBS Unified Numerical Framework - Local Server Launcher
echo  URL: http://localhost:8000/index.html
echo  Close this window to stop the server.
echo ============================================================
echo.
start "" /b python -m http.server 8000
timeout /t 2 >nul
start "" "http://localhost:8000/index.html"
echo.
echo Server is running. A browser tab should have opened.
echo Keep this window open while using the app.
pause
