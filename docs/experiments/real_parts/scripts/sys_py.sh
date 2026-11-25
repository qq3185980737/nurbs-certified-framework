#!/bin/bash
# sys_py.sh — 用系统 python3 运行脚本（纯文本处理，不需要 OCC）
set -u
exec python3 "$@"
