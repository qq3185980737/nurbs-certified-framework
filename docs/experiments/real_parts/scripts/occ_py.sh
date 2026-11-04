#!/bin/bash
# occ_py.sh — 用 OCCT 8.0.1 pythonocc 环境运行脚本
set -u
exec /root/miniconda3/envs/occ/bin/python "$@"
