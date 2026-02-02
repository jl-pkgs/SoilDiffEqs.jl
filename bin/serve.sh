#!/bin/bash

# 后台运行 opencode serve
# 终端关闭后不影响进程

nohup opencode serve > opencode.log 2>&1 &

echo "后台进程已启动: $!"
echo "日志: opencode.log"
