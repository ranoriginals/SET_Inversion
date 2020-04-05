#!/bin/bash
# 相速度反演需要编译两个程序，一个是srchwave591.SET.f,一个是rdsetupsimul.SET.f 
# -o后跟的是可执行文件的名字，这里我就分别命名为srchwave591 & rdsetupsimul 
# -lm后跟libsacio.a所在的路径
# 编译时会出现如下报错：
#srchwave591.SET.f(2157): warning #7960: The floating overflow condition was detected while evaluating this operation; the result is an Infinity.   [1.0E+50]
#      bmisfit = 1.0E+50
#----------------^
# 需要将srchwave591.SET.f中2157行中1.0E+50改为1.0E+38
ifort -132 -check all -fpe0 -f77rtl -o srchwave591 srchwave591.SET.f
ifort -132 -check all -fpe0 -f77rtl -o rdsetupsimul rdsetupsimul.SET.f -lm /usr/local/sac/lib/libsacio.a
