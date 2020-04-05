#!/bin/bash

ifort -132 -check all -fpe0 -f77rtl -o 222_sr srchwave591.JdF.f
ifort -132 -check all -fpe0 -f77rtl -o 111_rd rdsetupsimul.1.JdF.f -lm /usr/local/sac/lib/libsacio.a