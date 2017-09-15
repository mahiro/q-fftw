#!/bin/sh
if [[ `uname` = Darwin ]]; then
    gcc -DKXVER=3 -bundle -undefined dynamic_lookup "$@" -lfftw3 -lm -o qfftw.so qfftw.c
else
    gcc -shared -fPIC "$@" -lfftw3 -lm -o qfftw.so qfftw.c
fi
