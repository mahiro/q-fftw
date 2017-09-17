#!/bin/bash
if [[ $# -lt 2 ]]; then
    echo "usage: $0 ARCH KDB_SRC FFTW_SRC"
    exit 1
fi

arch=$1; shift
kdb_src="$1"; shift
fftw_src="$1"; shift

mflag=
if [[ `uname -m` = *_64 && $arch = 32 ]]; then
    mflag=-m32
fi

if [[ `uname` = 'Darwin' ]]; then
    set -x
    gcc -bundle -undefined dynamic_lookup -DKXVER=3 $mflag -O3 \
        -I"$kdb_src/c/c" -I"$fftw_src/api" -L"$fftw_src/.libs" \
        qfftw.c -lfftw3 -lm -o qfftw.so "$@"
elif [[ `uname` = MINGW* ]]; then
    if [[ ! -e libq.def ]]; then
	    echo 'LIBRARY q.exe' > libq.def
	    echo 'EXPORTS' >> libq.def
	    nm -p "$kdb_src/w$arch/q.lib" | grep ' T _' | sed 's/.* T _//' >> libq.def
    fi
    if [[ ! -e libq.a ]]; then
	    dlltool -v -l libq.a -d libq.def
    fi
    set -x
    gcc -shared -static-libgcc -DKXVER=3 $mflag -O3 \
        -I"$kdb_src/c/c" -I"$fftw_src/api" -L"$fftw_src/.libs" -L"." \
        qfftw.c -lq -lfftw3 -lm -o qfftw.dll "$@"
else
    set -x
    gcc -shared -fPIC -DKXVER=3 -DKXVER=3 $mflag -O3 \
        -I"$kdb_src/c/c" -I"$fftw_src/api" -L"$fftw_src/.libs" \
        qfftw.c -lfftw3 -lm -o qfftw.so "$@"
fi
