# q-FFTW

`q-FFTW` enables FFTW (http://www.fftw.org/) in KDB+/q.

FFTW is a C library to compute Fourier transforms efficiently in `O(n log n)` time using FFT (Fast Fourier Transform) algorithms.

`q-FFTW` provides two main components:

- dynamic library `qfftw.so` (which statically links FFTW `libfftw3.a`) and
- utility q script `fftw.q` to load functions into the `.fftw` directory namespace.

## Usage

See the [installation](#installation) section to install q-FFTW into your KDB+/q system.

Here is an example usage:

```
q)\l fftw.q
q)re: 0 1 0 -1 0f
q)im: 0 0 0 0 0f
q).fftw.dft (re;im)
0 1.118034  -1.118034 -1.118034  1.118034
0 -1.538842 0.3632713 -0.3632713 1.538842
q).fftw.dct re
0 2.351141 -4.440892e-16 -3.804226 1.110223e-16
```

- `.fftw.dft` computes DFT (discrete Fourier transform), complext-to-complex.
- `.fftw.dct` computes DCT (discrete cosine transform), real-to-real.

A complex number vector is represented as a mixed-list pair of the real-part vector and the imaginary-part vector, both of which are float lists of the same length.
An input vector can be of any numerical type as long as it can be converted to a float list.

The current version supports only 1-dimentional input/output vectors (rank 1).
Future versions can probably support arbitrary dimentions (ranks) as FFTW already supports.

## Installation

This section describes the installation procedure, which is unfortunately not (yet) as simple as `./configure && make && make install`.

### FFTW

You need to first install (or build) FFTW.
It can be either installed at the system level or simply compiled in the source directory.

Go to http://www.fftw.org/download.html and download the latest source (for Linux/MaxOS) or pre-compiled package (for Windows).

From the source directory, run these commands:

```
./configure
make
```

`[sudo] make install` is optional.

If your q architecture is 32-bit (e.g. `$QHOME/l32/q`) while the host architecture is 64-bit (which is most likely to be the case if you are using the free version of KDB+/q), run `./configure CFLAGS=-m32 LDFLAGS=-m32` instead, in order to build 32-bit version of fftw.

```
./configure CFLAGS=-m32 LDFLAGS=-m32
make
```

Note: You may also need `gcc-multilib` to use `-m32` flags. (See https://stackoverflow.com/a/17748092)

### KDB k.h

The header file `k.h` is required to build q-FFTW.

Retrieve the file from the KDB source (`./c/c/k.h`):

https://github.com/KxSystems/kdb

The header file should be saved in the q-FFTW directory.
Alternatively, specify the `-I` option for gcc to locate the directory where `k.h` is saved, as described in the next step.

### q-FFTW

Clone the `q-FFTW` repository:

```
git clone https://github.com/mahiro/q-fftw.git
```

Save the `k.h` file (either copy or symlink) in the q-FFTW directory.

See the below URL for environment-specific compile commands.

http://code.kx.com/q//interfaces/using-c-functions/#compiling-extension-modules

- Replace `bar.c` and `bar.so` with `qfftw.c` and `qfftw.so`.
- If FFTW was built from the source, reference the header and library path with `-I` and `-L`. These are not required if FFTW has been installed at the system level.
- Add `-lfftw3` and `-lm`.
- If your q architecture is 32-bit (e.g. `$QHOME/l32/q`) while the host architecture is 64-bit, add `-m32`.

Below is an example command for Linux:

```
FFTW_SRC=/path/to/source/fftw-3.x.x
OPTS="-I$FFTW_SRC/api -L$FFTW_SRC/.libs -lfftw3 -lm -m32"
gcc -shared -fPIC qfftw.c -o qfftw.so $OPTS
```

Once it is successful, `qfftw.so` file should be generated in the directory. Copy or symlink this file so it is available right next to the `q` executable file (e.g. `$QHOME/l32/q`).

Also, copy or symlink `fftw.q` into `$QHOME` directory (right next to `q.q` file etc.). Alternatively, this file can be placed anywhere you can load via the `\l` system command in q.

### Verification

Once `qfftw.so` and `fftw.q` are set up, try the following commands to verify the installation.

```
q)\l fftw.q
q).fftw.dft
code
```

## Functions

### Complex-to-Complex

DFT (discrete Fourier transform)

- .fftw.dft: DFT
- .fftw.idft: Inverse DFT

### Real-to-Real

DCT (discrete cosine transform)

- .fftw.dct: DCT -- REDFT10, a.k.a. DCT-II
- .fftw.idct: Inverse DCT -- REDFT01, a.k.a. DCT-III

DST (discrete sine transform)

- .fftw.dst: DST -- RODFT10, a.k.a. DST-II
- .fftw.idst: Inverse DCT -- RODFT01, a.k.a. DST-III

DHT (discrete Hartley transform)

- .fftw.dht: DHT -- same function for inverse

DFT Half Complex

- .fftw.r2hc: DFT (real to half-complex)
- .fftw.hc2r: Inverse DFT (half-complex to real)

DCT I-IV (Real Even DFTs)

- .fftw.dct1: REDFT00, a.k.a. DCT-I -- same function for inverse
- .fftw.dct2: REDFT10, a.k.a. DCT-II -- inverse is .fftw.dct3
- .fftw.dct3: REDFT01, a.k.a. DCT-III -- inverse is .fftw.dct2
- .fftw.dct4: REDFT11, a.k.a. DCT-IV -- same function for inverse

DST I-IV (Real Odd DFTs)

- .fftw.dst1: RODFT00, a.k.a. DST-I -- same function for inverse
- .fftw.dst2: RODFT10, a.k.a. DST-II -- inverse is .fftw.dst3
- .fftw.dst3: RODFT01, a.k.a. DST-III -- inverse is .fftw.dst2
- .fftw.dst4: RODFT11, a.k.a. DST-IV -- same function for inverse

