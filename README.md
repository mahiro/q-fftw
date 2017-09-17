# q-FFTW

`q-FFTW` enables FFTW 3.x (http://www.fftw.org/) for KDB+/q 3.x.

FFTW is a C-compatible library to compute Fourier transforms efficiently in `O(n log n)` time using FFT (Fast Fourier Transform) algorithms.

`q-FFTW` consists of two main files:

- dynamic library `qfftw.so`/`qfftw.dll` (which is statically linked to FFTW `libfftw3.a`) and
- utility q script `fftw.q` that loads functions into the `.fftw` directory namespace.

It can be installed in either Linux, MacOSX or Windows.

## Usage

See the [installation](#installation) section to install q-FFTW into your KDB+/q system.

Here is an example usage:

```
q)\l fftw.q

q)re: 1 0 1 0 1f
q)im: 0 1 0 1 0f

q).fftw.dft (re;im)
3 0.8632713  2.038842 -1.038842 0.1367287 
2 -0.1367287 1.038842 -2.038842 -0.8632713

q).fftw.dct re
6 1.110223e-16 1.236068 1.110223e-16 3.236068
```

- `.fftw.dft` computes DFT (discrete Fourier transform), complext-to-complex.
- `.fftw.dct` computes DCT (discrete cosine transform), real-to-real.

A complex number vector is represented as a mixed-list pair of the real-part vector and the imaginary-part vector, both of which are float lists of the same length.
An input vector can be of any numerical type as long as it can be converted to a float list.

The current version supports only 1-dimentional input/output vectors (rank 1).
Future versions can probably support arbitrary dimentions (ranks) as FFTW already supports.

## Installation

This section describes the installation procedure, which is unfortunately not (yet) as simple as `./configure && make && make install`.

For Windows, MinGW (http://www.mingw.org/) should be used to build library code.
On the other hand, the final `q-fftw` binary can be used in the native KDB+/q environment without MinGW.

### FFTW

You need to build FFTW first.

Go to http://www.fftw.org/download.html and download the latest source.

From the source directory, run these commands:

```
./configure
make
```

`[sudo] make install` is optional.

If your `q` architecture is 32bit (e.g. `$QHOME/l32/q`) while the host architecture is 64bit (which is most likely to be the case if you are using the free version of KDB+/q), run `./configure CFLAGS=-m32 LDFLAGS=-m32` instead, in order to build 32bit version of fftw.

```
./configure CFLAGS=-m32 LDFLAGS=-m32
make
```

Note: You may also need `gcc-multilib` to use `-m32` flags. (See https://stackoverflow.com/a/17748092)

### KDB

Clone the `KDB` source respository:

```
git clone https://github.com/KxSystems/kdb.git
```

### q-FFTW

Clone the `q-FFTW` source repository:

```
git clone https://github.com/mahiro/q-fftw.git
```

In the q-fftw directory, run `./build.sh` with three arguments:

- either "32" or "64" indicating the architecture of your KDB+/q installation
- KDB source directory path
- FFTW source directory path

The command should look something like this:

```
./build.sh 32 /path/to/source/kdb /path/to/source/fftw-3.x
```

For reference, this URL lists environment-specific compile commands:
http://code.kx.com/q/interfaces/using-c-functions/#compiling-extension-modules

If the build is successful, `qfftw.so` (Linux/MacOSX) or `qfftw.dll` (Windows) is generated in the same directory.
Copy or symlink this file so it is available right next to the `q` executable file (e.g. `$QHOME/l32/q`).

Also, copy or symlink `fftw.q` into `$QHOME` directory (right next to `q.q` file etc.). Alternatively, this file can be placed anywhere you can load via the `\l` system command in q.

### Verification

Once `qfftw.so` and `fftw.q` are set up, try the following commands to verify the installation.

```
q)\l fftw.q
q).fftw.dft
code
```

## API Functions

### Complex-to-Complex

DFT (discrete Fourier transform)

- `.fftw.dft`: DFT
- `.fftw.idft`: Inverse DFT

### Real-to-Real

DCT (discrete cosine transform)

- `.fftw.dct`: DCT -- REDFT10, a.k.a. DCT-II
- `.fftw.idct`: Inverse DCT -- REDFT01, a.k.a. DCT-III

DST (discrete sine transform)

- `.fftw.dst`: DST -- RODFT10, a.k.a. DST-II
- `.fftw.idst`: Inverse DCT -- RODFT01, a.k.a. DST-III

DHT (discrete Hartley transform)

- `.fftw.dht`: DHT -- same function for inverse

DFT Half Complex

- `.fftw.r2hc`: DFT (real to half-complex)
- `.fftw.hc2r`: Inverse DFT (half-complex to real)

DCT I-IV (Real Even DFTs)

- `.fftw.dct1`: REDFT00, a.k.a. DCT-I -- same function for inverse
- `.fftw.dct2`: REDFT10, a.k.a. DCT-II -- inverse is .fftw.dct3
- `.fftw.dct3`: REDFT01, a.k.a. DCT-III -- inverse is .fftw.dct2
- `.fftw.dct4`: REDFT11, a.k.a. DCT-IV -- same function for inverse

DST I-IV (Real Odd DFTs)

- `.fftw.dst1`: RODFT00, a.k.a. DST-I -- same function for inverse
- `.fftw.dst2`: RODFT10, a.k.a. DST-II -- inverse is .fftw.dst3
- `.fftw.dst3`: RODFT01, a.k.a. DST-III -- inverse is .fftw.dst2
- `.fftw.dst4`: RODFT11, a.k.a. DST-IV -- same function for inverse

