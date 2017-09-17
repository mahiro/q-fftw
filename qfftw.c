#include <fftw3.h>
#include <k.h>

#ifdef __cplusplus
extern "C"{
#endif

int _type(K x) {
    char t = x->t;
    return (t == 1) || (4 <= t && t <= 10) || (12 <= t && t <= 19);
}

void _copy_k2r(K x, double* c) {
    int i, n = x->n;

    switch (x->t) {
    case KB: /*  1: q boolean   (c char) */
    case KG: /*  4: q byte      (c char) */
        for (i = 0; i < n; i++) c[i] = (double) kG(x)[i];
        break;
    case KH: /*  5: q short     (c short) */
        for (i = 0; i < n; i++) c[i] = (double) kH(x)[i];
        break;
    case KI: /*  6: q int       (c int) */
    case KM: /* 13: q month     (c int) */
    case KD: /* 14: q date      (c int) */
    case KU: /* 17: q minute    (c int) */
    case KV: /* 18: q second    (c int) */
    case KT: /* 19: q time      (c int) */
        for (i = 0; i < n; i++) c[i] = (double) kI(x)[i];
        break;
    case KJ: /*  7: q long      (c int64_t) */
    case KP: /* 12: q timestamp (c int64_t) */
    case KN: /* 16: q timespan  (c int64_t) */
        for (i = 0; i < n; i++) c[i] = (double) kJ(x)[i];
        break;
    case KE: /*  8: q real      (c float) */
        for (i = 0; i < n; i++) c[i] = (double) kE(x)[i];
        break;
    case KF: /*  9: q float     (c double) */
    case KZ: /* 15: q datetime  (c double) */
        for (i = 0; i < n; i++) c[i] = (double) kF(x)[i];
        break;
    case KC: /* 10: q char      (c char) */
        for (i = 0; i < n; i++) c[i] = (double) kC(x)[i];
        break;
    }
}

void _copy_k2c(K x, fftw_complex* c, int imag) {
    int i, n = x->n;

    switch (x->t) {
    case KB: /*  1: q boolean   (c char) */
    case KG: /*  4: q byte      (c char) */
        for (i = 0; i < n; i++) c[i][imag] = (double) kG(x)[i];
        break;
    case KH: /*  5: q short     (c short) */
        for (i = 0; i < n; i++) c[i][imag] = (double) kH(x)[i];
        break;
    case KI: /*  6: q int       (c int) */
    case KM: /* 13: q month     (c int) */
    case KD: /* 14: q date      (c int) */
    case KU: /* 17: q minute    (c int) */
    case KV: /* 18: q second    (c int) */
    case KT: /* 19: q time      (c int) */
        for (i = 0; i < n; i++) c[i][imag] = (double) kI(x)[i];
        break;
    case KJ: /*  7: q long      (c int64_t) */
    case KP: /* 12: q timestamp (c int64_t) */
    case KN: /* 16: q timespan  (c int64_t) */
        for (i = 0; i < n; i++) c[i][imag] = (double) kJ(x)[i];
        break;
    case KE: /*  8: q real      (c float) */
        for (i = 0; i < n; i++) c[i][imag] = (double) kE(x)[i];
        break;
    case KF: /*  9: q float     (c double) */
    case KZ: /* 15: q datetime  (c double) */
        for (i = 0; i < n; i++) c[i][imag] = (double) kF(x)[i];
        break;
    case KC: /* 10: q char      (c char) */
        for (i = 0; i < n; i++) c[i][imag] = (double) kC(x)[i];
        break;
    }
}

void _copy_r2k(K x, double* c) {
    int i, n = x->n;

    for (i = 0; i < n; i++) {
        kF(x)[i] = c[i];
    }
}

void _copy_c2k(K x, fftw_complex* c, int imag) {
    int i, n = x->n;

    for (i = 0; i < n; i++) {
        kF(x)[i] = c[i][imag];
    }
}

K _dft_r2r(K x, fftw_r2r_kind kind) {
    double *in, *out;
    int n = x->n;
    fftw_plan p;
    K y;

    if (!_type(x)) {
        return krr("type");
    }

    in  = fftw_alloc_real(n);
    out = fftw_alloc_real(n);

    _copy_k2r(x, in);

    p = fftw_plan_r2r(1, &n, in, out, &kind, FFTW_ESTIMATE);
    fftw_execute(p);

    y = ktn(KF, n);

    _copy_r2k(y, out);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return y;
}

K _dft_c2c(K x, int sign) {
    int i, n;
    fftw_complex *in, *out;
    fftw_plan p;
    K x_re, x_im;
    K y_re, y_im, y;

    if (x->t != 0) {
        return krr("type");
    }

    if (x->n != 2) {
        return krr("rank");
    }

    for (i = 0; i < 2; i++) {
        if (!_type(kK(x)[i])) {
            return krr("type");
        }
    }

    x_re = kK(x)[0];
    x_im = kK(x)[1];

    if (x_re->n != x_im->n) {
        return krr("rank");
    }

    n = x_re->n;

    in  = fftw_alloc_complex(n);
    out = fftw_alloc_complex(n);

    _copy_k2c(x_re, in, 0);
    _copy_k2c(x_im, in, 1);

    p = fftw_plan_dft(1, &n, in, out, sign, FFTW_ESTIMATE);
    fftw_execute(p);

    y_re = ktn(KF, n);
    y_im = ktn(KF, n);
    y = knk(2, y_re, y_im);

    _copy_c2k(y_re, out, 0);
    _copy_c2k(y_im, out, 1);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    
    return y;
}

K  dft(K x) {return _dft_c2c(x, FFTW_FORWARD);}
K idft(K x) {return _dft_c2c(x, FFTW_BACKWARD);}

K  dct(K x) {return _dft_r2r(x, FFTW_REDFT10);}
K idct(K x) {return _dft_r2r(x, FFTW_REDFT01);}

K  dst(K x) {return _dft_r2r(x, FFTW_RODFT10);}
K idst(K x) {return _dft_r2r(x, FFTW_RODFT01);}

K  dht(K x) {return _dft_r2r(x, FFTW_DHT);}

K r2hc(K x) {return _dft_r2r(x, FFTW_R2HC);}
K hc2r(K x) {return _dft_r2r(x, FFTW_HC2R);}

K dct1(K x) {return _dft_r2r(x, FFTW_REDFT00);}
K dct2(K x) {return _dft_r2r(x, FFTW_REDFT10);}
K dct3(K x) {return _dft_r2r(x, FFTW_REDFT01);}
K dct4(K x) {return _dft_r2r(x, FFTW_REDFT11);}

K dst1(K x) {return _dft_r2r(x, FFTW_RODFT00);}
K dst2(K x) {return _dft_r2r(x, FFTW_RODFT10);}
K dst3(K x) {return _dft_r2r(x, FFTW_RODFT01);}
K dst4(K x) {return _dft_r2r(x, FFTW_RODFT11);}

#ifdef __cplusplus
}
#endif
