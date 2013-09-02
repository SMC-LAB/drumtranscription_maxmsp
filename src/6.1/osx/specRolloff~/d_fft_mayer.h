/*
 ** FFT and FHT routines
 **  Copyright 1988, 1993; Ron Mayer
 **  
 **  mayer_fht(fz,n);
 **      Does a hartley transform of "n" points in the array "fz".
 **  mayer_fft(n,real,imag)
 **      Does a fourier transform of "n" points of the "real" and
 **      "imag" arrays.
 **  mayer_ifft(n,real,imag)
 **      Does an inverse fourier transform of "n" points of the "real"
 **      and "imag" arrays.
 **  mayer_realfft(n,real)
 **      Does a real-valued fourier transform of "n" points of the
 **      "real" array.  The real part of the transform ends
 **      up in the first half of the array and the imaginary part of the
 **      transform ends up in the second half of the array.
 **  mayer_realifft(n,real)
 **      The inverse of the realfft() routine above.
 **      
 **      
 ** NOTE: This routine uses at least 2 patented algorithms, and may be
 **       under the restrictions of a bunch of different organizations.
 **       Although I wrote it completely myself, it is kind of a derivative
 **       of a routine I once authored and released under the GPL, so it
 **       may fall under the free software foundation's restrictions;
 **       it was worked on as a Stanford Univ project, so they claim
 **       some rights to it; it was further optimized at work here, so
 **       I think this company claims parts of it.  The patents are
 **       held by R. Bracewell (the FHT algorithm) and O. Buneman (the
 **       trig generator), both at Stanford Univ.
 **       If it were up to me, I'd say go do whatever you want with it;
 **       but it would be polite to give credit to the following people
 **       if you use this anywhere:
 **           Euler     - probable inventor of the fourier transform.
 **           Gauss     - probable inventor of the FFT.
 **           Hartley   - probable inventor of the hartley transform.
 **           Buneman   - for a really cool trig generator
 **           Mayer(me) - for authoring this particular version and
 **                       including all the optimizations in one package.
 **       Thanks,
 **       Ron Mayer; mayer@acuson.com
 **
 */

/* This is a slightly modified version of Mayer's contribution; write
 * msp@ucsd.edu for the original code.  Kudos to Mayer for a fine piece
 * of work.  -msp
 */
/* These pragmas are only used for MSVC, not MinGW or Cygwin <hans@at.or.at> */
#ifdef _MSC_VER
#pragma warning( disable : 4305 )  /* uncast const double to float */
#pragma warning( disable : 4244 )  /* uncast double to float */
#pragma warning( disable : 4101 )  /* unused local variables */
#endif

/* the following is needed only to declare pd_fft() as exportable in MSW */
//#include "m_pd.h"
#include <stdlib.h>
#include <string.h>

#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects	

#define REAL t_sample
#define GOOD_TRIG

#ifdef GOOD_TRIG
#else
#define FAST_TRIG
#endif

#if defined(GOOD_TRIG)
#define FHT_SWAP(a,b,t) {(t)=(a);(a)=(b);(b)=(t);}
#define TRIG_VARS                                                \
int t_lam=0;
#define TRIG_INIT(k,c,s)                                         \
{                                                           \
int i;                                                     \
for (i=2 ; i<=k ; i++)                                     \
{coswrk[i]=costab[i];sinwrk[i]=sintab[i];}             \
t_lam = 0;                                                 \
c = 1;                                                     \
s = 0;                                                     \
}
#define TRIG_NEXT(k,c,s)                                         \
{                                                           \
int i,j;                                                \
(t_lam)++;                                              \
for (i=0 ; !((1<<i)&t_lam) ; i++);                      \
i = k-i;                                                \
s = sinwrk[i];                                          \
c = coswrk[i];                                          \
if (i>1)                                                \
{                                                    \
for (j=k-i+2 ; (1<<j)&t_lam ; j++);                 \
j         = k - j;                                  \
sinwrk[i] = halsec[i] * (sinwrk[i-1] + sinwrk[j]);  \
coswrk[i] = halsec[i] * (coswrk[i-1] + coswrk[j]);  \
}                                                    \
}
#define TRIG_RESET(k,c,s)
#endif

#if defined(FAST_TRIG)
#define TRIG_VARS                                        \
REAL t_c,t_s;
#define TRIG_INIT(k,c,s)                                 \
{                                                    \
t_c  = costab[k];                                   \
t_s  = sintab[k];                                   \
c    = 1;                                           \
s    = 0;                                           \
}
#define TRIG_NEXT(k,c,s)                                 \
{                                                    \
REAL t = c;                                         \
c   = t*t_c - s*t_s;                                \
s   = t*t_s + s*t_c;                                \
}
#define TRIG_RESET(k,c,s)
#endif


/* --------------- fft -------------------- */
void mayer_fht(REAL *fz, int n);
void mayer_fft(int n, REAL *real, REAL *imag);
void mayer_ifft(int n, REAL *real, REAL *imag);
void mayer_realfft(int n, REAL *real);
void mayer_realifft(int n, REAL *real);