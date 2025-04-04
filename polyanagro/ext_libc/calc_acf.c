#ifndef __RG_H
#define __RG_H

#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
  #include <omp.h>
  #define USED_OPENMP 1
#else
  #define USED_OPENMP 0
#endif

static void c_acf_e2e(int ndim, int nchains, int ndumps, float uree[ndim][nchains][ndumps], float rEE_ACF[])
{

    /*
    !This subroutine calculates the end to end autocorrelation
    !function
    !
    !  EEACF = <u(t+to)*u(to)> where u(t) = R(t)/|R(t)| (Normalize end2end vector)
    */

    int i, j, i1, i2, k, n_it, maxiter;
    float dot_result;

    for (i=0; i<ndumps; i++) rEE_ACF[i] = 0.0;

    // Dump interactions
    maxiter = ndumps - 1;
    rEE_ACF[0] = 1.0;
    for (i=1; i<maxiter; i++) {
        n_it = 0;
        // Value of the t0 time
        for (j=0; j<(maxiter-i); j++) {
            i1 = j;
            i2 = j+i;
            for (k=0; k<nchains; k++) {
                dot_result = 0.0;
                for (int idim = 0; idim < ndim; idim++) {
                    dot_result += uree[idim][k][i1]*uree[idim][k][i2];
                }
                rEE_ACF[i] = rEE_ACF[i] + dot_result;
                n_it += 1;
            }
        }
        if (n_it != 0) rEE_ACF[i] = rEE_ACF[i] / (float) n_it;
    }

}

static void c_acf2_e2e(int ndim, int nchains, int ndumps, float ree[ndim][nchains][ndumps], float rEE_ACF[])
{

    /*
    !This subroutine calculates the end to end autocorrelation
    !function
    !
    !  EEACF = <u(t+to)*u(to)> where u(t) = R(t)/|R(t)| (Normalize end2end vector)
    */

    int i, j, i1, i2, k, n_it, maxiter;
    float dot_result;

    for (i=0; i<ndumps; i++) rEE_ACF[i] = 0.0;

    // Dump interactions
    maxiter = ndumps - 1;
    //rEE_ACF[0] = 10.0;
    for (i=0; i<maxiter; i++) {
        n_it = 0;
        // Value of the t0 time
        for (j=0; j<(maxiter-i); j++) {
            i1 = j;
            i2 = j+i;
            for (k=0; k<nchains; k++) {
                dot_result = 0.0;
                for (int idim = 0; idim < ndim; idim++) {
                    dot_result += ree[idim][k][i1]*ree[idim][k][i2];
                }
                rEE_ACF[i] = rEE_ACF[i] + dot_result;
                n_it += 1;
            }
        }
        if (n_it != 0) rEE_ACF[i] = rEE_ACF[i] / (float) n_it;
    }

    // Normalize
    for (i=0; i<maxiter; i++) {
        rEE_ACF[i] = rEE_ACF[i]/rEE_ACF[0];
    }

}

#endif