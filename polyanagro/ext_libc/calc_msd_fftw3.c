#ifndef __MSD_FFTW3_H
#define __MSD_FFTW3_H

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <time.h>

#ifdef _OPENMP
  #include <omp.h>
  #define USED_OPENMP 1
#else
  #define USED_OPENMP 0
#endif

typedef double real;  // Replace with your actual type if needed (e.g., double)
// ======================================================================
void substractMean(real *signal, int nSignals, int signalSize) {

    real mean;
    int i, j;
    unsigned long idx;

    for (i=0; i<nSignals; i++) {
        // Calculate the mean of the dignal block
        // s(signal[i * signalSize] to signal[(i + 1) * signalSize - 1]

        mean = 0.0;
        for (j = 0; j < signalSize; j++) {
            idx = (unsigned long) i * signalSize + j;
            mean += signal[idx];
        }
        mean /= signalSize;

        // Subtract the mean from each element of the signal
        for (j = 0; j < signalSize; j++) {
            idx = (unsigned long) i * signalSize + j;
            signal[idx] -= mean;
        }
    }

    // DEBUG
    // Print the result for verification
//    for (i=0; i<signalSize*nSignals; i++) {
//        printf("%d: %f\n", i, signal[i]);
//    }
//    printf("////////\n");

}

// ======================================================================
real* autocorrCPU(real *signalOriginal, int signalSize, int nSignals) {

    // Pad with zeros
    int signalSizePad = signalSize * 2;
    unsigned long N = (unsigned long) 2 * (signalSizePad / 2 + 1);
    unsigned long totalSize =  (unsigned long) nSignals * N + 2;
    size_t n = (size_t) totalSize*sizeof(real);
    clock_t start, end;


    real *signal = (real *)fftw_malloc(n);
    if (!signal) {
        fprintf(stderr, "Memory allocation failed. Requested %zu bytes.\n", n);
        exit(EXIT_FAILURE);
    }
    memset(signal, 0, sizeof(real) * totalSize);  // if zero-initialization is needed

    start = clock();
    // Copy and zero-pad signal
    for (int i = 0; i < nSignals; i++) {
        for (int j = 0; j < signalSize; j++) {
            signal[j + i * N] = signalOriginal[j + i * signalSize];
        }
    }
    end = clock();
    printf("\t\t\t[Timing partial] AutoCorr (copy and zero pad) : %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    // FFT plan for forward transform
    start = clock();
    fftw_plan plan_r2c = fftw_plan_many_dft_r2c(
        1, &signalSizePad, nSignals,
        signal, NULL, 1, N,
        (fftw_complex *)signal, NULL, 1, signalSizePad / 2 + 1,
        FFTW_ESTIMATE
    );
    fftw_execute(plan_r2c);
    fftw_destroy_plan(plan_r2c);
    end = clock();
    printf("\t\t\t[Timing partial] AutoCorr (FFT forward) : %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    // Power spectrum
    start = clock();
    fftw_complex *signal_complex = (fftw_complex *)signal;
    for (int i = 0; i < (nSignals * (signalSizePad / 2 + 1)); i++) {
        real re = signal_complex[i][0];
        real im = signal_complex[i][1];
        signal_complex[i][0] = (re * re + im * im) / signalSizePad;
        signal_complex[i][1] = 0.0;
    }
    end = clock();
    printf("\t\t\t[Timing partial] AutoCorr (FFT Power spectrum) : %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    // FFT plan for inverse transform
    fftw_plan plan_c2r = fftw_plan_many_dft_c2r(
        1, &signalSizePad, nSignals,
        signal_complex, NULL, 1, signalSizePad / 2 + 1,
        signal, NULL, 1, N,
        FFTW_ESTIMATE
    );
    fftw_execute(plan_c2r);
    fftw_destroy_plan(plan_c2r);
    end = clock();
    printf("\t\t\t[Timing partial] AutoCorr (FFT inverse) : %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    // Normalize and store result
    real *res = (real *)malloc(sizeof(real) * signalSize * nSignals);
    for (int i = 0; i < nSignals; i++) {
        for (int j = 0; j < signalSize; j++) {
            res[i * signalSize + j] = signal[i * N + j] / (signalSize - j);
            //printf("%d %d %d %d %f\n", i, j, i * signalSize + j, i * signalSize + j, res[i * signalSize + j]);
        }
    }
    end = clock();
    printf("\t\t\t[Timing partial] AutoCorr (Normalize) : %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    free(signal);

    return res;
}

// ======================================================================
real* computeS2(real *signal, int signalSize, int nSignals, int dimensions) {

    unsigned long totalSignals = (unsigned long) nSignals * dimensions;
    unsigned long j;

    // Step 1: Get autocorrelation
    real *S2 = autocorrCPU(signal, signalSize, totalSignals);

    // Step 2: Allocate output and temporary buffer
    real *s2Averaged = (real *)malloc(sizeof(real) * signalSize * dimensions);
    real *sumAux = (real *)malloc(sizeof(real) * totalSignals);

    // Step 3: Loop through all time lags
    for (int i = 0; i < signalSize; i++) {

        // Progress bar
        if (i % (signalSize / 100 + 1) == 0 || i == signalSize - 1) {
            printf("\r\t\t[computeS2] Progress: %3d%%", (i * 100) / signalSize);
            fflush(stdout);
        }

        // sumAux reshuffle and normalize
        for (j = 0; j < totalSignals; j++) {
            int signalIndex = j / dimensions;
            int dimIndex = j % dimensions;
            sumAux[signalIndex + nSignals * dimIndex] = S2[i + j * signalSize] / nSignals;
        }

        // Step 4: For each dimension, sum across nSignals
        for (int d = 0; d < dimensions; d++) {
            real sum = 0.0;
            for (int k = 0; k < nSignals; k++) {
                sum += sumAux[k + nSignals * d];
            }
            s2Averaged[i + signalSize * d] = sum;
            //printf("%d %d %d %f\n", i, d, i + signalSize * d, s2Averaged[i + signalSize * d]);
        }
    }

    printf("\n");

    free(S2);
    free(sumAux);

    return s2Averaged;
}

// ======================================================================
real* computeS1(real *signal, int signalSize, int nSignals, int dimensions) {

    unsigned long totalSignals = (unsigned long) nSignals * dimensions;
    unsigned long j;

    real *D = (real *)malloc(sizeof(real) * signalSize);
    real *S1 = (real *)calloc(signalSize * dimensions, sizeof(real));  // Initialized to 0

    for (j = 0; j < totalSignals; j++) {
        int dimIndex = j % dimensions;

        // Step 1: Copy signal[j] and compute square into D
        for (int i = 0; i < signalSize; i++) {
            real val = signal[j * signalSize + i];
            D[i] = val * val;
        }

        // Step 2: Compute total Q
        real Q = 0.0;
        for (int i = 0; i < signalSize; i++) {
            Q += D[i];
        }
        Q *= 2.0;

        // Step 3: Slide through time lags
        for (int i = 0; i < signalSize; i++) {
            if (i > 0) {
                Q -= D[i - 1] + D[signalSize - i];
            }
            S1[i + signalSize * dimIndex] += Q / ((real)(signalSize - i) * nSignals);
            //printf("%d %d %d %f\n", j, i, i + signalSize * d, s2Averaged[i + signalSize * d]);

        }

        // Progress bar
        if (j % (totalSignals / 100 + 1) == 0 || j == totalSignals - 1) {
            printf("\r\t\t[computeS1] Progress: %3ld%%", (j * 100) / totalSignals);
            fflush(stdout);
        }
    }

    printf("\n");

//    // DEBUG
//    for (int i=0; i<signalSize * dimensions; i++) {
//        printf("%d %f\n", i, S1[i]);
//    }

    free(D);
    return S1;
}

// ======================================================================
void estimateMemoryUsage(int signalSize, int nSignals, int dimensions) {

     // Used in substractMean
    unsigned long elements_signal = (unsigned long) signalSize * nSignals * dimensions;
    size_t bytes_signal = (unsigned long) signalSize * nSignals * dimensions * sizeof(real);

    // Used in autocorrCPU
    unsigned long signalSizePad = (unsigned long) signalSize * 2;
    unsigned long N = (unsigned long) 2 * (signalSizePad / 2 + 1);
    unsigned long totalSize = nSignals * N + 2;
    size_t bytes_totalSize = (unsigned long) totalSize * sizeof(real);

    printf("[signalSize size]           %d elements\n", signalSize);
    printf("[nSignals size]             %d elements\n", nSignals);
    printf("[Dimensions size]           %d elements\n", dimensions);
    printf("[signal array elements]     %lu    elements\n", elements_signal);
    printf("[totalSize array elements]  %lu    elements\n", totalSize);
    printf("[signal array size]         ~%.2f GB\n", (double)bytes_signal/(1024*1024*1024));
    printf("[totalSize array size]      ~%.2f GB\n", (double)bytes_totalSize/(1024*1024*1024));
    printf("[Total memory size]         ~%.2f GB\n", (double)(bytes_signal+bytes_totalSize)/(1024*1024*1024));

}

// MAIN FUNCTION
double* c_msd_fftw3_fast(real *signal, int nSignals, int signalSize, int dimensions) {

    int totalSize = signalSize * dimensions;
    clock_t start, end;

    // DEBUG
    //printf("[Signal size]      %d elements\n", signalSize);
    //printf("[nSignals size]    %d elements\n", nSignals);
    //printf("[Dimensions size]  %d elements\n", dimensions);
    //estimateMemoryUsage(nSignals, signalSize, dimensions);
    // DEBUG

    printf("\t\t ========= TIME ========\n");

    // Step 1: Subtract the mean from the signal
    start = clock();
    printf("\t\t[Step] substractMean (Step 1 of 4)\n");
    substractMean(signal, nSignals*dimensions, signalSize);
    end = clock();
    printf("\t\t[Timing] substractMean (Step 1 of 4): %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    // Step 2: Compute S1 and S2
    start = clock();
    printf("\t\t[Step] ComputeS2 (Step 2 of 4)\n");
    real *S2 = computeS2(signal, signalSize, nSignals, dimensions);
    end = clock();
    printf("\t\t[Timing] computeS2 (Step 2 of 4): %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    real *S1 = computeS1(signal, signalSize, nSignals, dimensions);
    end = clock();
    printf("\t\t[Timing] computeS1 (Step 3 of 4): %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    // Step 4: Allocate and compute MSD = S1 - 2*S2
    start = clock();
    real *msd = (real *)malloc(sizeof(real) * totalSize);
    for (int i = 0; i < totalSize; i++) {
        msd[i] = S1[i] - 2.0 * S2[i];
        //printf("%d %f\n", i, msd[i]);
    }
    end = clock();
    printf("\t\t[Timing] compute MSD (Step 4 of 4): %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    free(S1);
    free(S2);

    return msd;
}


#endif