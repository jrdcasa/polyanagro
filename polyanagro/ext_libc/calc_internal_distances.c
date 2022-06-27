#ifndef __INTERNALDIST_H
#define __INTERNALDIST_H

#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
  #include <omp.h>
  #define USED_OPENMP 1
#else
  #define USED_OPENMP 0
#endif

float *rsq_kremer;
int *rsqcount_kremer;
int maxnbondsperchain;

static void c_setup_internal_distances(int maxnbondsperch) {

    // Maximum number of backbone bonds per chain. TODO: Test if it works with different length chains
    maxnbondsperchain = maxnbondsperch; //nbin

    if ( ((rsq_kremer = (float*) malloc(sizeof(float) * (maxnbondsperchain+1)) ) == NULL) ||
         ((rsqcount_kremer = (int*) malloc(sizeof(int) * (maxnbondsperchain+1)) ) == NULL)) {
        fprintf(stderr, "c_internal_distances.pyx: out of memory\n");
        exit(0);
    }

    for (int i=0;i<=maxnbondsperchain;i++) {
        rsq_kremer[i] = 0.0;
        rsqcount_kremer[i] = 0;
    }


}

static void c_internal_distances_iframe(int nchains, int* head_array, int* ichatbb, float* coords) {

    // Calculate:
    // (0) -- (1) -- (2) -- (3) ....... (n)
    //
    //  rsq(j-i) = r(j)**2 - r(i)**2

    int iatom_o, iatom_l;
    int ilength;
    int ndim = 3;
    float rsq;
    float* r;
    float acc;
    FILE *javi = fopen("abb.txt","a");

    if ( ((r = (float*) malloc(sizeof(float) * (ndim)) ) == NULL) ) {
        fprintf(stderr, "c_internal_distances.pyx: out of memory\n");
        exit(0);
    }

    // Initialize
    for (int i=0; i<ndim; i++) r[i] = 0.0;

    // For each chain
    acc = 0.0;
    for (int ich=0; ich<nchains; ich++) {
        iatom_o = head_array[ich];
        iatom_l = ichatbb[iatom_o];
        // Pivot atom
        while (1) {
            ilength = 1;
            // next atoms to pivot
            while (iatom_l >= 0) {
                r[0] = coords[iatom_l*ndim+0] - coords[iatom_o*ndim+0];
                r[1] = coords[iatom_l*ndim+1] - coords[iatom_o*ndim+1];
                r[2] = coords[iatom_l*ndim+2] - coords[iatom_o*ndim+2];
                rsq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
                rsq_kremer[ilength] += rsq;
                rsqcount_kremer[ilength] += 1;
                if (ilength == 1) {
                   acc += rsq;
                   FILE *fptr = fopen("a.txt","a");
                   fprintf(fptr, "%d %f %f %d %d %f %f %d\n", ich, rsq, sqrt(rsq), iatom_l+1, iatom_o+1, acc, rsq_kremer[ilength],  rsqcount_kremer[ilength]);
                   fclose(fptr);

                   fprintf(javi, "%d %d %d \n", ich, iatom_l+1, iatom_o+1);
                   //fclose(fptr2);
                }
                iatom_l = ichatbb[iatom_l];
                ilength += 1;
            }
            if (ichatbb[iatom_o]<=0) break;
            iatom_o = ichatbb[iatom_o];
            iatom_l = ichatbb[iatom_o];
        }
    }
    fprintf(javi, "==========================\n");
    fclose(javi);
}

static int  c_internal_distances_return (float* rsq, int* rsqcount, float* rsqavg_intdist ) {

    for (int i=1; i<=maxnbondsperchain; i++) {
        rsq[i] = rsq_kremer[i];
        rsqcount[i] = rsqcount_kremer[i];
        rsqavg_intdist[i] = rsq_kremer[i]/(float)rsqcount_kremer[i];
        printf("%d %f %d %f %f\n", i, rsq[i],  rsqcount[i],(float)rsqcount_kremer[i],  rsqavg_intdist[i]);
    }

    return maxnbondsperchain;
}

#endif