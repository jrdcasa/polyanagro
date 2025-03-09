#ifndef __RDF_H
#define __RDF_H

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
  #include <omp.h>
  #define USED_OPENMP 1
#else
  #define USED_OPENMP 0
#endif


#define MAX(a, b) ((a)>(b)) ? (a) : (b)
#define MIN(a, b) ((a)<(b)) ? (a) : (b)
#define ANINT(a) (double) (int)  ( ((a) < 0.0) ? ((a) - 0.5) : ((a) + 0.5) )

// Maximum distances in the RDF histograms
float maxDist;
// Number of bins in the RDF histograms
int nbins;

void c_setup_rdf(int nbins, float* total_rdf)  {

    // float* total_rdf array to store the total RDF

    // Allocate arrays
    if ( ((total_rdf  = (float*)  malloc(sizeof(float)      * (nbins)) ) == NULL) ) {
        fprintf(stderr, "c_setup_rdf: out of memory\n");
        exit(0);
    }

    // Initialize arrays
    for (int ibin=0; ibin<nbins; ibin++) {
        total_rdf[ibin] = 0.0;
    }

}

float distance_minimum_image(float x1, float y1, float z1, float x2, float y2, float z2,
                             float* box, float* boxi) {

    float delx, dely, delz;
    float rsq, r;

    delx = x1 - x2;
    dely = y1 - y2;
    delz = z1 - z2;

    delx = delx - box[0]*ANINT(boxi[0]*delx);
    dely = dely - box[1]*ANINT(boxi[1]*dely);
    delz = delz - box[2]*ANINT(boxi[2]*delz);

//    delx = delx - box[0]*round(delx/box[0]);
//    dely = dely - box[1]*round(dely/box[1]);
//    delz = delz - box[2]*round(delz/box[2]);

    //delx = delx - round(delx);
    //dely = dely - round(dely);
    //delz = delz - round(delz);

    rsq = delx*delx+dely*dely+delz*delz;
    r = sqrt(rsq);

    return r;

}


void c_rdf_hist(int natoms, int nat_A, int nat_B, int nbins, float delta_r, float cutoff,
                int* atindexA, int* atindexB,
                float* coords, float* box, int* iatch,
                int* hist_total_rdf, int* hist_intra_rdf,
                int* hist_inter_rdf, int* hist_self_rdf, int* npairs_rdf_total,
                int* npairs_rdf_intra,
                int* npairs_rdf_inter) {


    /*
    Accumulate the histograms for each step
    */

    float dij;
    int dim2 = 3;
    float xA, yA, zA;
    float xB, yB, zB;
    int k;
    int iat;
    int jat;
    int ich, jch;
    float cutoffi = 1.0 / cutoff;
    float boxi[3];
    float dmax;
    int nocut = 0;
    // Mask to skip unnecessary distance calculations.
    // For example, avoid redundant pairs like [0, 26] and [26, 0].
    int mask[nat_A][nat_B];

    // Initialize all elements to zero.
    for (int i = 0; i < nat_A; i++) {
        for (int j = 0; j < nat_B; j++) {
            mask[i][j] = 0;
        }
    }

    // Determine the maximum allowed cutoff value.
    // If the cutoff exceeds half the box size,
    // set it to the smaller half-box value.
    dmax   = MIN(cutoff, box[0] * 0.5);
    dmax   = MIN(dmax, box[1] * 0.5);
    dmax   = MIN(dmax, box[2] * 0.5);
    if (nocut == 1) dmax = cutoff;

    for (int i = 0; i<dim2; i++) {
        boxi[i]  = 1.0 / box[i];
    }

    // Atoms in setA. Coordinates are converted to box=1 units
    for (int idx=0; idx<nat_A; idx++) {
          int count_valid = 0;  // Track number of valid pairs for this atom A_i
          iat = atindexA[idx];
          ich = iatch[idx];
          xA = coords[iat*dim2+0];
          yA = coords[iat*dim2+1];
          zA = coords[iat*dim2+2];
        // Atoms in setB. Coordinates are converted to box=1 units
        for (int jdx=0; jdx<nat_B; jdx++) {
            jat = atindexB[jdx];
            jch = iatch[jdx];
            if (iat == jat) continue;
            // Avoid double counting
            if (mask[iat][jat] != 0 || mask[jat][iat] != 0) continue;
            xB = coords[jat*dim2+0];
            yB = coords[jat*dim2+1];
            zB = coords[jat*dim2+2];
            dij = distance_minimum_image(xA, yA, zA, xB, yB, zB, box, boxi);
            mask[iat][jat] = 1;
            mask[jat][iat] = 1;
            k = (int) (nbins*dij*cutoffi);
            //printf("%d %d %f %d %f\n",iat, jat, dij, k, dmax);

            if (dij<=dmax) hist_total_rdf[k] = hist_total_rdf[k] + 1;
            if ((dij<=dmax) && (ich == jch)) hist_intra_rdf[k] = hist_intra_rdf[k] + 1;
            if ((dij<=dmax) && (ich != jch)) hist_inter_rdf[k] = hist_inter_rdf[k] + 1;

            npairs_rdf_total[0] += 1;
            if (ich == jch) npairs_rdf_intra[0] += 1;
            if (ich != jch) npairs_rdf_inter[0] += 1;

        }
    }
}


void c_rdf_gr(int nframes, int nat_A, int nat_B, int nbins,
              float cutoff, float volume_avg, int* hist_total_rdf,
              int* hist_intra_rdf,
              int* hist_inter_rdf, int* hist_self_rdf,
              float* total_rdf, int npairs_rdf_total,
              float* intra_rdf, int npairs_rdf_intra,
              float* inter_rdf, int npairs_rdf_inter) {

    float pi = 3.14159265359;
    float bin_width;
    double shell_volume;
    double radii[nbins];

    bin_width = cutoff / nbins;

    for (int i=0; i<nbins; i++) {

        radii[i] = (i + 0.5) * bin_width;
        shell_volume = 4.0 * pi * radii[i] * radii[i] * bin_width;

        total_rdf[i] = (volume_avg / (npairs_rdf_total)) *
                       (hist_total_rdf[i]/shell_volume);
        intra_rdf[i] = (volume_avg / (npairs_rdf_intra)) *
                       (hist_intra_rdf[i]/shell_volume);
        inter_rdf[i] = (volume_avg / (npairs_rdf_inter)) *
                       (hist_inter_rdf[i]/shell_volume);

    }

    for (int i=0; i<nbins; i++) {
        total_rdf[i] = total_rdf[i]/nframes;
        intra_rdf[i] = intra_rdf[i]/nframes;
        inter_rdf[i] = inter_rdf[i]/nframes;
    }

}


void c_rdf_hist_openmp(int natoms, int nat_A, int nat_B, int nbins, float delta_r, float cutoff,
                int* atindexA, int* atindexB,
                float* coords, float* box, int* iatch,
                int* hist_total_rdf, int* hist_intra_rdf,
                int* hist_inter_rdf, int* hist_self_rdf, int* npairs_rdf_total,
                int* npairs_rdf_intra,
                int* npairs_rdf_inter) {

    float dij;
    int dim2 = 3;
    float xA, yA, zA;
    float xB, yB, zB;
    int k, iat, jat, ich, jch;
    float cutoffi = 1.0 / cutoff;
    float boxi[3];
    float dmax;
    int nocut = 0;

    dmax = MIN(cutoff, box[0] * 0.5);
    dmax = MIN(dmax, box[1] * 0.5);
    dmax = MIN(dmax, box[2] * 0.5);
    if (nocut == 1) dmax = cutoff;

    for (int i = 0; i < dim2; i++) {
        boxi[i] = 1.0 / box[i];
    }

    // OpenMP parallelization
    #pragma omp parallel for private(iat, ich, xA, yA, zA, jat, jch, xB, yB, zB, dij, k) \
                             reduction(+: hist_total_rdf[:nbins], hist_intra_rdf[:nbins], hist_inter_rdf[:nbins], \
                                        npairs_rdf_total[0], npairs_rdf_intra[0], npairs_rdf_inter[0])
    for (int idx = 0; idx < nat_A; idx++) {
        iat = atindexA[idx];
        ich = iatch[idx];
        xA = coords[iat * dim2 + 0];
        yA = coords[iat * dim2 + 1];
        zA = coords[iat * dim2 + 2];

        for (int jdx = 0; jdx < nat_B; jdx++) {
            jat = atindexB[jdx];
            jch = iatch[jdx];

            if (iat == jat) continue;

            xB = coords[jat * dim2 + 0];
            yB = coords[jat * dim2 + 1];
            zB = coords[jat * dim2 + 2];

            dij = distance_minimum_image(xA, yA, zA, xB, yB, zB, box, boxi);

            k = (int)(nbins * dij * cutoffi);
            if (dij <= dmax) hist_total_rdf[k]++;
            if ((dij <= dmax) && (ich == jch)) hist_intra_rdf[k]++;
            if ((dij <= dmax) && (ich != jch)) hist_inter_rdf[k]++;

            npairs_rdf_total[0]++;
            if (ich == jch) npairs_rdf_intra[0]++;
            if (ich != jch) npairs_rdf_inter[0]++;
        }
    }
}

void c_rdf_hist_excl_openmp(int natoms, int nat_A, int nat_B, int nbins, float delta_r, float cutoff,
                int* atindexA, int* atindexB, int* excl_array,
                float* coords, float* box, int* iatch,
                int* hist_total_rdf, int* hist_intra_rdf,
                int* hist_inter_rdf, int* hist_self_rdf, int* npairs_rdf_total,
                int* npairs_rdf_intra,
                int* npairs_rdf_inter) {

    float dij;
    int dim2 = 3;
    int dimarray2 = 40;
    float xA, yA, zA;
    float xB, yB, zB;
    int k, iat, jat, kat, ich, jch, isexcl;
    float cutoffi = 1.0 / cutoff;
    float boxi[3];
    float dmax;
    int nocut = 0;

    dmax = MIN(cutoff, box[0] * 0.5);
    dmax = MIN(dmax, box[1] * 0.5);
    dmax = MIN(dmax, box[2] * 0.5);
    if (nocut == 1) dmax = cutoff;

    for (int i = 0; i < dim2; i++) {
        boxi[i] = 1.0 / box[i];
    }

    // OpenMP parallelization
//    #pragma omp parallel for private(iat, ich, xA, yA, zA, jat, jch, xB, yB, zB, dij, k, isexcl, kat) \
//                             reduction(+: hist_total_rdf[:nbins], hist_intra_rdf[:nbins], hist_inter_rdf[:nbins], \
//                                        npairs_rdf_total[0], npairs_rdf_intra[0], npairs_rdf_inter[0])

    for (int idx = 0; idx < nat_A; idx++) {
        iat = atindexA[idx];
//        printf( "******* %d \n", iat);
        ich = iatch[idx];
        isexcl = 0;
        xA = coords[iat * dim2 + 0];
        yA = coords[iat * dim2 + 1];
        zA = coords[iat * dim2 + 2];

        for (int jdx = 0; jdx < nat_B; jdx++) {
            jat = atindexB[jdx];
            jch = iatch[jdx];

            if (iat == jat) continue;

            //printf("\n%d ::: ", iat);
            // Check if jat should be excluded
            for (int kdx = 0; kdx < dimarray2; kdx++) {
                kat = excl_array[iat * dimarray2 + kdx];
//                if (iat == 10) {
//                    if ( kat != -1 ) {
//                        printf("======== %d %d %d \n", iat, kdx, kat);
//                    }
//                }
                if (jat == kat) {
                    isexcl = 1;
                    break;
                }

                //printf("%d ", kat);
            }
            if (isexcl == 1) {
//                printf("\n* %d %d %d *\n", iat, jat, kat);
                continue;
            }

            xB = coords[jat * dim2 + 0];
            yB = coords[jat * dim2 + 1];
            zB = coords[jat * dim2 + 2];

            dij = distance_minimum_image(xA, yA, zA, xB, yB, zB, box, boxi);

            k = (int)(nbins * dij * cutoffi);
            if (dij <= dmax) hist_total_rdf[k]++;
            if ((dij <= dmax) && (ich == jch)) hist_intra_rdf[k]++;
            if ((dij <= dmax) && (ich != jch)) hist_inter_rdf[k]++;

            npairs_rdf_total[0]++;
            if (ich == jch) npairs_rdf_intra[0]++;
            if (ich != jch) npairs_rdf_inter[0]++;
        }
    }
}

#endif