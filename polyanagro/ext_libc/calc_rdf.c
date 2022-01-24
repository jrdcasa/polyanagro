#ifndef __RDF_H
#define __RDF_H

#include <stdio.h>
#include <math.h>

/*#ifdef _OPENMP
  #include <omp.h>
  #define USED_OPENMP 1
#else
  #define USED_OPENMP 0
#endif
*/

// Maximum distances in the RDF histograms
float maxDist;
// Number of bins in the RDF histograms
int nbins;
// Array to store the total RDF
float *gr;

void c_setup_rdf(float max_box, float delta_r)  {

    // Setup some variables for the histogram
    maxDist = max_box/2.0;
    nbins = (int) ceil((maxDist+1.0)/(2.0*delta_r)) + 1;

    // Allocate arrays
    if ( ((gr  = (float*)  malloc(sizeof(float)      * (nbins)) ) == NULL) ) {
        fprintf(stderr, "c_setup_rdf: out of memory\n");
        exit(0);
    }

    // Initialize arrays
    for (int ibin=0; ibin<nbins; ibin++) {
        gr[ibin] = 0.0;
    }

    printf("%f %d\n", maxDist, nbins);
}

//void c_calc_rdf_equal(int nat_A, int nat_B, int* index_A, int *index_B, float* coords_A, float* coords_B, float* box)  {
//    for (int iat=0; iat<nat_A; iat++) {
//        for (int jat=iat+1; jat<nat_B; jat++) {
//
//        printf("%d %d\n", index_A[iat], index_B[jat]);
//
//    }}
//
//
//
//
//    printf("%d %d\n", nat_A, nat_B);
//    for (int iat=0; iat<nat_A; iat++) {
//        printf("%d %f %f %f\n", iat, coords_A[iat+0], coords_A[iat+1], coords_A[iat+2]);
//    }
//    printf("*****\n");
//
//    for (int iat=0; iat<nat_B; iat++) {
//        printf("%d %f %f %f\n", iat, coords_B[iat+0], coords_B[iat+1], coords_B[iat+2]);
//    }
//
//    printf("*****\n");
//
//    printf("%f %f %f\n", box[0], box[1], box[2]);
//    printf("=============================================\n");


void c_calc_rdf_equal(int nat_A, long* index_A, float* coords_A, float* box)  {

     // If sets A and B contain the same number of atoms, only the half of interactions
     // have to be calculated
     //
     //  nat_A (int)      : Number of atoms in set A (and B)
     //  index_A (long)   : Array of Atom indices in the set,  index_A[nat_A]
     //  coords_A (float) : Coordinates in the set,  coords_A[nat_A,3]
     //  box (float)      : Box dimension,  box[3]


    float* xyz;
    float hbox[3];
    float tmp;
    float* rij;

    // Allocate arrays
    if ( ((rij  = (float*)  malloc(sizeof(float)      * (nat_A*nat_A)) ) == NULL) ||
         ((xyz  = (float*)  malloc(sizeof(float)      * (3)) ) == NULL))  {
        fprintf(stderr, "c_calc_rdf_equal: out of memory\n");
        exit(0);
    }

    // Half box
    hbox[0] = box[0]/2.0;
    hbox[1] = box[1]/2.0;
    hbox[2] = box[2]/2.0;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    // Main loop for distances
    for (int iat=0; iat<nat_A; iat++) {
        for (int jat=iat+1; jat<nat_A; jat++) {
            for (int k=0; k<3; k++){
                xyz[k]  = coords_A[iat+k]-coords_A[jat+k];
                if (fabs(tmp) > hbox[k]) {
                    if(tmp>0.0) {
                        xyz[k] = tmp - box[k];
                    } else {
                        xyz[k] = tmp + box[k];
                    }
                }
            }
            rij[iat*nat_A+jat] = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
        }
    }




//    for (int iat=0; iat<nat_A; iat++) {
//        printf("%d %f %f %f\n", iat, coords_A[iat+0], coords_A[iat+1], coords_A[iat+2]);
//    }
//    printf("*****\n");
//
//    printf("*****\n");
//
//    printf("%f %f %f\n", box[0], box[1], box[2]);
//    printf("=============================================\n");

}

#endif