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

float distance_minimum_image(float x1, float y1, float z1, float x2, float y2, float z2, float* box) {

    float delx, dely, delz;
    float rsq, r;

    delx = x1 - x2;
    dely = y1 - y2;
    delz = z1 - z2;

    //delx = delx - box[0]*round(delx/box[0]);
    //dely = dely - box[1]*round(dely/box[1]);
    //delz = delz - box[2]*round(delz/box[2]);

    delx = delx - round(delx);
    dely = dely - round(dely);
    delz = delz - round(delz);

    rsq = delx*delx+dely*dely+delz*delz;
    r = sqrt(rsq);

    return r;

}

void c_rdf_hist(int natoms, int nat_A, int nat_B, int nbins, float delta_r,
                int* atindexA, int* atindexB,
                float* coords, float* box, int* hist_rdf) {

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
    float blmax = -1.0;
    float dr_reduced;

    for (int i = 0; i<dim2; i++) {
        if (box[i] > blmax) blmax = box[i];
    }
    dr_reduced = delta_r/blmax;

    // Atoms in setA. Coordinates are converted to box=1 units
    for (int idx=0; idx<nat_A; idx++) {
          iat = atindexA[idx];
          xA = coords[iat*dim2+0]/box[0];
          yA = coords[iat*dim2+1]/box[1];
          zA = coords[iat*dim2+2]/box[2];
        // Atoms in setB
        for (int jdx=0; jdx<nat_B; jdx++) {
            jat = atindexB[jdx];
            if (iat == jat) continue;
            xB = coords[jat*dim2+0]/box[0];
            yB = coords[jat*dim2+1]/box[1];
            zB = coords[jat*dim2+2]/box[2];
            dij = distance_minimum_image(xA, yA, zA, xB, yB, zB, box);

            k = floor(dij / dr_reduced) + 1;
            printf("%d %d %f %d %f\n",iat, jat, dij, k, dr_reduced);
            if (k<=nbins) hist_rdf[k] = hist_rdf[k] + 1;
        }
    }

}

void c_rdf_gr() {

    float prefact, cons, vol;
    float pi = 3.14159265359;

    cons = 4.0 * pi / 3.0;
    vol = 1.0;

    // Calculate gr for an ideal gas of the same density
    prefact = cons / vol;

    printf("HJKKKJJJJ\n");
}

#endif