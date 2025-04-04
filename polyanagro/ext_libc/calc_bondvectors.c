#ifndef __BONDVECTORS_H
#define __BONDVECTORS_H

#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
  #include <omp.h>
  #define USED_OPENMP 1
#else
  #define USED_OPENMP 0
#endif

double *ODF, *P2, *ODFl, *P2l, *pxy, *pxy2;
double *pxx, *pyy, *pzz, *px4, *py4, *pz4, *pxy4, *pxz4, *pyz4;
int    *N;
int ngroups;
int nbondsperchain;

static void c_setup_odf_intra(int maxnbondsperch) {

    // Initialize arrays

    // Number of bonds per chain. TODO: Test if it works with different length chains
    int nbondsperchain = maxnbondsperch; //nbin

    if ( ((N  = (int*)  malloc(sizeof(int)      * (nbondsperchain+1)) ) == NULL) ||
       ((ODFl = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
       ((P2l  = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
       ((ODF  = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
       ((P2   = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ) {
        fprintf(stderr, "c_odf_intra: out of memory\n");
        exit(0);
    }

    if ( ((pxx = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
           ((pyy = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
           ((pzz = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
           ((px4 = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
           ((py4 = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
           ((pz4 = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
           ((pxy4 = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
           ((pxz4 = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
           ((pyz4 = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
           ((pxy = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ||
           ((pxy2 = (double*) malloc(sizeof(double) * (nbondsperchain+1)) ) == NULL) ) {
        fprintf(stderr, "c_odf_intra: out of memory\n");
        exit(0);
    }

    for (int i=0;i<=nbondsperchain;i++) {
        N[i] = 0;
        ODFl[i] = 0.0;
        P2l[i] = 0.0;
        ODF[i] = 0.0;
        P2[i] = 0.0;
        pxx[i] = 0.0; pyy[i] = 0.0; pzz[i] = 0.0;
        px4[i] = 0.0; py4[i] = 0.0; pz4[i] = 0.0;
        pxy4[i] = 0.0; pxz4[i] = 0.0; pyz4[i] = 0.0;
        pxy[i] = 0.0; pxy2[i] = 0.0;

    }

}


static void unit_vector(float c1[3], float c2[3], float dist[4]) {

   // Calculate the distance between at1 and at2 (float[3])
   // It returns a float[4] containing the distance and the three components of the
   // unit vector.

    float dx, dy, dz;

    // Distance between the two points
    dx = c2[0]-c1[0];
    dy = c2[1]-c1[1];
    dz = c2[2]-c1[2];
    dist[0] = sqrt((dx*dx)+(dy*dy)+(dz*dz));

    // Normalize the distance vector
    dist[1] = dx/dist[0];
    dist[2] = dy/dist[0];
    dist[3] = dz/dist[0];


}

static float c_unit_bond_vectors(int natoms, int nchains, int nbonds, int maxnbondsperch,
                                int* bonds, float* coords,
                                int* iatch, float* uux, float* uuy, float* uuz)
{

    /*
        Calculate the unit vectors. Calculate Cn and data for persistence of lenth (bond-bond correlation decay)

    natoms [int]         --> Total number of atoms in the system
    nchains [int]        --> Total number of chains in the system
    nbonds [int]         --> Number of backbone bonds
    maxnbondsperch [int] --> Maximum number of backbone bonds in a chain
    bonds [nbonds,2]     --> Each row is a bond. Example bonds[2,0] = 2 and bonds[2,1] =3 means a bond 2-3
    coords[natoms,3]     --> Coordinates of the ith-atom
    iatch[natoms]        --> iatch[i]=j --> The i-th atom is in the j-ch chain
    uux[nchains,maxnbondsperch] --> x-component of the unit vectors per chain
    uuy[nchains,maxnbondsperch] --> y-component of the unit vectors per chain
    uuz[nchains,maxnbondsperch] --> z-component of the unit vectors per chain

    Calculate several things:

        1) corr[0:maxnbondsperch-2] = (1/maxnbondsperch) sum(i=0, maxnbondsperch-1)sum(j=i+1, maxnbondsperch) <ui*uj>.
           The average <> is taken into account all chains in the system.
           The index of corr array is (j-i-1)
        2) s = sum(i=0, maxnbondsperch-2])corr[i]
        3) Cn using the formula from (2.9) to (2.11) in Polymer Physics (M.Rubinstein-RH Colby)
           A similar formula (Eq.8) can be found in C.A. Helfer, W.L. Mattice / Polymer 46 (2005) 4361–4367
           I have check that both formulas give the same results
           I have implemented the second one.

    */

    int dim1 = 2;
    int dim2 = 3;
    int at1, at2, ich1, ich2;
    float uu_vector[4]; //= malloc(4*sizeof(float));
    float c1[3], c2[3];
    int ilocalbond = -1;
    int ich_prev;
    float cn[nchains];
    float cn_avg;
    float corr;
    float tmp;

    // Bonds =============== Calculate the uux, uuy and uuz vectors
    ich_prev = 0;
    for (int ibond=0; ibond<nbonds; ibond++) {
        at1 = bonds[ibond*dim1+0];
        at2 = bonds[ibond*dim1+1];
        ich1 = iatch[at1];
        ich2 = iatch[at2];

        if (ich1 != ich2) {
            printf("Something is wrong in c_unit_bond_vectors()!!!!\n");
            printf("Expected %d and %d be in the same chains", at1, at2);
            printf("Chains %d and %d", ich1, ich2);
            continue;
        }

        if (ich1 != ich_prev) {
            ilocalbond=0;
            ich_prev = ich1;
        } else {
            ilocalbond+=1;
        }

        for (int i=0; i<dim2; i++) {
            c1[i] = coords[at1*dim2+i];
            c2[i] = coords[at2*dim2+i];
        }
        unit_vector(c1, c2, uu_vector);

        uux[ich1*maxnbondsperch+ilocalbond] = uu_vector[1];
        uuy[ich1*maxnbondsperch+ilocalbond] = uu_vector[2];
        uuz[ich1*maxnbondsperch+ilocalbond] = uu_vector[3];

//        printf("%d %f %f %f\n", at1, c1[0], c1[1], c1[2]);
//        printf("%d %f %f %f\n", at2, c2[0], c2[1], c2[2]);
//        printf("%d %d %d %d %f\n",at1, at2, ich1, ich2, uu_vector[0]);
//        printf("%f %f %f\n",uu_vector[1], uu_vector[2], uu_vector[3]);
//        printf("============\n");
    }

    // Cn calculation
    //  corr[j-i] = ui*uj
    cn_avg = 0.0;
    for (int ich=0; ich<nchains; ich++) {
        corr = 0.0;
        for (int i=0; i<maxnbondsperch-1; i++) {
            for (int j=i+1; j<maxnbondsperch; j++) {
                tmp =  uux[ich*maxnbondsperch+i]*uux[ich*maxnbondsperch+j]+
                       uuy[ich*maxnbondsperch+i]*uuy[ich*maxnbondsperch+j]+
                       uuz[ich*maxnbondsperch+i]*uuz[ich*maxnbondsperch+j];
                corr += tmp;

            }
        }

        cn[ich] = 1. +  2./maxnbondsperch*corr;
        //printf("%d %f %f\n",ich, corr, cn[ich]);
        cn_avg += cn[ich];
    }

    cn_avg = cn_avg/nchains;

    return cn_avg;

}

static float c_odf_intra(int iframe, int natoms, int nchains, int nbonds,
                         int maxnbondsperch,
                         int* bonds, float* coords,
                         int* iatch, float* uux, float* uuy, float* uuz)
{

    /*
        Calculate intra chain orientation correlation

    iframe [int]         --> Number of the current frame
    natoms [int]         --> Total number of atoms in the system
    nchains [int]        --> Total number of chains in the system
    nbonds [int]         --> Number of backbone bonds
    maxnbondsperch [int] --> Maximum number of backbone bonds in a chain
    bonds [nbonds,2]     --> Each row is a bond. Example bonds[2,0] = 2 and bonds[2,1] =3 means a bond 2-3
    coords[natoms,3]     --> Coordinates of the ith-atom
    iatch[natoms]        --> iatch[i]=j --> The i-th atom is in the j-ch chain
    uux[nchains,maxnbondsperch] --> x-component of the unit vectors per chain
    uuy[nchains,maxnbondsperch] --> y-component of the unit vectors per chain
    uuz[nchains,maxnbondsperch] --> z-component of the unit vectors per chain

    TODOOOOOOO


    */

    int dim1 = 2;
    int dim2 = 3;
    int at1, at2, ich1, ich2;
    float uu_vector[4]; //= malloc(4*sizeof(float));
    float c1[3], c2[3];
    int ilocalbond = -1;
    int ich_prev;
    int jj;
    double costhe, axx, ayy, azz;

    // Bonds =============== Calculate the uux, uuy and uuz vectors
    ich_prev = 0;
    for (int ibond=0; ibond<nbonds; ibond++) {
        at1 = bonds[ibond*dim1+0];
        at2 = bonds[ibond*dim1+1];
        ich1 = iatch[at1];
        ich2 = iatch[at2];

        if (ich1 != ich2) {
            printf("Something is wrong in c_odf_intra()!!!!\n");
            printf("Expected %d and %d be in the same chains", at1, at2);
            printf("Chains %d and %d", ich1, ich2);
            continue;
        }

        if (ich1 != ich_prev) {
            ilocalbond=0;
            ich_prev = ich1;
        } else {
            ilocalbond+=1;
        }

        for (int i=0; i<dim2; i++) {
            c1[i] = coords[at1*dim2+i];
            c2[i] = coords[at2*dim2+i];
        }

        unit_vector(c1, c2, uu_vector);

        uux[ich1*maxnbondsperch+ilocalbond] = uu_vector[1];
        uuy[ich1*maxnbondsperch+ilocalbond] = uu_vector[2];
        uuz[ich1*maxnbondsperch+ilocalbond] = uu_vector[3];

//        printf("%d %f %f %f\n", at1, c1[0], c1[1], c1[2]);
//        printf("%d %f %f %f\n", at2, c2[0], c2[1], c2[2]);
//        printf("%d %d %d %d %f\n",at1, at2, ich1, ich2, uu_vector[0]);
//        printf("%f %f %f\n",uu_vector[1], uu_vector[2], uu_vector[3]);
//        printf("============\n");
    }

    // Number of bonds per chain. TODO: Test if it works with different lenght chains
    nbondsperchain = maxnbondsperch; //nbin
    // Total number of bonds
    ngroups = nbondsperchain*nchains;

    // Write test file for DEBUG Reasons====================================
    // This file can be used directly in the original odfintra.c (Meyer) ===
//    FILE *fp;
//    fp = fopen("./intra_fake.inp", "a");
//    fprintf(fp, "10.0 10.0 10.0\n");
//    for (int igr=0; igr<nbondsperchain; igr++)  {
//        fprintf(fp, "0.000 0.00 0.000 %f %f %f\n",uux[igr], uuy[igr] , uuz[igr]  );
//    }
//    fclose(fp);

    /*Example: 3 chains

                 j
    i: ()-()-(ii)-()-()-(jj)-()-()-()   nbondsperchain=8

       ()-()-(  )-()-()-(  )-()-()-()

       ()-()-(  )-()-()-(  )-()-()-()

                                     ngroups = 8*3=24
    */
    /* for reasons of precision, sum up every configuration independently */
    for (int i=0;i<=nbondsperchain;i++) {
        ODFl[i] = 0.0; P2l[i] = 0.0;
    }

    for (int i = 0; i < ngroups/nbondsperchain; i++) {


      for (int j=1;j<nbondsperchain;j++) {
        //ii bond vector initial
        //ii = i*nbondsperchain+j;
        for (int ii = i*nbondsperchain; ii<(i+1)*nbondsperchain-j; ii++) {
           //jj bond vector end
          jj = ii+j;
          ++N[j];
          axx = uux[ii] * uux[jj];  // u(0)u(n)
          ayy = uuy[ii] * uuy[jj];
          azz = uuz[ii] * uuz[jj];
          costhe = axx+ayy+azz;
          ODFl[j] += costhe;
          P2l[j] += costhe*costhe;
          pxx[j] += axx;
          pyy[j] += ayy;
          pzz[j] += azz;
          px4[j] += axx*axx;
          py4[j] += ayy*ayy;
          pz4[j] += azz*azz;
          pxy4[j] += axx*ayy;
          pxz4[j] += axx*azz;
          pyz4[j] += ayy*azz;
          costhe = uux[ii] * uuy[jj] + uuy[ii] * uuz[jj] + uux[ii] * uuz[jj]
                 + uuy[ii] * uux[jj] + uuz[ii] * uuy[jj] + uuz[ii] * uux[jj];
          pxy[j]  += costhe;
          pxy2[j] += costhe*costhe;  /* probably not a useful quantity */

        }
      }
    }

    for (int i=0; i<=nbondsperchain; i++) {
      ODF[i] += ODFl[i];
      P2[i] += P2l[i];
    }

    return 0;
}

void c_avg_write_odf_intra(int nframes, int maxnbondsperch, char* filename) {

    // Write the ODF INTRA results
    int frame = nframes;

    /* normalise */
    for (int i = 1; i < maxnbondsperch; i++) {
        if (N[i] > 0) {
          ODF[i] /=  N[i];
          P2[i]  = 1.5*P2[i]/ N[i] - 0.5;
          pxx[i] /= N[i];
          pyy[i] /= N[i];
          pzz[i] /= N[i];
          px4[i] /= N[i];
          py4[i] /= N[i];
          pz4[i] /= N[i];
          pxy4[i] /= N[i];
          pxz4[i] /= N[i];
          pyz4[i] /= N[i];
          pxy[i] /= N[i]*0.5;
          pxy2[i] /= N[i]*0.5;
      }
    }

    FILE *fp;
    fp = fopen(filename, "w");
    fprintf(fp, "#%s\n", filename);
    fprintf(fp,"# odfintra Version: intra-chain ODF with P2 and pxy\n");
    fprintf(fp,"# %d vectors per frame, %d vectors per chain, %d chains\n",ngroups, nbondsperchain,ngroups/nbondsperchain);
    fprintf(fp,"# P1(n) = <u(0)u(n)> = <costhe>  \n");
    fprintf(fp,"# P2(n) = 1/2*<3*(u(0)u(n))^2 - 1> = 1/2<3*costhe^2 - 1> 2nd Legendre polynomial \n");
    fprintf(fp,"# P2 values +1: Parallel, -0.5: orthogonal, 0: No correlation ");
    fprintf(fp,"# u(0) u(n) are unit vectors at distance n bonds.  \n");
    fprintf(fp,"# n  P1  P2  pxy  pxy2  xx yy zz  x4 y4 z4  xxyy xxzz yyzz\n");
    fprintf(fp, "# frames %d\n 0 1.0 1.0 1.0 1.0  1.0 1.0 1.0  1.0 1.0 1.0  1.0 1.0 1.0\n",frame);
    for (int i = 1; i < maxnbondsperch; i++) {
        fprintf(fp, "%d %.8lg %.8lg %lg %lg",i,ODF[i],P2[i],pxy[i], pxy2[i]);
        fprintf(fp, " %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", pxx[i],pyy[i],pzz[i],
           px4[i],py4[i],pz4[i],pxy4[i],pxz4[i],pyz4[i]);
    }

    fclose(fp);

}

static void c_bond_bond_orientation(int natoms, int nchains, int nbonds,
                                     int maxnbondsperch, int* bonds, float* coords,
                                     int* iatch, float *cbb)
{

    /*
        Calculate the unit vectors. Calculate Cn and data for persistence of lenth (bond-bond correlation decay)

    natoms [int]         --> Total number of atoms in the system
    nchains [int]        --> Total number of chains in the system
    nbonds [int]         --> Number of backbone bonds
    maxnbondsperch [int] --> Maximum number of backbone bonds in a chain
    bonds [nbonds,2]     --> Each row is a bond. Example bonds[2,0] = 2 and bonds[2,1] =3 means a bond 2-3
    coords[natoms,3]     --> Coordinates of the ith-atom
    iatch[natoms]        --> iatch[i]=j --> The i-th atom is in the j-ch chain

    Calculate several things:

        1) The bond-bond correlation decay is calculates using the formula (3) in
            Salerno et al. J. Chem. Theory Comput. 2018, 14, 2219-2229
            cbb(m) = <bi*bi+m/|bi||bi+m|> = <ui*ui+m>
    */

    int dim1 = 2;
    int dim2 = 3;
    int at1, at2, ich1, ich2;
    float uu_vector[4]; //= malloc(4*sizeof(float));
    float c1[3], c2[3];
    int ilocalbond = -1;
    int ich_prev;
    int m;
    float tmp;
    float uux[nchains*maxnbondsperch];
    float uuy[nchains*maxnbondsperch];
    float uuz[nchains*maxnbondsperch];
    float cbb_nelem[maxnbondsperch-1];

    // Bonds =============== Calculate the uux, uuy and uuz vectors
    ich_prev = 0;
    for (int ibond=0; ibond<nbonds; ibond++) {
        at1 = bonds[ibond*dim1+0];
        at2 = bonds[ibond*dim1+1];
        ich1 = iatch[at1];
        ich2 = iatch[at2];

        if (ich1 != ich2) {
            printf("Something is wrong in c_unit_bond_vectors()!!!!\n");
            printf("Expected %d and %d be in the same chains", at1, at2);
            printf("Chains %d and %d", ich1, ich2);
            continue;
        }

        if (ich1 != ich_prev) {
            ilocalbond=0;
            ich_prev = ich1;
        } else {
            ilocalbond+=1;
        }

        for (int i=0; i<dim2; i++) {
            c1[i] = coords[at1*dim2+i];
            c2[i] = coords[at2*dim2+i];
        }
        unit_vector(c2, c1, uu_vector);

        uux[ich1*maxnbondsperch+ilocalbond] = uu_vector[1];
        uuy[ich1*maxnbondsperch+ilocalbond] = uu_vector[2];
        uuz[ich1*maxnbondsperch+ilocalbond] = uu_vector[3];

    }

    // Bond-to-bond correlation
    for (int ibond=0; ibond<maxnbondsperch-1; ibond++) {
        cbb[ibond] = 0.0;
        cbb_nelem[ibond] = 0;
    }

    for (int ich=0; ich<nchains; ich++) {
        for (int i=0; i<maxnbondsperch-1; i++) {
            for (int j=i+1; j<maxnbondsperch; j++) {
                tmp =  uux[ich*maxnbondsperch+i]*uux[ich*maxnbondsperch+j]+
                       uuy[ich*maxnbondsperch+i]*uuy[ich*maxnbondsperch+j]+
                       uuz[ich*maxnbondsperch+i]*uuz[ich*maxnbondsperch+j];
                m = j-i;
                cbb[m] += tmp;
                cbb_nelem[m] += 1;
            }
        }
    }

    for (int m=1; m<maxnbondsperch-1; m++) {
//        printf("%d %f %f\n", m, cbb[m], cbb_nelem[m] );
        cbb[m] = cbb[m]/cbb_nelem[m];
    }

}

#endif