#ifndef __MSD_H
#define __MSD_H

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
  #include <omp.h>
  #define USED_OPENMP 1
#else
  #define USED_OPENMP 0
#endif

// ==============================================================================
void shiftcm(double *positions, int iframe, int natoms,
             double *xraw, double *yraw, double *zraw )
{
     /* subtract total center of mass from coordinates*/
      int i;
      double xcm, ycm, zcm;
      int dim = 3;

      xcm = ycm = zcm = 0.0;
      for (i=0; i<natoms; i++) {
        xcm += positions[iframe*natoms*dim + i*dim + 0];
        ycm += positions[iframe*natoms*dim + i*dim + 1];
        zcm += positions[iframe*natoms*dim + i*dim + 2];
      }
      xcm /= natoms;
      ycm /= natoms;
      zcm /= natoms;

      for (i=0; i<natoms; i++) {
        xraw[i] = positions[iframe*natoms*dim + i*dim + 0] - xcm;
        yraw[i] = positions[iframe*natoms*dim + i*dim + 1] - ycm;
        zraw[i] = positions[iframe*natoms*dim + i*dim + 2] - zcm;
      }


} /* shiftcm() */

// ==============================================================================
void copyfrm(double *positions, int iframe, int natoms,
             double *xraw, double *yraw, double *zraw)
{

  int i;
  int dim = 3;
  for (i=0; i<natoms; i++) {
    xraw[i] = positions[iframe*natoms*dim + i*dim + 0];
    yraw[i] = positions[iframe*natoms*dim + i*dim + 1];
    zraw[i] = positions[iframe*natoms*dim + i*dim + 2];
  }

} /* copyfrm() */

// ==============================================================================
void calccm(int ns, double *xraw, double *yraw, double *zraw,
            int nmol, double *xcm, double *ycm, double *zcm, int g2flag)
{
    /* calculate center of mass of molecules */
    int i, jj, nsit;
    double xx,yy,zz;

    nsit = ns/nmol;  /* number of sites per molecule */
    for (i=0; i<nmol; i++) {
        xcm[i] = 0.0;
        ycm[i] = 0.0;
        zcm[i] = 0.0;
    }

    for (i=0; i<ns; i++) {
        jj = i/nsit;
        xcm[jj] += xraw[i];
        ycm[jj] += yraw[i];
        zcm[jj] += zraw[i];
    }

    for (i=0; i<nmol; i++) {
        xcm[i] /= (double) nsit;
        ycm[i] /= (double) nsit;
        zcm[i] /= (double) nsit;
    }

    if (g2flag==1) {
        /* msdmolg2:  subtract individual cm of each molecule: */
        for (i=0; i<ns; i++) {
            jj = i/nsit;
            xraw[i] -= xcm[jj];
            yraw[i] -= ycm[jj];
            zraw[i] -= zcm[jj];
        }
    }

} /* calccm() */


// ==============================================================================
void c_msd_all (double *positions, int nframes, int nmol,
                int natoms, double timestep, int weedf,
                const char *filename) {

    /*
        positions(nframes, natoms, 3) -> Positions of all atoms in the trayectory
        nmol -> Number of molecules
        natoms -> Number of atoms
        weedf -> Number of starting points of independient time series
    */

    // Declaration
    int i, j, k, ii, jj, kk, ij, nsit;
    int tframes, frames;
    double aa, xx, yy, zz;
    double *time, *times;
    double *dx, *dy, *dz;                             /* actual positions (2d array atom,frm) */
    double *x2, *y2, *z2, *r2;                        /* msd of components and total 2nd+4th moment */
    double *x2av, *y2av, *z2av, *r2av, *r4av, *r6av;  /* total averages */
    double *x4av, *y4av, *z4av;                       /* xyz */
    double *r2a, *r2cm, *r4cm;                        /* partial averages */
    double dtime;
    int normcorrect = 0;
    clock_t start, end;

    // Declaration and initialization
    int np = natoms + nmol;  //Allocate the atoms and com in the same array
    int ns = natoms;

    // Default values can be modified in the input command
    //option -m
    int meanf = 99999;
    int nfrm = nframes;
    // option -abs
    int noshiftcm = 1;          /* flag: 0=subtract total cm; 1=do nothing */
    // option -nosym
    int nosym = 0;  /* flag =0: apply chain symmetry (2 ends) 1: no symmetry */
    // option: -g2 -> calculate in center-of-mass
    // coordinates of each individual chain for g2 (monomerresolved)
    int g2flag = 0;             /* flag =0: standard MSD in cm-frame of simulation box */
                                /*       1: calculate in cm-frame of each chain (g2) */
    // option "-nocorrect"
    int timecorrect = 1;        /* flag=0: calculate time differences to first frame */
                                /*     =1: correct for jumps between trajectories */
    int innerav = 16;   /* number of inner monomers for averaging */
    // option --xz
    int filmoutput = 0;

    /* maxmean: maximum number of frames, value needed for array
       allocation - real number of frames is (max-weedf)*weedf.
       default: 600  needed: (tframes/weedf)+weedf
    nfrm =
    meanf =
    */
    if (nfrm < meanf) meanf = nfrm;
    /* restrict max number of frames to size of allocated arrays. */
//    fprintf(stderr,"msdmol info: nfrm=%d, trajectories may contain total of %d frames [with g6].\n",
//	        nfrm,(nfrm-weedf)*weedf);
    nsit = ns/nmol;
    if (innerav > nsit/8) innerav = (nsit-1)/8 +1;  /* will lead to max 1/4 of chain */
    if (( (time = (double*) malloc(sizeof(double) * nfrm)) == NULL) ||
        ( (times = (double*) malloc(sizeof(double) * weedf)) == NULL) ||
        ( (dx = (double*) malloc(sizeof(double) * np*nfrm)) == NULL) ||
        ( (dy = (double*) malloc(sizeof(double) * np*nfrm)) == NULL) ||
        ( (dz = (double*) malloc(sizeof(double) * np*nfrm)) == NULL) ) {
        fprintf (stderr, "%s: out of memory (dxyz)\n", "c_msd_all");
        exit(2);
    }

    if (( (x2 = (double*) malloc(sizeof(double) * np*nfrm)) == NULL) ||
        ( (y2 = (double*) malloc(sizeof(double) * np*nfrm)) == NULL) ||
        ( (z2 = (double*) malloc(sizeof(double) * np*nfrm)) == NULL) ||
        ( (r2 = (double*) malloc(sizeof(double) * np*nfrm)) == NULL)  ) {
        fprintf (stderr, "%s: out of memory (msd)\n", "c_msd_all");
        exit(2);
    }

    if (( (x2av = (double*) malloc(sizeof(double) * nfrm)) == NULL) ||
        ( (y2av = (double*) malloc(sizeof(double) * nfrm)) == NULL) ||
        ( (z2av = (double*) malloc(sizeof(double) * nfrm)) == NULL) ||
        ( (r2av = (double*) malloc(sizeof(double) * nfrm)) == NULL) ||
        ( (r4av = (double*) malloc(sizeof(double) * nfrm)) == NULL) ||
        ( (x4av = (double*) malloc(sizeof(double) * nfrm)) == NULL) ||
        ( (y4av = (double*) malloc(sizeof(double) * nfrm)) == NULL) ||
        ( (z4av = (double*) malloc(sizeof(double) * nfrm)) == NULL) ||
        ( (r2a = (double*) malloc(sizeof(double) * nsit*nfrm)) == NULL) ) {
      fprintf (stderr, "%s: out of memory (average)\n", "c_msd_all");
      exit(2);
    }

    // Initialize
    for (i=0; i<natoms*nfrm; i++) {
        x2[i] = 0.0;
        y2[i] = 0.0;
        z2[i] = 0.0;
        r2[i] = 0.0;
    }

    tframes = 0;     /* Count the number of frames */
    frames = weedf;  /* Count the frames for averaging */

    start = clock();
    printf("\t\t[Step 1 of 3] \n");
    for (int iframe = 0; iframe<nframes; iframe++) {

        // Progress bar
        if (iframe % (nframes / 100 + 1) == 0 || iframe == nframes - 1) {
            printf("\r\t\t[Step 1 of 3] Progress: %3d%%", (iframe * 100) / nframes);
            fflush(stdout);
        }

        ii = tframes%weedf;  //Module

        if ((ii==0) && (tframes>0)) {

            if (tframes==weedf) {  /* first series -> save time-array */
	            for (k=0;k<weedf; k++) times[k] = time[k];
            }

            /* Calculate MSD for last weedf frames (short time average) */
            for (k=1; k<weedf; k++) {
                kk = k * np;
                /* Loop over atoms */
                for (j=0; j<np; j++) {
                    jj = kk + j;
                    xx = dx[jj] - dx[j]; xx = xx*xx;
                    yy = dy[jj] - dy[j]; yy = yy*yy;
                    zz = dz[jj] - dz[j]; zz = zz*zz;
	                aa = xx+yy+zz;
                    x2[jj] += xx;
                    y2[jj] += yy;
                    z2[jj] += zz;
                    r2[jj] += aa;
                }
            }

            for (j=0; j<np; j++) { /* copy frame 0 to save list */
                jj = frames*np+j;
                dx[jj] = dx[j];
                dy[jj] = dy[j];
                dz[jj] = dz[j];
            }

            time[frames] = time[0];  /* copy time also */
            frames++;    /* count frame */

        }

        /* Process new frame */
        if (noshiftcm == 0) {
            shiftcm(positions, tframes, natoms, dx+(ii*np), dy+(ii*np), dz+(ii*np));
        } else {
            copyfrm(positions, tframes, natoms, dx+(ii*np), dy+(ii*np), dz+(ii*np));
        }

        calccm(ns, dx+(ii*np), dy+(ii*np), dz+(ii*np), nmol,
                   dx+(ii*np+ns), dy+(ii*np+ns), dz+(ii*np+ns), g2flag);

        /* Update the frames count */
        tframes++;
        time[ii] = tframes * timestep;
    }
    printf("\n");
    end = clock();
    printf("\t\t[Step 1 of 3] Timing : %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    if (frames == meanf) {
        fprintf(stderr, "msdmol warning: reached max number of frames - increase -m nfrm ?");
    } else {
       ii = tframes%(weedf);
       if ((ii==0) && (tframes > 0)){  /* if last series just terminated */
           /* calculate MSD for last weedf frames (short time average) */
           /* debuginfo fprintf(stderr,"count short time (last)\n"); */
          for (k=1; k<weedf; k++) {  /* loop over time differences */
              kk = k*np;
              for (j=0; j<np; j++) {  /* loop over atoms */
                  jj=kk  +j;
                  xx = dx[jj] - dx[j]; xx = xx*xx;
                  yy = dy[jj] - dy[j]; yy = yy*yy;
                  zz = dz[jj] - dz[j]; zz = zz*zz;
                  aa = xx+yy+zz;
                  x2[jj] += xx;
                  y2[jj] += yy;
                  z2[jj] += zz;
                  r2[jj] += aa;
              }
          }
          normcorrect=-1;
       }
       normcorrect++; /* there may exist one more interval than sub time series */
       for (j=0; j<np; j++) {  /* copy frame 0 to save list */
           jj = frames*np+j;
           dx[jj] = dx[j];
           dy[jj] = dy[j];
           dz[jj] = dz[j];
       }
       time[frames] = time[0];  /* copy time also */
       frames++;    /* count frame */
    }

    /* *** calculate averages */

    /* average over starting times */
    start = clock();
    printf("\t\t[Step 2 of 3] \n");
    for (i=weedf; i<frames-1; i++) {  /* loop over frames */

        // Progress bar
        if (i % ((frames -1) / 1000 + 1) == 0 || i == (frames - 1) - 1) {
            printf("\r\t\t[Step 2 of 3] Progress: %5.1f%%", (i * 1000) / (frames-1)/10.0);
            fflush(stdout);
        }

        /* first weedf frames are already calculated during reading */
        ii = i*np;
        for (k=1; k<frames-i; k++) {  /* loop over time differences */
            kk = (k+weedf)*np;
            jj = (i+k)*np;
            for (j=0; j<np; j++) {  /* loop over atoms */
                xx = dx[ii+j] - dx[jj+j]; xx = xx*xx;
                yy = dy[ii+j] - dy[jj+j]; yy = yy*yy;
                zz = dz[ii+j] - dz[jj+j]; zz = zz*zz;
	            aa = xx+yy+zz;
                x2[kk+j] += xx;
                y2[kk+j] += yy;
                z2[kk+j] += zz;
                r2[kk+j] += aa;
	       }
       }
   }
    printf("\n");
    end = clock();
    printf("\t\t[Step 2 of 3] Timing : %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

   /* average over atoms  */
    start = clock();
    printf("\t\t[Step 3 of 3] \n");
    ii=1+nsit/12;     /* number of inner monomers used for averaging of r6 */
    ij = (nsit)/2;    /* center of molecule */
    for (k=1; k<frames; k++) {  /* loop over time differences */
       /* normalization factors */
        if ( k<weedf) {  /* first weedf frames have less starting points */
            aa = 1.0/ (double) (ns*(frames-weedf-normcorrect));
            xx = 1.0/ (double) (nmol*(frames-weedf-normcorrect));
            yy = 1.0/ (double) (2*ii*(frames-weedf-normcorrect));
        } else {
            aa = 1.0/ (double) (ns*(frames-k));
            xx = 1.0/ (double) (nmol*(frames-k));
            yy = 1.0/ (double) (2*ii*(frames-k));
        }

        x2av[k] = y2av[k] = z2av[k] = r2av[k] = 0.0;
        r4av[k] = 0.0;
        x4av[k] = y4av[k] = z4av[k] = 0.0;  /* xyz */
        for (j=0; j<nsit; j++) {  /* loop over atoms per molecule */
            kk=k*nsit+j;
            r2a[kk] = 0.0;
            dx[kk] = 0.0;
            dy[kk] = 0.0;
            dz[kk] = 0.0;
        }
       /* average over all monomers */
        for (j=0; j<ns; j++) {  /* loop over atoms */
            kk=k*np+j;
            x2av[k] += x2[kk];
            y2av[k] += y2[kk];
            z2av[k] += z2[kk];
            r2av[k] += r2[kk];
            jj = k*nsit+(j%nsit);
            r2a[jj] += r2[kk];   /* average of monomer position inside molecule */
            dx[jj] += x2[kk];
            dy[jj] += y2[kk];
            dz[jj] += z2[kk];
        }

        /* average of center of mass */
        for (j=ns; j<np; j++) {  /* loop over molecules */
            /* center of mass coordinates are stocked after particles in r2 */
            kk=k*np+j;
            r4av[k] += r2[kk];
            /* xyz */
            x4av[k] += x2[kk];
            y4av[k] += y2[kk];
            z4av[k] += z2[kk];
        }

        /* average of r6: use 2*ii inner monomers */
        /*    for (j=ij-ii; j < ij+ii; j++) /* average over 2*ii inner monomers */
        /*  r6av[k] += r6[k*np+j];  */
        /* now normalize everything */
        x2av[k] *= aa;
        y2av[k] *= aa;
        z2av[k] *= aa;
        r2av[k] *= aa;
        r4av[k] *= xx;  /* center of mass: normalize with number of molecules */
        x4av[k] *= xx;
        y4av[k] *= xx;
        z4av[k] *= xx;

        for (j=0; j<nsit; j++) {  /* loop over atoms per molecule */
            kk=k*nsit+j;
            r2a[kk] *= xx;
            dx[kk] *= xx;
            dy[kk] *= xx;
            dz[kk] *= xx;
        }
    }
    printf("\n");
    end = clock();
    printf("\t\t[Step 3 of 3] Timing : %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    /* restore initial times of first weedf frames */
    for (j=0;j<weedf; j++)  time[j] = times[j];

    /* check time differences (recognize time jumps between trajectories) */
    dtime = time[weedf+1]-time[0];
    j=0;
    for (k=weedf+2; k<frames; k++) {
        if (fabs((time[k] - time[k-1]) -dtime)> 1e-6) {
            fprintf(stderr,"msdmol warning: irregular time interval at virtual frame %d: %lg -> %lg.\n",k-weedf,time[k-1],time[k]);
            j +=1;  /* count number of wrong intervals */
        }
    }

    if ((timecorrect==1) && (j>0) ) {  /* correct time differences */
        fprintf(stderr,"msdmol info: correct times under assumption of equidistance dt=%lg with respect to starttime %lg\n",dtime,time[0]);
        for (k=weedf+2; k<frames; k++) {
            /* assume equidistant trajectory */
            time[k] = time[0] + (k-weedf)*dtime;
        }
    }

    /* *** output */
//    FILE *fp = fopen("msd.dat", "w");
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        perror("Error opening file");
        return 1;
    }

    if (g2flag == 1) {
        fprintf(fp, "# calculated MSD in cm frame of each chain (monomer resolved g2)\n");
    }
    if (noshiftcm == 1) {
        fprintf(fp, "# calculate in absolute frame: total CM is NOT considered\n");
    }
    fprintf(fp, "# weed=%d maxmean=%d #frames=%d/%d(+%d) #particles=%d nmol=%d\n",
            weedf, meanf, tframes, frames - weedf - normcorrect, normcorrect, ns, nmol);

    if (filmoutput == 0) {
        fprintf(fp, "# dt x2 y2 z2  g0(all) g3(cm) g3_xz  g4(end) ... g1(inner) (7+%d)\n", ij);
        for (k = 1; k < frames; k++) {
            if (k != weedf) {
                fprintf(fp, "%lf %lg %lg %lg %lg %lg %lg ",
                        time[k] - time[0],
                        x2av[k], y2av[k], z2av[k],
                        r2av[k], r4av[k], x4av[k] + z4av[k]);
                if (nsit > 2) {
                    if (nosym == 1) {
                        for (j = 0; j < nsit; j++) {
                            fprintf(fp, " %lg", r2a[k * nsit + j]);
                        }
                    } else {
                        for (j = 0; j < (nsit + 1) / 2; j++) {
                            fprintf(fp, " %lg", (r2a[k * nsit + j] + r2a[k * nsit + nsit - j - 1]) * 0.5);
                        }
                    }
                }
                fprintf(fp, "\n");
            }
        }
    } else {
        fprintf(fp, "# dt x2 y2 z2  g0(all) g3(cm) g3_xz g0_xz g3_x g3_y g3_z  g4(end) g1(inner) g1av(%d inner) g1av_x g1av_y g1av_z  ...components xyz:g_e,g_2,4,8...\n", 2 * innerav);
        for (k = 1; k < frames; k++) {
            if (k != weedf) {
                fprintf(fp, "%lf %lg %lg %lg %lg %lg %lg",
                        time[k] - time[0],
                        x2av[k], y2av[k], z2av[k],
                        r2av[k], r4av[k], x4av[k] + z4av[k]);

                fprintf(fp, " %lg %lg %lg %lg",
                        x2av[k] + z2av[k],
                        x4av[k], y4av[k], z4av[k]);

                if (nsit > 2) {
                    if (nosym == 1) {
                        for (j = 0; j < nsit; j++) {
                            fprintf(fp, " %lg", r2a[k * nsit + j]);
                        }
                    } else {
                        fprintf(fp, " %lg", (r2a[k * nsit] + r2a[k * nsit + nsit - 1]) * 0.5);
                        fprintf(fp, " %lg", (r2a[k * nsit + (nsit - 1) / 2] + r2a[k * nsit + nsit - 1 - (nsit - 1) / 2]) * 0.5);

                        aa = xx = yy = zz = 0.0;
                        for (j = (nsit + 1) / 2 - innerav; j < (nsit + 1) / 2; j++) {
                            aa += (r2a[k * nsit + j] + r2a[k * nsit + nsit - j - 1]) * 0.5;
                            xx += (dx[k * nsit + j] + dx[k * nsit + nsit - j - 1]) * 0.5;
                            yy += (dy[k * nsit + j] + dy[k * nsit + nsit - j - 1]) * 0.5;
                            zz += (dz[k * nsit + j] + dz[k * nsit + nsit - j - 1]) * 0.5;
                        }

                        fprintf(fp, " %lg %lg %lg %lg", aa / innerav, xx / innerav, yy / innerav, zz / innerav);

                        for (j = 1; j <= (nsit + 1) / 2; j *= 2) {
                            fprintf(fp, " %lg", (dx[k * nsit + j - 1] + dx[k * nsit + nsit - j]) * 0.5);
                            fprintf(fp, " %lg", (dy[k * nsit + j - 1] + dy[k * nsit + nsit - j]) * 0.5);
                            fprintf(fp, " %lg", (dz[k * nsit + j - 1] + dz[k * nsit + nsit - j]) * 0.5);
                        }
                    }
                }
                fprintf(fp, "\n");
            }
        }
    }

}

//// ==============================================================================
//void c_com_chain(int mols[], float mass[], float coords_unwrap[],
//                 int nchains, int maxatomsch, float rgsq_ich_iframe[])
//{
//
//    // Variables
//    int ich, inx;
//    int iat, iatom;
//    int dim=3;
//    int dimrg = 4;
//    float* com = malloc(nchains*dim*sizeof(float));
//    float* mass_ich = malloc(nchains*sizeof(float));
//    int* nat_ich = malloc(nchains*sizeof(float));
//
//    // Initializate center of mass (com)
//    for (ich = 0; ich < nchains; ich++) {
//        mass_ich[ich] = 0.0;
//        nat_ich[ich] = 0;
//        for (inx = 0; inx < dim; inx++) {
//            com[ich*dim+inx] = 0.0;
//        }
//    }
//
//    // DEBUG
//    // for (ich = 0; ich < nchains; ich++) {
//    //      printf("ich:%d %f %f %f\n", ich, com[ich*dim+0], com[ich*dim+1], com[ich*dim+2] );
//    //  }
//
//    // Calculate radius of gyration. The three first columns are Rgx, Rgy and Rgz.
//    // The fourth column is the Rg^2 = Rgx^2 + Rgy^2 + Rgz^2
//    #ifdef PARALLEL
//    #pragma omp parallel for private(ich, iat, iatom, inx) \
//    shared(nchains, com, rgsq_ich_iframe, mass_ich)
//    #endif
//   for (ich = 0; ich < nchains; ich++) {
//        // Calculate center of mass (com)
//        for (iat = 0; iat < maxatomsch; iat++) {
//            iatom = mols[ich*maxatomsch+iat];
//            mass_ich[ich] += mass[iatom];
//            for (inx = 0; inx < dim; inx++) {
//                com[ich*dim+inx] += mass[iatom] * coords_unwrap[iatom*dim+inx];
//            }
//        }
//        for (inx = 0; inx < dim; inx++) {
//            com[ich*dim+inx] = com[ich*dim+inx]/mass_ich[ich];
//        }




#endif