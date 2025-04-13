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

static void c_rg_chain(int mols[], float mass[], float coords_unwrap[],
                       int nchains, int maxatomsch, float rgsq_ich_iframe[])
{

    // Variables
    int ich, inx;
    int iat, iatom;
    int dim=3;
    int dimrg = 4;
    float* com = malloc(nchains*dim*sizeof(float));
    float* mass_ich = malloc(nchains*sizeof(float));
    int* nat_ich = malloc(nchains*sizeof(float));

    // Initializate center of mass (com)
    for (ich = 0; ich < nchains; ich++) {
        mass_ich[ich] = 0.0;
        nat_ich[ich] = 0;
        for (inx = 0; inx < dim; inx++) {
            com[ich*dim+inx] = 0.0;
        }
    }

    // DEBUG
    // for (ich = 0; ich < nchains; ich++) {
    //      printf("ich:%d %f %f %f\n", ich, com[ich*dim+0], com[ich*dim+1], com[ich*dim+2] );
    //  }

    // Calculate radius of gyration. The three first columns are Rgx, Rgy and Rgz.
    // The fourth column is the Rg^2 = Rgx^2 + Rgy^2 + Rgz^2
    #ifdef PARALLEL
    #pragma omp parallel for private(ich, iat, iatom, inx) \
    shared(nchains, com, rgsq_ich_iframe, mass_ich)
    #endif
   for (ich = 0; ich < nchains; ich++) {
        // Calculate center of mass (com)
        for (iat = 0; iat < maxatomsch; iat++) {
            iatom = mols[ich*maxatomsch+iat];
            mass_ich[ich] += mass[iatom];
            for (inx = 0; inx < dim; inx++) {
                com[ich*dim+inx] += mass[iatom] * coords_unwrap[iatom*dim+inx];
            }
        }
        for (inx = 0; inx < dim; inx++) {
            com[ich*dim+inx] = com[ich*dim+inx]/mass_ich[ich];
        }

        // Calculate the Rg
        for (iat = 0; iat < maxatomsch; iat++) {
            iatom = mols[ich*maxatomsch+iat];
            nat_ich[ich] += 1;
            for (inx = 0; inx < (dimrg-1); inx++) {
               rgsq_ich_iframe[ich*dimrg+inx] = coords_unwrap[iatom*dim+inx] - com[ich*dim+inx];
            }
            rgsq_ich_iframe[ich*dimrg+(dimrg-1)] += rgsq_ich_iframe[ich*dimrg+0]*rgsq_ich_iframe[ich*dimrg+0] +
                                                    rgsq_ich_iframe[ich*dimrg+1]*rgsq_ich_iframe[ich*dimrg+1] +
                                                    rgsq_ich_iframe[ich*dimrg+2]*rgsq_ich_iframe[ich*dimrg+2];
        }

   }

   for (ich = 0; ich < nchains; ich++) {
        rgsq_ich_iframe[ich*dimrg+(dimrg-1)] = rgsq_ich_iframe[ich*dimrg+(dimrg-1)] / nat_ich[ich];
        //DEBUG
        //printf("ich: %d, nat_ich: %d %f\n",  ich, nat_ich[ich], rgsq_ich_iframe[ich*dimrg+(dimrg-1)]);
        //printf("ich: %d, nat_ich: %d %f \n", ich, nat_ich[ich], sqrt( rgsq_ich_iframe[ich*dimrg+(dimrg-1)]));
   }
   free(com);
   free(mass_ich);
}


static void c_rg_chain_massweigth(int mols[], float mass[], float coords_unwrap[],
                                  int nchains, int maxatomsch, float rgsq_ich_iframe[])
{

    // Variables
    int ich, inx;
    int iat, iatom;
    int dim=3;
    int dimrg = 4;
    float* com = malloc(nchains*dim*sizeof(float));
    float* mass_ich = malloc(nchains*sizeof(float));
    int* nat_ich = malloc(nchains*sizeof(float));

    // Initializate center of mass (com)
    for (ich = 0; ich < nchains; ich++) {
        mass_ich[ich] = 0.0;
        nat_ich[ich] = 0;
        for (inx = 0; inx < dim; inx++) {
            com[ich*dim+inx] = 0.0;
        }
    }

    // DEBUG
    // for (ich = 0; ich < nchains; ich++) {
    //      printf("ich:%d %f %f %f\n", ich, com[ich*dim+0], com[ich*dim+1], com[ich*dim+2] );
    //  }

    // Calculate radius of gyration. The three first columns are Rgx, Rgy and Rgz.
    // The fourth column is the Rg^2 = Rgx^2 + Rgy^2 + Rgz^2
   #ifdef PARALLEL
   #pragma omp parallel for private(ich, iat, iatom, inx) \
   shared(nchains, com, rgsq_ich_iframe, mass_ich)
   #endif
   for (ich = 0; ich < nchains; ich++) {
        // Calculate center of mass (com)
        for (iat = 0; iat < maxatomsch; iat++) {
            iatom = mols[ich*maxatomsch+iat];
            mass_ich[ich] += mass[iatom];
            for (inx = 0; inx < dim; inx++) {
                com[ich*dim+inx] += mass[iatom] * coords_unwrap[iatom*dim+inx];
            }
        }
        for (inx = 0; inx < dim; inx++) {
            com[ich*dim+inx] = com[ich*dim+inx]/mass_ich[ich];
        }

        // Calculate the Rg
        for (iat = 0; iat < maxatomsch; iat++) {
            iatom = mols[ich*maxatomsch+iat];
            nat_ich[ich] += 1;
            for (inx = 0; inx < (dimrg-1); inx++) {
               rgsq_ich_iframe[ich*dimrg+inx] = coords_unwrap[iatom*dim+inx] - com[ich*dim+inx];
            }
            rgsq_ich_iframe[ich*dimrg+(dimrg-1)] += (rgsq_ich_iframe[ich*dimrg+0]*rgsq_ich_iframe[ich*dimrg+0] +
                                                    rgsq_ich_iframe[ich*dimrg+1]*rgsq_ich_iframe[ich*dimrg+1] +
                                                    rgsq_ich_iframe[ich*dimrg+2]*rgsq_ich_iframe[ich*dimrg+2]) * mass[iatom];


        }
   }

   for (ich = 0; ich < nchains; ich++) {
        rgsq_ich_iframe[ich*dimrg+(dimrg-1)] = rgsq_ich_iframe[ich*dimrg+(dimrg-1)] / mass_ich[ich];
        //DEBUG
        //printf("ich: %d, nat_ich: %d %f\n",  ich, nat_ich[ich], rgsq_ich_iframe[ich*dimrg+(dimrg-1)]);
        //printf("ich: %d, nat_ich: %d %f \n", ich, nat_ich[ich], sqrt( rgsq_ich_iframe[ich*dimrg+(dimrg-1)]));
   }

   free(com);
   free(mass_ich);
}

#endif