01-PVL4_C0_LOPLS/03-PREPARE_LIBRARY_POLYPLY_40Mon/04-MD_NPT_0000_0200ns_500K

p "./msd_fftw3.dat" u ($1*20):($2+$3+$4) w p, "./02-YASP_MSD_ALLATOM/out_msd.dat" u 1:5, "./01-GROMACS_MSD_ALLATOM/msd_allinternal_0000-0010ns.xvg" u 1:($2*100)
