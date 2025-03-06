
set PSFFILE confout_unwrap_labeled.psf
set PDBTRJ  01-PVL4_C0_LOPLS_trj.pdb
set LAMMPSTRJ 01-PVL4_C0_LOPLS_trj.lammpstrj

#LOAD MOLECULES
mol new $PSFFILE type psf
animate read pdb $PDBTRJ waitfor all

animate write lammpstrj $LAMMPSTRJ 
