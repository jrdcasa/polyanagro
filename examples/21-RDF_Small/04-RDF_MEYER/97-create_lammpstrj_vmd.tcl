
set PSFFILE PVLMon_residues.psf
set PDBTRJ  01-PVMon.pdb
set LAMMPSTRJ 01-PVMon.lammpstrj

#LOAD MOLECULES
mol new $PSFFILE type psf
animate read pdb $PDBTRJ waitfor all

animate write lammpstrj $LAMMPSTRJ 
