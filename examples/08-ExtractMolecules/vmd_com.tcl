proc newRep { sel type color rep imol} {
    mol selection $sel
    mol representation $type
    mol addrep $imol
    mol showrep $imol $rep on
    mol modcolor $rep $imol $color
}

set dir ""

display projection orthographic
axes location off
color Display Background white
display depthcue off


mol new ./00-PREPARE/EVOH_residues.pdb type pdb
set imol1 [molinfo top]
mol delrep 0 $imol1
set rep1 0
newRep "all" "CPK" "Name" $rep1 $imol1

mol new coords_com.pdb type pdb
set imol2 [molinfo top]
mol delrep 0 $imol2
set rep2 0
newRep "all" "vdw 0.6" "Index" $rep2 $imol2

pbc box

