# ===============================================================
proc newRep { sel type color rep imol} {
  
    mol selection $sel
    mol representation $type
    mol addrep $imol
    mol showrep $imol $rep on
    mol modcolor $rep $imol $color

}
# ===============================================================
proc assign_mass_UA_model {} {
  
    set allatoms [atomselect top "all"]
    foreach index [$allatoms list] name [$allatoms get name] {
		if {$name == "C3"} {
			set m 15.03
		} elseif {$name == "C2"}  {
			set m 14.03
		} elseif {$name == "C1"}  {
			set m 13.02
		}
		set iatom [atomselect top "index $index"]
		$iatom set mass $m
    }
}

# ===============================================================
proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        #Set mass for the resid 1
        foreach index [$selection list] coord [$selection get {x y z} ] {
          if {$index !=0 && $index != 119} {
            set m 14.02658
          } else {
            set m 15.03452
          }
          set mass [ expr {$mass + $m} ]
          set com  [ vecadd $com [vecscale $m $coord] ]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can’t be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}

# ===============================================================
proc gyr_radius {sel} {
  # make sure this is a proper selection and has atoms
  if {[$sel num] <= 0} {
    error "gyr_radius: must have at least one atom in selection"
  }
  # gyration is sqrt( sum((r(i) - r(center_of_mass))^2) / N)
  set com [center_of_mass $sel]
  set sum 0
  foreach coord [$sel get {x y z}] {
    set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
  }

  return [expr sqrt($sum / ([$sel num] + 0.0))]
}

# =====================  MAIN  ==================================
display projection orthographic
axes location off
color Display Background white
display depthcue off

mol new ../RUN-001/confout.gro
mol addfile ../RUN-001/traj_comp.xtc waitfor -1
set imol1 [molinfo top]
animate goto start

pbc box
mol delrep 0 $imol1
pbc join fragment -bondlist -all
#pbc box
set rep 0
newRep "index 0 to 23" "CPK" "ColorId 0" $rep $imol1
set rep [expr $rep + 1]
newRep "index 24 to 47" "CPK" "ColorId 1" $rep $imol1
set rep [expr $rep + 1]
newRep "index 48 to 71" "CPK" "ColorId 4" $rep $imol1

assign_mass_UA_model

set nframes [molinfo top get numframes]
for {set ifr 1} {$ifr < $nframes} {incr ifr} {

   animate goto $ifr
   set ch1 [atomselect top "index 0 to 23"]
   set ch2 [atomselect top "index 24 to 47"]
   set ch3 [atomselect top "index 48 to 71"]

   set com1 [center_of_mass $ch1]
   set rg1 [gyr_radius $ch1]
   puts "Chain 0 com: $com1"
   puts "Chain 0 Rg2: $rg1"

   set com2 [center_of_mass $ch2]
   set rg2 [gyr_radius $ch2]
   puts "Chain 1 com: $com2"
   puts "Chain 1 Rg2: $rg2"

   set com3 [center_of_mass $ch3]
   set rg3 [gyr_radius $ch3]
   puts "Chain 2 com: $com3"
   puts "Chain 2 Rg2: $rg3"
   puts "===================="

}

