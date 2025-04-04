import os

label_natoms = "ITEM: NUMBER OF ATOMS"
label_coords = "ITEM: ATOMS"
filename="../02-YASP_MSD_ALLATOM/01-PVL4_C0_LOPLS_trj.lammpstrj"

# Open file and read lines
with open(filename, 'r') as flammps:
    lines = flammps.readlines()

# Locate the number of atoms
for idx, iline in enumerate(lines):

    if iline.find(label_natoms) != -1:
        natoms = int(lines[idx+1])
        break

# Locate and write the coordinates to output file
filename_out = os.path.splitext(os.path.split(filename)[-1])[0]+".spunto"
with open(filename_out, 'w') as fout:
    idx = 0
    iframe = 0
    fout.writelines("---pos.dat---\n")
    while idx < len(lines):
        if lines[idx].find(label_coords) != -1:
            print("FRAME {}\n".format(iframe))
            fout.writelines("#\n")
            for iat in range(idx+1, idx+natoms+1):
                i, _, x, y, z = lines[iat].split()
                line_write = "{0:s} {1:s} {2:s}\n".format(str(x), str(y), str(z))
                fout.writelines(line_write)
            idx = iat
            iline = lines[idx]
            iframe += 1
        else:
            idx += 1
    fout.writelines("----\n")


