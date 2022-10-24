# info_trj

---

## **Syntax**

---

```bash
usage: info_trj [-h] -t TRAJ [TRAJ ...] [--log LOG] (--topo TPR|DATA|PDB)

Get information from a MD trajectory. This is part of the polyanagro library

optional arguments:
  -h, --help            show this help message and exit
  -t TRAJ [TRAJ ...], --traj TRAJ [TRAJ ...]
                        A list of trajectories from MD simulations. Allowed trajectories are XTC and DCD.
  --log LOG             Name of the file to write logs from this command
  --topo TPR|DATA|PDB   A topology file in tpr, data or pdb format. tpr --> GROMACS, dat --> LAMMPS, pdb --> OTHERS

```


---

## **Description**

---

The **info_trj** command will print some useful information about the trajectory and topology. The command produce a log file (**info_trj.log**)

---

## **Examples**

---

* **Analysis from Gromacs output.**

A (or a list of) xtc trajectory and a tpr topology file from GROMACS are required.

   ```bash
  cd examples/01-Info_trj_gromacs

  info_trj -t ../../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc ../../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc --topo ../../data/0003Ch-C020-002br04/RUN-001/topol.tpr

  ```

* **Analysis from Lammps output.**

A dcd trajectory and a data topology file from LAMMPS are required.

  ```bash
    cd examples/02-Info_trj_lammps

    info_trj -t  ../../data/PE_OPLS_40Ch_102C_LAMMPS_1ns/dump.dcd --topo ../../data/PE_OPLS_40Ch_102C_LAMMPS_1ns/restart.data.500000.data
  ```

* **Analysis from a PDB file used as topology**

A dcd or xtc trajectory and a data topology file from LAMMPS are required.

  ```bash
    cd examples/03-Info_trj_pdb

    info_trj -t  ../../data/PE_OPLS_40Ch_102C_LAMMPS_1ns/dump.dcd --topo ../../data/PE_OPLS_40Ch_102C_LAMMPS_1ns/PE102_40Ch_residues_replicate.pdb

  ```
