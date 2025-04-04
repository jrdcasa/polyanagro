# energy_analysis

---

## **Syntax**

---

```bash
usage: energy_analysis info [-h] -e ENERGY_LIST [ENERGY_LIST ...] [--log LOG]

optional arguments:
  -h, --help            show this help message and exit
  -e ENERGY_LIST [ENERGY_LIST ...], --energy ENERGY_LIST [ENERGY_LIST ...]
                        Energy file from MD package. The package is detected by the extension of the file
  --log LOG             Name of the file to write logs from this command
```

```bash
usage: energy_analysis calc [-h] -e ENERGY_LIST [ENERGY_LIST ...] [--log LOG] [--tbegin TBEGIN] [--tend TEND]
                            [--joinpath JOINPATH] [--groupterms] [--avg] [--acf ACF_LIST [ACF_LIST ...]]

optional arguments:
  -h, --help            show this help message and exit
  -e ENERGY_LIST [ENERGY_LIST ...], --energy ENERGY_LIST [ENERGY_LIST ...]
                        Energy file from MD package. The package is detected by the extension of the file
  --log LOG             Name of the file to write logs from this command
  --tbegin TBEGIN       Starting time in ps to perform the analysis. Example: A value of 10, start the analysis at 10ps.
  --tend TEND           Ending time in ps to perform the analysis. Example: A value of 30, end the analysis at 30ps.
  --joinpath JOINPATH   Path to the program to join energy files. Example: For GROMACS --> /usr/bin/gmx
  --groupterms          Group energy terms (bond, nonbond, ...).
  --avg                 Calculate averages from <begin> to <end>.
  --acf ACF_LIST [ACF_LIST ...]
                        Autocorrelation function (ACF) of the time series. A list with the labels of the parameter to
                        calculate the ACF

```

## **Description**

---

The **energy_analysis** command will (...)

---

---

## **Examples**

---

