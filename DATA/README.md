## `DATA/EOS/`

- This folder should contain **all** the EoS files to be processed. Each file must contain at least two comma-separated columns which are the energy densities and corresponding pressures. The code accepts data in either $MeV/fm^3$ or CGS units. The column ordering and units can be controlled using flags in `Params.txt`.

- Quite a few equations of state are already provided, most of them taken from [CompOSE](https://compose.obspm.fr/). 

- For hybrid star EoSs with a first-order quark-hadron phase transition, a **dummy third column** is required. This column must contain a **single zero** at the beginning of the discontinuity. For example, in `DD2_HQ.csv`, _line 935_ contains a zero. Although the third column may represent number density, its values are **not used** — they can be arbitrary non-zero values, as long as the start of the discontinuity is marked with a single zero. There must **not** be any pressure or energy data between the start and end of the discontinuity. To verify that the discontinuity has been correctly detected, uncomment _lines 143–144_ in `Main.f90` and build again.
---

## `DATA/EOS_DATA/`

- This folder contains output CSVs of processed input EOSs. The outputs can be controlled using flags in `Params.txt`. 

- Output units are mentioned in the column headers.