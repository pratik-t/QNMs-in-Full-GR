# Quasi-Normal Oscillation Mode Solver for Neutron Stars

This Fortran code computes the $l = 2$ quasi-normal oscillation modes (QNMs) — $g_1$, $f$, and $p_1$ — of neutron stars by solving the perturbed Einstein field equations. It takes as input an equation of state (EoS) relating pressure to energy density above nuclear saturation. Multiple EoSs can be solved together using CPU parallelization. The code outputs mass, radius, dimensionless tidal deformability, and QNMs using both the Cowling approximation and a fully relativistic formalism, for a list of central energy densities and pressures. Also outputs the damping times from the full GR calculation. Information about the formalisms and results can be found in the following publications:

- [Phys. Rev. D 110, 103045 (2024)](https://doi.org/10.1103/PhysRevD.110.103045)  
- [MNRAS stae834 (2024)](https://doi.org/10.1093/mnras/stae834)  
- [Eur. Phys. J. C 10052, 13066 (2024)](https://doi.org/10.1140/epjc/s10052-024-13066-0)  
- [Phys. Rev. D 107, 103054 (2023)](https://doi.org/10.1103/PhysRevD.107.103054)
---

## Repository Structure

```text
├── BUILD/             # Compiled objects and module files
├── DATA/              # Output data and EoS folders
│   ├── EOS/           # Folder containing all input equations of state
│   ├── EOS_DATA/      # Folder containing all output data
│   └── README.md      # Notes about input and output files
├── MATH/              # LSODA ODE solver source files
├── SRC/               # Fortran source files
│   ├── Main.f90       # Main driver
│   ├── IO.f90         # I/O utilities
│   └── ...
├── EOS_inputs.txt     # List of EOS filenames to process
├── Params.txt         # Control parameters
├── Makefile           # Build script (GNU Make)
└── README.md          # This documentation
```
---

## User Inputs

### `DATA/EOS/`

- This folder should contain **all** the EoS files to be processed. Each file must contain at least two comma-separated columns which are the energy densities and corresponding pressures. The code accepts data in either $MeV/fm^3$ or CGS units. The column ordering and units can be controlled using flags in `Params.txt`.

- Quite a few equations of state are already provided, most of them taken from [CompOSE](https://compose.obspm.fr/). 

- For hybrid star EoSs with a first-order quark-hadron phase transition, a **dummy third column** is required. This column must contain a **single zero** at the beginning of the discontinuity. For example, in `DD2_HQ.csv`, _line 935_ contains a zero. Although the third column may represent number density, its values are **not used** — they can be arbitrary non-zero values, as long as the start of the discontinuity is marked with a single zero. There must **not** be any pressure or energy data between the start and end of the discontinuity. To verify that the discontinuity has been correctly detected, uncomment _lines 143–144_ in `Main.f90` and build again.

### `EOS_inputs.txt`

This file should contain the **names** of the EoS files to be solved in the current run. Each name should appear on a new line. The code processes each EoS using a separate CPU thread. Some EoS names are provided by default.

### `Params.txt`

This file specifies input and output parameters using logical flags. Each flag must be set to `.TRUE.` or `.FALSE.` :

| Flag       | Description |
|------------|-------------|
| `IS_PE`    | `.TRUE.` - if EoS columns are Pressure, Energy and `.FALSE.` - if columns are Energy, Pressure |
| `IS_MEV`   | `.TRUE.` - if units are in $MeV/fm^3$ and `.FALSE.` - if in CGS |
| `IS_STOP`  | `.TRUE.` - to stop when maximum TOV mass is reached and `.FALSE.` - to continue even with unphysical masses |
| `IS_TIDE`  | `.TRUE.` - to calculate the dimensionless tidal deformability |
| `IS_COWL`  | `.TRUE.` - to calculate QNM frequencies ($g_1$, $f$, $p_1$) using the Cowling approximation |
| `IS_GRG`   | `.TRUE.` - to calculate the fully relativistic $g_1$ mode |
| `IS_GRF`   | `.TRUE.` - to calculate the fully relativistic $f$ mode |
| `IS_GRP1`  | `.TRUE.` - to calculate the fully relativistic $p_1$ mode |

---

## Requirements

This code was developed in a Linux environment. It also works on Windows, but requires changes to the `Makefile` and the `make_filename` subroutine in `IO.f90` to ensure Windows-compatible file paths.

A Fortran compiler is needed. By default the compilation is performed using `gfortran`. For other compilers, adjust the `Makefile` accordingly. The code also requires LAPACK and BLAS libraries. On Linux, they can be install with:

```bash
sudo apt install liblapack-dev libblas-dev
```
---

## Building and Console Output

Make sure the terminal is opened in the **parent directory** of the repository. For the first run, a fresh build is needed:

```bash
make clean
make
```

This will compile and execute the code. For subsequent runs (without cleaning), you can simply use:

```bash
make
```

During execution, currently processing file and total files will be displayed in the terminal. Central energy densities (in $MeV/fm^3$) will also be displayed next to `RUNNING:`. This provides a rough indication of progress through the input EoS file.

To display additional quantities such as mass, radius, and QNM frequencies in the console output, uncomment and modify the relevant variables at _line 302_ in `Main.f90`. Note that some results may require unit conversions to match expected physical units.

---

## Notes

- The $g_1$ mode in full GR often requires fine‑tuning of the integrator to function reliably. It may work for certain central energy densities but fail quite frequently.
- Output files may occasionally contain `NaN` values. These typically result from failed integration or root‑finding steps.
- If the full GR calculations are set to `.FALSE.`, the code will run significantly faster. GR calculations are computationally intensive and sensitive to numerical settings, and may occasionally produce unphysical results. Such data points should be filtered out when analysing or plotting results.

---

## Integrator Source

The `MATH/` directory contains the thread‑safe LSODA and LSODAR source files, adapted from the ODEPACK collection.

---

## References

- Wogan, N. (2020). *ODEPACK: A collection of Fortran solvers for the solution of ordinary differential equation initial value problems* (thread‑safe LSODA/LSODAR version). GitHub repository. https://github.com/Nicholaswogan/odepack  
- **License**: ODEPACK has been declared to be in the **Public Domain**. 