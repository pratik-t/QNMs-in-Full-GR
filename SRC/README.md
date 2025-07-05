This folder contains all Fortran source files used to build the QNM solver. Key files include:

- **Main.f90**  
  - Entry point of the program.  
  - Parses `EOS_inputs.txt` and `Params.txt`, initializes threads, and orchestrates the TOV integration and QNM routines.  
  - Line 302 can be modified to print additional diagnostics (mass, radius, QNM frequencies). Units might not be as expected at this point.
  - Lines 143–144 can be uncommented to check if the EoS discontinuity detected is correct.

- **IO.f90**  
  - Handles all file input/output.  
  - Contains `make_filename` for platform‑independent path generation.  
  - Utilities for reading comma‑separated EoS tables and writing output columns.
  - Utilities used by other source files such as interpolation, matrix multiplication and curve fitting.

- **TOV.f90**  
  - Implements the Tolman–Oppenheimer–Volkoff (TOV) equations integration.

- **Cowling.f90**  
  - Solves the perturbed fluid equations with the Cowling approximation to obtain $g_1$, $f$, and $p_1$ modes.

- **FullGR.f90**  
  - Fully relativistic QNM solver.
  - Comnputes QNMs using Brent's method to find minimas in the frequency spectrum.