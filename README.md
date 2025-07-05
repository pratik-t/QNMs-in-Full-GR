This code was developed for a linux environment. I have run it in windows as well, but it needs some reconfiguring of the makefile, and in the `make_filename` subroutine in `IO.f90` to ensure windows compatible file paths.

Make sure you have a fortran compiler. I use gfortran. Reconfigure the makefile for other compilers. 
lapack and lblas linear algebra packages are required, can install via: `sudo apt install liblapack-dev libblas-dev` in linux

To run, just type 'make'

Hybrid star eos with first orderquark-hadron PT should have a dummy third column where the value should be 0 at the start of discontinuity. EG: DD2_HQ has 0 in line 935. Although the third column there is number density, these values are never used, so the third column can be all ones, with a singluar 0 at the start of the discontinuity. There shouldn't be any P/E data between start and end of discont. As a check of if the discontinuity detected is as expected, uncomment lines 143-144 in Main.f90

In console, central energy densities [in Mev/fm3] are displayed in RUN to get an idea of how far along the code has run. to display M,R,P etc uncomment and add relevant variables in line 302 of Main.f90

G1-mode in full GR acts up a lot sometimes and needs some fine-tuning for the integrator. Works sometimes, but fails for a lot of central energy densities.

The output might sometimes have some NaN values. this is often due to integration/root-finding failure

IS_PE : .TRUE. if EOS columns are Pressure-Energy; .FALSE. if columns are Energy-Pressure
IS_MEV : .TRUE. if Pressure and Energy units are MeV/fm^3; .FALSE. if units are in CGS
IS_STOP : .TRUE. to stop code upon reaching maximum TOV mass; .FALSE. will also give unphysical masses
IS_TIDE : .TRUE. to calculate the dimensionless tidal deformability
IS_COWL  .TRUE. to calculate the QNM frequencies (g1,f,p1) within the Cowling approximation
IS_GRG : .TRUE. to calculate the fully relativistic g1 mode
IS_GRF : .TRUE. to calculate the fully relativistic f mode
IS_GRP1 : .TRUE. to calculate the fully relativistic p1 mode