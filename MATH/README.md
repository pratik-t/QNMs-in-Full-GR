This folder contains the thread‑safe LSODA and LSODAR integration routines adapted from ODEPACK. These routines are critical for both TOV and QNM integrations.

`Makefile` compiles these files as follows:

```bash
gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack_common.f90 -o odepack_common.o
gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack_interface.f90 -o odepack_interface.o
gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack_mod.f90 -o odepack_mod.o
gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack.f -o odepack.o
gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack_sub1.f -o odepack_sub1.o
gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack_sub2.f -o odepack_sub2.o
ar rcs libodepack.a odepack_mod.o odepack_interface.o odepack_common.o odepack.o odepack_sub1.o odepack_sub2.o
```

`LSODA Comments.md` contains details about the integrator setup and parameters. 