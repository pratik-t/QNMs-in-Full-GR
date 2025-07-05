Make sure to have odepack_common.mod, odepack_interface.mod and odepack_mod.mod to the BUILD directory of your code
Go to MATH folder and compile as follows:

gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack_common.f90 -o odepack_common.o
gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack_interface.f90 -o odepack_interface.o
gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack_mod.f90 -o odepack_mod.o
gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack.f -o odepack.o
gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack_sub1.f -o odepack_sub1.o
gfortran -march=native -flto -O3 -floop-nest-optimize -c odepack_sub2.f -o odepack_sub2.o
ar rcs libodepack.a odepack_mod.o odepack_interface.o odepack_common.o odepack.o odepack_sub1.o odepack_sub2.o

Copy the files libodepack.a and odepack_mod.o to the MATH directory of your code