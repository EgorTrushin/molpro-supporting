## The contents of the directories in the repository
- **docs** Markdown documentation for our programs/codes
- **examples** example input files for calculations using our programs/codes
- **files** various auxiliary files
- **pyscripts** python scripts useful for working with Molpro and Cluster
- **testjobs** testjobs for our programs/codes

## Molpro on various Linux systems

### Molpro on OpenSUSE with Intel oneAPI (tcsv020)
**Install Global Arrays:**
```
I_MPI_CXX=icpc MPICC=mpiicc MPIF77=mpiifort ./configure --prefix=/home/trushin/libs/ga-5.8.2
make
make check
make install
```

**Install eigen3**  
eigen3 does not need to be compiled, but needs to be downloaded and unpacked into a directory.

**Install Molpro**
```
I_MPI_CXX=icpc CC=icc FC=ifort FOPT=-O2 CPPFLAGS=-I/home/trushin/libs/eigen-3.4.0/include/eigen3 PATH=$PATH:/home/trushin/libs/ga-5.8.2/bin ./configure --prefix=/home/trushin/Molpro/molpro-ksinv --bindir=/home/trushin/Molpro/molpro-ksinv --disable-gfortran-check
make -j 28 symtrans_FLAGS=-O0
make quicktest
```
The following two files should be replaced to compile (see directory 'files' in repo):  
src/util/molpro_main.cpp  
src/mrci/kext.F90

It is necessary to have /home/Tools/progs/intel/oneapi/mkl/2021.1.1/lib/intel64 in LD_LIBRARY_PATH to use compiled molpro.exe. Add this, e.g., to .bashrc:
```
export LD_LIBRARY_PATH=/home/Tools/progs/intel/oneapi/mkl/2021.1.1/lib/intel64:$LD_LIBRARY_PATH
```

### Molpro on Ubuntu with gcc/gfortran
gcc, mpich, make, cmake, git, etc. are required.

**Install Global Arrays:**
```
./configure --prefix=/home/trushin/libs/ga-5.8.2  
make
make check
make install
```
**Install BLAS, lapack, lapacke:**
```
sudo apt-get install liblapack-dev liblapack-doc liblapack-pic liblapack3 liblapack-test liblapacke liblapacke-dev
```
**Install eigen3**:
```
sudo apt-get install libeigen3-dev
```
**Install Molpro**
```
PATH=$PATH:/home/trushin/libs/ga-5.8.2/bin ./configure FOPT=-O2
make -j 16
make quicktest
```
Copy a Molpro token to /home/trushin/.molpro/ before making quicktest.

### Molpro on Xubuntu with gcc/gfortran
With Xubuntu everything worked as with Ubuntu, except that eigen was manually installed/copied and the configuration was as follows:
```
PATH=$PATH:/home/trushin/libs/ga-5.8.2/bin ./configure FOPT=-O2 CPPFLAGS=-I/home/trushin/libs/eigen-3.4.0/ --disable-xml2
```

## Automatic formatting for Fortran code
[fprettify](https://github.com/pseewald/fprettify) is a great tool. My choice for formatting:
```
fprettify -i 2 -l 80 -w 1 -s
```

## Molpro with SLURM running
One can replace LAUNCHER in bin/molpro (bin/molpro.sh)
```
LAUNCHER="/usr/bin/srun %x”
```
This helped to run properly hybrid MPI+OpenMP calculations on multiple nodes using cluster machines.
