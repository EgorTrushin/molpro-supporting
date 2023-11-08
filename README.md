**Table of content:**
- [The contents of the directories in the repository](#item_dirs)
- [Compiling Molpro on various Linux systems](#item_compile)
- [Automatic formatting for Fortran code with fprettify](#item_fprettify)
- [Compile Molpro with warnings](#item_warnings)
- [Molpro with SLURM launcher](#item_slurm)

<a id="item_dirs"></a>
## The contents of the directories in the repository
- **docs** Markdown documentation for our programs/codes
- **examples** example input files for calculations using our programs/codes
- **files** various auxiliary files
- **pyscripts** python scripts useful for working with Molpro and Cluster
- **testjobs** testjobs for our programs/codes

<a id="item_compile"></a>
## Compiling Molpro on various Linux systems

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
clone ThC-AG molpro version from GitHub:

```
git clone https://github.com/ThCErlangen/molpro
```

Navigate to the created directory and configure the installation by:

```
I_MPI_CXX=icpc CC=icc FC=ifort FOPT=-O2 CPPFLAGS=-I/home/trushin/libs/eigen-3.4.0/include/eigen3 PATH=$PATH:/home/trushin/libs/ga-5.8.2/bin ./configure --prefix=/home/trushin/Molpro/molpro --bindir=/home/trushin/Molpro/molpro --disable-gfortran-check
```

now replace the files:

src/util/molpro_main.cpp  
src/mrci/kext.F90

by files of the same name from [*files*](./files) directory
then run, from:

```
make -j 28 symtrans_FLAGS=-O0
make quicktest
```

It is necessary to have /home/Tools/progs/intel/oneapi/mkl/2021.1.1/lib/intel64 in LD_LIBRARY_PATH to use compiled molpro.exe. Add this, e.g., to .bashrc or your slurm file:

```
export LD_LIBRARY_PATH=/home/Tools/progs/intel/oneapi/mkl/2021.1.1/lib/intel64:$LD_LIBRARY_PATH
```

### Molpro on Ubuntu with gcc/gfortran
gcc, mpich, make, cmake, git, etc. are required.

**Install Global Arrays:**
```
./configure  
make
make check
make install
```
**Install BLAS, lapack, lapacke, einen3:**
```
sudo apt-get install liblapack-dev liblapack-doc liblapack-pic liblapack3 liblapack-test liblapacke liblapacke-dev libeigen3-dev
```
**Install Molpro**
```
./configure FOPT=-O2
make -j 16
make quicktest
```
Copy a Molpro token to /home/trushin/.molpro/ before making quicktest.

<a id="item_fprettify"></a>
## Automatic formatting for Fortran code with fprettify
[fprettify](https://github.com/pseewald/fprettify) is a great tool. My choice for formatting:
```
fprettify -i 2 -l 80 -w 1 -s
```

<a id="item_warnings"></a>
## Compile Molpro with warnings
### ifort
For Intel Fortran, if you want, e.g., to see warnings for unused variables configure with FCFLAGS="-warn unused". For instance
```
I_MPI_CXX=icpc CC=icc FC=ifort FOPT=-O2 CPPFLAGS=-I/home/trushin/libs/eigen-3.4.0/include/eigen3 FCFLAGS="-warn unused" PATH=$PATH:/home/trushin/libs/ga-5.8.2/bin ./configure --prefix=/home/trushin/Molpro/molpro --bindir=/home/trushin/Molpro/molpro --disable-gfortran-check
```
For other possible flags see [corresponding documentation](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-0/warn.html)

### gfortran
For gfotran just to configure with FCFLAGS=-Wall:
```
./configure FCFLAGS=-Wall
```

<a id="item_slurm"></a>
## Molpro with SLURM launcher
One can replace LAUNCHER in bin/molpro (bin/molpro.sh)
```
LAUNCHER="/usr/bin/srun %x”
```
This helped to run properly hybrid MPI+OpenMP calculations on multiple nodes using cluster machines.
