**Table of content:**
- [The contents of the directories in the repository](#item_dirs)
- [Markdown documentation for programs/codes](#item_docs)
- [Compiling Molpro on various Linux systems](#item_compile)
- [Automatic formatting for Fortran code with fprettify](#item_fprettify)
- [Compile Molpro with warnings](#item_warnings)
- [Molpro with SLURM launcher](#item_slurm)

<a id="item_dirs"></a>
## The contents of the directories in the repository
- **docs** Markdown documentation for programs/codes
- **examples** example input files for calculations using our programs/codes
- **files** various auxiliary files
- **testjobs** testjobs for our programs/codes

<a id="item_docs"></a>
## Markdown documentation for programs/codes
- [RIRPA/URIRPA program](https://www.molpro.net/manual/doku.php?id=kohn-sham_random-phase_approximation#rirpa_program)
- [SCEXX/USCEXX program](/docs/scexx_uscexx.md)
- [KSINV/UKSINV program](/docs/ksinv_uksinv.md)
- [SCRPA/USCRPA program](/docs/scrpa_uscrpa.md)

<a id="item_compile"></a>
## Compiling Molpro on various Linux systems

### Molpro on OpenSUSE with Intel oneAPI (tcsv020)

<details><summary><b>Install Global Arrays</b></summary>
  
For **ifx**:
```
I_MPI_CXX=icpx MPICC=mpiicx MPIF77=mpiifx ./configure --prefix=/home/trushin/libs/ga-5.8.2_ifx
make
make check
make install
```
For **ifort**:
```
I_MPI_CXX=icpc MPICC=mpiicc MPIF77=mpiifort ./configure --prefix=/home/trushin/libs/ga-5.8.2
make
make check
make install
```
</details>

<details><summary><b>Install eigen3</b></summary>

eigen3 does not need to be compiled, but needs to be downloaded and unpacked into a directory.
</details>

<details><summary><b>Install Molpro</b></summary>

Clone Molpro from GitHub. e.g.:

```
git clone https://github.com/molpro/molpro.git
```

Navigate to the created directory and configure the installation by:

```
FC=ifx CC=icx CXX=mpiicpx COPT=-O3 FOPT=-O1 CPPFLAGS=-I/home/trushin/libs/eigen-3.4.0/include/eigen3 LDFLAGS=-lstdc++fs PATH=$PATH:/home/trushin/libs/ga-5.8.2_ifx/bin ./configure --disable-gfortran-check --disable-aims --disable-slater --without-hdf5
```

now replace the files:

src/util/molpro_main.cpp  
src/util/remove_all.cpp  

by the file of the same name from [*files*](./files) directory
then run, from:

```
make -j 28 symtrans_FLAGS=-O0
make quicktest
```

It is necessary to have /home/Tools/progs/intel/oneapi/mkl/2021.1.1/lib/intel64 in LD_LIBRARY_PATH to use compiled molpro.exe. Add this, e.g., to .bashrc or your slurm file:

```
export LD_LIBRARY_PATH=/home/Tools/progs/intel/oneapi/mkl/2021.1.1/lib/intel64:$LD_LIBRARY_PATH
```
</details>

### Molpro on Ubuntu with gcc/gfortran
gcc, mpich, make, cmake, git, etc. are required.

<details><summary><b>Install Global Arrays</b></summary>

```
./configure  
make
make check
make install
```
</details>

<details><summary><b>Install BLAS, lapack, lapacke, einen3</b></summary>

```
sudo apt-get install liblapack-dev liblapack-doc liblapack-pic liblapack3 liblapack-test liblapacke liblapacke-dev libeigen3-dev
```
</details>

<details><summary><b>Install Molpro</b></summary>

```
./configure FOPT=-O2
make -j 16
make quicktest
```
Copy a Molpro token to /home/trushin/.molpro/ before making quicktest.
</details>

<a id="item_fprettify"></a>
## Automatic formatting for Fortran code

<details><summary><b>fprettify</b></summary>

[fprettify](https://github.com/pseewald/fprettify) is a great tool. My choice for formatting:
```
fprettify -i 2 -l 80 -w 1 -s
```
</details>

<a id="item_warnings"></a>
## Compile Molpro with warnings

<details><summary><b>ifort</b></summary>

For Intel Fortran, if you want, e.g., to see warnings for unused variables configure with FCFLAGS="-warn unused". For instance
```
I_MPI_CXX=icpc CC=icc FC=ifort FOPT=-O2 FCFLAGS="-warn unused" CPPFLAGS=-I/home/trushin/libs/eigen-3.4.0/include/eigen3 PATH=$PATH:/home/trushin/libs/ga-5.8.2/bin  ./configure --disable-gfortran-check
```
For other possible flags see [corresponding documentation](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-0/warn.html)
</details>

<details><summary><b>gfortran</b></summary>

For gfotran just to configure with FCFLAGS=-Wall:
```
./configure FCFLAGS=-Wall
```
</details>

<a id="item_slurm"></a>
## Molpro with SLURM launcher

<details><summary><b>SLURM launcher</b></summary>

One can replace LAUNCHER in bin/molpro (bin/molpro.sh)
```
LAUNCHER="/usr/bin/srun %x”
```
This helped to run properly hybrid MPI+OpenMP calculations on multiple nodes using cluster machines.
</details>
