## Molpro on Ubuntu with gcc/gfortran
gcc, mpich, make, cmake, git, etc. are required.

**Install Global Arrays:**  
./configure --prefix=/home/trushin/libs/ga-5.8.2  
make  
make check  
make install

**Install BLAS, lapack, lapacke:**  
sudo apt-get install liblapack-dev liblapack-doc liblapack-pic liblapack3 liblapack-test liblapacke liblapacke-dev

**Install eigen3**:  
sudo apt-get install libeigen3-dev

**Install Molpro**  
PATH=$PATH:/home/trushin/libs/ga-5.8.2/bin ./configure FOPT=-O2 --disable-xml2  
make -j 16  
make quicktest

Copy a Molpro token to /home/trushin/.molpro/ before making quicktest
## Molpro on Xubuntu with gcc/gfortran
With Xubuntu eigen was installed/copied manually and configuration was different:  
PATH=$PATH:/home/trushin/libs/ga-5.8.2/bin ./configure FOPT=-O2 CPPFLAGS=-I/home/trushin/libs/eigen-3.4.0/ --disable-xml2