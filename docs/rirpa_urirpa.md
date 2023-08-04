## The RIRPA and URIRPA programs
The RIRPA and URIRPA programs allow non-self-consistent spin-restricted and spin-unrestricted resolution of identity (RI) random phase approximation (RPA) [1-3] and σ-functional [4-6] calculations. These methods should be used in conjunction with conventional Kohn-Sham (KS) density functional theory (DFT) calculations, i.e. data from a preceding KS DFT calculation should be provided. Conventional KS DFT means calculations with LDA, GGA and hybrid exchange-correlation functionals.

**Bibilography:**  
**RPA:**  
[1] F. Furche, [Phys. Rev. B 64](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.64.195120), 195120 (2001)  
[2] A. Heßelmann and A. Görling, [Mol. Phys.](https://www.tandfonline.com/doi/abs/10.1080/00268976.2011.614282) 109, 2473 (2011).  
[3] X. Ren, P. Rinke, C. Joas, and M. Scheffler, [J. Mater. Sci.](https://link.springer.com/article/10.1007/s10853-012-6570-4) 47, 7447 (2012).  
**σ-functionals:**  
[4] E. Trushin, A. Thierbach, A. Görling, [J. Chem. Phys.](https://aip.scitation.org/doi/full/10.1063/5.0026849) 154, 014104 (2021)  
[5] S. Fauser, E. Trushin, C. Neiss, A. Görling, [J. Chem. Phys.](https://aip.scitation.org/doi/abs/10.1063/5.0059641) 155, 134111 (2021)  
[6] C. Neiss, S. Fauser, A. Görling, [J. Chem. Phys.](https://aip.scitation.org/doi/full/10.1063/5.0129524) 158, 044107 (2023)  
**Other papers cited in the documentation:**  
[7] E. Trushin, A. Görling, [J. Chem. Phys.](https://aip.scitation.org/doi/full/10.1063/5.0056431) 155, 054109 (2021)

We kindly ask you to cite the original publications of the corresponding methods in the publications that result from these programs.

Example input file for spin-restricted calculations for the CO molecule:
```
gthresh,twoint=1d-20,energy=1d-10,orbital=1d-8 ! tighter thresholds are recommended
gdirect ! integral-direct mode

basis={
default,aug-cc-pwCVQZ ! orbital basis
set,ri; default,aug-cc-pwCVQZ/mp2fit ! RI basis
}

symmetry,nosym ! RIRPA does not use symmetry

angstrom
geometry={
2

C 0.000000 0.000000 -0.646514
O 0.000000 0.000000 0.484886
}

df-ks,pbe,df_basis=aug-cc-pwCV5Z/mp2fit,maxit=200 ! DFT calculation
{cfit,basis_coul=aug-cc-pwCV5Z/mp2fit,basis_exch=aug-cc-pwCV5Z/mp2fit}

acfd;rirpa ! RPA/σ-functional calculation
```
As well as an example of spin-unrestricted calculation for the NH<sub>2</sub> molecule:
```
gthresh,twoint=1d-20,energy=1d-10,orbital=1d-8 ! tighter thresholds are recommended
gdirect ! integral-direct mode

basis={
default,aug-cc-pwCVQZ ! orbital basis
set,ri; default,aug-cc-pwCVQZ/mp2fit ! RI basis
}

symmetry,nosym ! URIRPA does not use symmetry

angstrom
geometry={
3

N 0.000000 0.000000 0.142235
H 0.000000 0.800646 -0.497821
H 0.000000 -0.800646 -0.497821
}

spin=1

df-uks,pbe,df_basis=aug-cc-pwCV5Z/mp2fit,maxit=200 ! DFT calculation
{cfit,basis_coul=aug-cc-pwCV5Z/mp2fit,basis_exch=aug-cc-pwCV5Z/mp2fit}
 
acfd;urirpa ! RPA/σ-functional calculation
```
The following options are available for the RIRPA and URIRPA programs:
* **orb** record number containing the orbital coefficients, eigenvalues, etc. from the preceding DFT calculation (default: ‘2100.2’ and ‘2200.2’ for RIRPA and URIRPA, respectively)
* **dfit** logical flag to enable density fitting during the reference energy calculation (default: ’1’)
* **sigma** logical flag to enable σ-functional calculation (default: ’1’)
* **sigma_param** string containing a name for the parametrization used. Choose 'PBE_S2' [6], 'PBE0_S2' [6], 'TPSS_W' 1[5], 'B3LYP_W1' [5] parameterisation in combination with a preceding DFT calculation with PBE, PBE0, TPSS or B3LYP exchange correlation functional, respectively (default: ‘PBE_S2’)
* **write_sigma** logical flag to enable writing of sigma.dat file with reference energy, frequency integration weights and σ-values (default: ’0’)
* **thr_overlap_ri** threshold for processing RI basis according to Section IIB2 in Ref. [7] (default: ‘1d-99’)
* **thr_fai_ri** threshold for processing RI basis according to Section IIB5 in Ref. [7] (default: ‘1d-14’)
* **thr_rpa** threshold for throwing out contributions corresponding to small eigenvalue differences during construction of the response function (default: ‘1d-6’)
* **nquadint** number of logarithmically spaced intervals for frequency integration (default ‘1’)
* **nquad** number of points per interval for frequency integration (default '50')
* **w0** caling factor for rational the function mapping the Gauss–Legendre quadrature for the interval [−1, 1] to the interval [0, ∞], see Eqs. 37-38 in Ref. [4] for details (default: ‘2.5’)
* **vc_scal** scaling factor for the Coulomb kernel, which can be used to mimic the effect of the inclusion of the exact-exchange kernel. In the special case of non-spin-polarized two-electron systems, the RPA calculation with a Coulomb kernel scaled by 1/2 is equivalent to including of the exact-exchange kernel. Implemented only in RIRPA (default: ‘1d0’)
* **verb** determines the level of verbosity in the output file, integer values of 0, 1, 3 provide different levels of verbosity (default ’0’)