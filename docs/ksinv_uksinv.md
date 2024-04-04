## KSINV and UKSINV programs
The programs `KSINV` and `UKSINV` solve the inverse Kohn-Sham (KS) problem in the spin-restricted and spin-unrestricted cases, respectively.
In the spin-restricted case, the KS inversion determines the KS potential for a given target density.
In the spin-unrestricted case, the KS potentials for the $\alpha$ and $\beta$ spin are determined for given target electron densities for $\alpha$ and $\beta$ electrons.

The program is designed to work with densities provided by many-body methods such as Coupled Cluster (CC) or Configuration Interaction (CI). The implementation for closed-shell systems was described and tested in Ref. [1] and follows the idea of Ref. [2] to use the KS response function to update the KS potential. The extension to open-shell systems was described in Ref. [X].

---

**Important:** Read carefully the information below about the selection of basis sets and the `thr_fai_oep` parameter:  
Kohn-Sham inversion calculations requires to specify two basis sets, namely orbital and OEP basis sets. The orbital basis set is the primary basis set for e.g. the representation of orbitals and the calculation of energy contributions. The OEP basis set is the auxiliary basis set for representing e.g. the exchange-correlation potential. 
Gaussian basis set KS inversion methods for a long time suffered from numerical instabilities resulting from problems in balancing the orbital basis sets with the auxiliary basis sets.
Present implementation of KS inversion in Molpro use a preprocessing of the auxiliary basis set that enables such balancing in an automatic fashion for standard Gaussian basis sets.
The preprocessing step removes linear combinations of auxiliary basis functions that couple poorly to products of occupied times unoccupied Kohn-Sham orbitals. It involves the threshold `thr_fai_oep`, which determines how many linear combinations of auxiliary basis functions are removed and varies with respect to the size of the orbital basis set used. This scheme was tested using Dunning correlation consistent basis sets aug-cc-pwCVXZ (X = T, Q, 5) and recommended thresholds are
| Orbital basis set | **thr_fai_oep** |
| :----: | :----: |
| aug-cc-pwCVTZ  |  5e-2   |
| aug-cc-pwCVQZ  | 1.7e-2  |
| aug-cc-pwCV5Z  |  5e-3   |

These thresholds are also expected to work for orbital basis sets without augmentation cc-pwCVXZ (X = T, Q, 5) and without core-polarization functions aug-cc-pVXZ (X = T, Q, 5).
 
As an auxiliary basis set (OEP), the aug-cc-pVXZ/mp2fit (X=T,Q,5) family of basis sets is recommended. In particular, it is recommended to use the aug-cc-pVDZ/mp2fit auxiliary basis sets for atoms up to neon and the aug-cc-pVTZ/mp2fit auxiliary basis sets for heavier atoms.
 
Def2 basis sets of Ahlrichs family of similar quality also might work, but were not tested extensively in practice.

---

Example input file for KS inversion of the reference density provided by the CCSD method for the CO molecule:
```
basis={
default,aug-cc-pwCVQZ ! orbital basis
set,oep;default,aug-cc-pVDZ/mp2fit ! OEP basis
}

symmetry,nosym ! KSINV does not use symmetry

angstrom
geometry={
2

C        0.000000    0.000000   -0.646514
O        0.000000    0.000000    0.484886
}

hf,maxit=30 ! HF calculations

{CCSD;core;dm,1325.1} ! reference CCSD calculation
e_ref_ccsd=energy ! save energy to provide as a reference to KSINV

! save relaxed density to the first record - required by KSINV
{matrop;load,den,den,1325.1;save,den,1325.1}

acfd;ksinv,refden=1325.1,e_ref=e_ref_ccsd,thr_fai_oep=1.7d-2 ! KSINV calculation
```
It may be useful to split the reference density and KS inversion calculations into two separate calculations. In this case, the reference density can be saved as:
```
basis={
default,aug-cc-pwCVQZ ! orbital basis
}

symmetry,nosym ! KSINV does not use symmetry

angstrom
geometry={
2

C        0.000000    0.000000   -0.646514
O        0.000000    0.000000    0.484886
}

hf,maxit=30 ! HF calculations

{CCSD;core;dm,1325.1} ! reference CCSD calculation
{matrop;load,den,den,1325.1;write,den,dm.dat,status=rewind,prec=sci}
```
And then the KS inversion can be performed:
```
basis={
default,aug-cc-pwCVQZ ! orbital basis
set,oep;default,aug-cc-pVDZ/mp2fit ! OEP basis
}

symmetry,nosym ! KSINV does not use symmetry

angstrom
geometry={
2

C        0.000000    0.000000   -0.646514
O        0.000000    0.000000    0.484886
}

hf,maxit=0 ! HF calculation with 0 iterations, KSINV uses this for initialization

{matrop;read,den,type=density,file=dm.dat;save,den,1325.1} ! read stored density matrix

acfd;ksinv,refden=1325.1,e_ref=-113.285493180105,thr_fai_oep=1.7d-2 ! KSINV calculation
```
The following options are available for the KSINV and UKSINV programs:
- **refden** record from which to read the reference density
- **orb** record from which the occupation numbers are read (default: ‘2100.2’)
- **save** record in which the resulting orbital coefficients, eigenvalues, etc. are written (default: '2101.2')
- **e_ref** total energy of calculations from which reference density is provided
- **thr_inversion** convergence threshold for self-consistent procedure (default '1d-10')
- **maxit** maximum number of iterations (default '30')
- **minit** minimum number of iterations (default '3')
- **mixing_a**  real parameter to control the update of the exchange correlation potential in the KS inversion procedure. The parameter corresponds to $a$ in Appendix A of Ref. [1] (default: '0.1')
- **mix_switch_iter** integer parameter to control the update of the exchange correlation potential in the KS inversion procedure. The parameter corresponds to $i$ in Appendix A of Ref. [1] (default: '3')
- **backfilter** if set to $\neq$ 0, Eq. (52) of Ref. [1] is applied at each iteration to remove contributions not contained in the space spanned by the current iteration's auxiliary basis set (default: '1')
- **thr_sym** threshold for symmetrization of the OEP basis set to enforce that OEP basis exhibits full symmetry of molecule. Set the threshold to 1d-10 to enable symmetrization. (default: ‘0d0’)
- **thr_fai_oep** threshold for processing OEP basis according to Section IIB5 in Ref. [3] (default: ‘1.7d-2’)
- **thr_overlap_oep** threshold for processing OEP basis according to Section IIB2 in Ref. [3] (default: ‘1d-99’)
- **vref_fa** if set to $\neq$ 0, enable the use of the Fermi-Amaldi potential as reference potential when applying the charge condition. Otherwise, the reference potential is constructed according to Eq. (45) of Ref. [3] (default: '1')
- **thr_oep** threshold for throwing out contributions corresponding to small eigenvalue differences appearing in the denominator when constructing the so-called lambda term $1/(\varepsilon_a - \varepsilon_i)$ of the static Kohn-Sham response function (default: ‘1d-6’)
- **solve** matrix inversion methods to solve the OEP equation. The different options are: GESV, TSVD, GTSVD. GESV corresponds to a direct solution without any regularization technique. TSVD and GTSVD correspond to two solutions with regularization according to Eqs. (55) and (56) of Ref. [4], respectively. (default: 'GTSVD')
- **thr_solve** threshold used during matrix inversion to solve the OEP equation with TSVD and GTSVD methods. Note that the default threshold of 1d-99 results in the absence of regularization (default: ‘1d-99’)
- **plot_always** if set to $\neq$ 0, enable writing of data-files for plotting for every iteration. Otherwise, only final results are written. (default: '0')
- **plot_vc** if set to $\neq$ 0, enable writing of file with correlation potential (default: '0')
- **plot_vx** if set to $\neq$ 0, enable writing of file with exchange potential (default: '0')
- **plot_vxc** if set to $\neq$ 0,enable writing of file with exchange-correlation potential (default: '0')
- **plot_vref** if set to $\neq$ 0, enable writing of file with reference potential (default: '0')
- **plot_rho_ks** if set to $\neq$ 0, enable writing of file with Kohn-Sham density (default: '0')
- **plot_rho_ref** if set to $\neq$ 0, enable writing of file with reference reference (default: '0')
- **plot_rho_diff** if set to $\neq$ 0, enable writing of file with difference between Kohn-Sham and reference density (default: '0')
- **plot_x** if set to $\neq$ 0, enable writing of file with plotting data along x-axis (default: '0')
- **plot_y** if set to $\neq$ 0, enable writing of file with plotting data along y-axis (default: '0')
- **plot_z** if set to $\neq$ 0, enable writing of file with plotting data along z-axis (default: '0')
- **plot_xy** if set to $\neq$ 0, enable writing of file with plotting data for xy-plane (default: '0')
- **plot_yz** if set to $\neq$ 0, enable writing of file with plotting data for yz-plane (default: '0')
- **plot_xz** if set to $\neq$ 0, enable writing of file with plotting data for xz-plane (default: '0')
- **gridsize** determine the number of grid points for stored data (default: '2048')
- **plotrange** determine the range for which plotting data will be evaluated as [-plotrange:plotrange] (default: '20d0')
- **verb** determines the level of verbosity in the output file, integer values of 0, 1, 2, and 3 provide different levels of verbosity (default ’0’)
- **noa** number of electron in $\alpha$ spin channel, required for calculations of open-shell systems with KSINV.
- **nob** number of electron in $\beta$ spin channel, required for calculations of open-shell systems with KSINV.
- **vhoep** if set to 0, enable the calculation of the Hartree potential from the representation in the OEP basis instead of the construction from the density matrix as in the Hartree-Fock calculation (default: ‘0’)
- **space_sym** if set to 0, enable the space-symmetrization. When active sets vhoep=1 thr_sym=1d-10. (default: '0d0')
- **vref_fa_sameab** if set to $\neq$ 0, force the Fermi-Amaldi reference potential to be the same for $\alpha$ and $\beta$ spin channels in UKSINV calculations (default: ‘0’)
- **homo** if set to $\neq$ 0, enable the use of the HOMO condition. Note that epsilon_major/epsilon_minor must be specified to get meaningful results. (default '0')
- **epsilon_major** when homo=1, specifies the energy of HOMO orbital in $\alpha$ spin channel.
- **epsilon_minor** when homo=1, specifies the energy of HOMO orbital in $\beta$ spin channel.
- **density_test** if set to $\neq$ 0, enable density difference test. $|\rho_{ref}(r) - \rho_{KS}(r)|$ is evaluated on real space grid and integrated. (default '0')
- **thr_den_intgr** threshold for density difference integration test (default '1d-6')

Since KS correlation and exchange potentials are important in KS inversion, we provide an illustration of how to plot these quantities. Let us assume that we have performed calculations for CO with the following options:
```
acfd;ksinv,refden=1325.1,e_ref=-113.285493180105,thr_fai_oep=1.7d-2,\
plot_vx=1,plot_vc=1,plot_vxc=1,plot_vref=1,plot_z=1
```
At the end one has the files vref-final.z, vx-final.z, vc-final.z and vxc-final.z with reference, exchange, correlation and exchange-correlation potentials. The potentials can be plotted using Python and matplotlib as follows:
```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

vref = pd.read_csv('vref-total-final-z.csv')
vx = pd.read_csv('vx-total-final-z.csv')
vc = pd.read_csv('vc-total-final-z.csv')
vxc = pd.read_csv('vxc-total-final-z.csv')

plt.plot(vx.z, vref.value+vx.value, color='orangered', label='$v_x$')
plt.plot(vc.z, vc.value, color='orange', label='$v_c$')
plt.plot(vxc.z, vref.value+vxc.value, color='dodgerblue', label='$v_{xc}$')

plt.ylabel('Potential (Hartree)', fontsize=16)
plt.xlabel('r (Bohr)', fontsize=16)
plt.ylim(-2, 0.35)
plt.xlim(-5, 5)
plt.legend(frameon=False, fontsize=14)
plt.show()
```
![](ksinv_co.png)

CCSD(T) and FCI reference densities can be prepared in the same way as CCSD.

For CCSD(T):
```
{CCSD(T);core;dm,1325.1}
{matrop;load,den,den,1325.1;write,den,dm.dat,status=rewind,prec=sci}
```
For FCI:
```
{FCI;core;dm,1325.1}
{matrop;load,den,den,1325.1;write,den,dm.dat,status=rewind,prec=sci}
```
**Bibilography:**  
[1] J. Erhard, E. Trushin, A. Görling [J. Chem. Phys.](https://aip.scitation.org/doi/full/10.1063/5.0087356) 156, 204124 (2022)  
[2] A. Görling, [Phys. Rev. A](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.46.3753) 46, 3753 (1992)  
[3] E. Trushin, A. Görling, [J. Chem. Phys.](https://aip.scitation.org/doi/full/10.1063/5.0056431) 155, 054109 (2021)  
[4] P. Bleiziffer, A. Heßelmann, A. Görling [J. Chem. Phys.](https://doi.org/10.1063/1.4818984) 139, 084113 (2013)
