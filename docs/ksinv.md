## KSINV program
The KSINV program allows ... spin-restricted

```
memory,2000,m ! memory specification

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

hf,maxit=30 ! HF calculation

{CCSD;core;dm,1325.1} ! reference CCSD calculation
e_ref_ccsd=energy ! save energy to provide as a reference to KSINV

acfd;ksinv,refden=1325.1,e_ref=e_ref_ccsd,thr_fai_oep=1.7d-2 ! KSINV calculation
```
... bla-bla ...
```
memory,2000,m ! memory specification

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
{matrop;export,1325.1,dm.dat,status=rewind,prec=sci}
```
... bla-bla ...
```
memory,2000,m ! memory specification

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

{matrop;import,1325.1,file=dm.dat} ! read stored density matrix

acfd;ksinv,refden=1325.1,e_ref=-113.285493180105,thr_fai_oep=1.7d-2 ! KSINV calculations
```

- **refden** record from which the reference density is read
- **orb** record from which the occupation numbers are read (default: ‘2100.2’)
- **save** record in which the resulting orbital coefficients, eigenvalues, etc. are written (default: '2101.2')
- **e_ref** total energy of calculations from which reference density is provided
- **thr_inversion**
- **maxit** maximum number of iterations (default '30')
- **minit** minimum number of iterations (default '3')
- **mixing_a**
- **mix_switch_iter**
- **backfilter** ... (default '1')
- **thr_sym** threshold for symmetrization of the OEP basis set to enforce OEP basis exhibits full symmetry of molecule. Set the threshold to 1d-10 to enable symmetrization. (default: ‘0d0’)
- **thr_fai_oep** threshold for processing OEP basis according to Section IIB5 in Ref. [2] (default: ‘1.7d-2’)
- **thr_overlap_oep** threshold for processing OEP basis according to Section IIB2 in Ref. [2] (default: ‘1d-99’)
- **vref_fa** if set to $\neq$ 0, enable the use of the Fermi-Amaldi potential as reference potential when applying the charge condition. Otherwise, the reference potential is constructed according to Eq. (45) of Ref. [2] (default: '1')
- **thr_oep** threshold for throwing out contributions corresponding to small eigenvalue differences appearing in the denominator when constructing the so-called lambda term $1/(\varepsilon_a - \varepsilon_i)$ of the static Kohn-Sham response function (default: ‘1d-6’)
- **solve** matrix inversion methods to solve the OEP equation. The different options are: GESV, TSVD, GTSVD. GESV corresponds to a direct solution without any regularization technique. TSVD and GTSVD correspond to two solutions with regularization according to Eqs. (55) and (56) of Ref. [3], respectively. (default: 'GTSVD')
- **thr_solve** threshold used during matrix inversion to solve the OEP equation with TSVD and TGSVD methods. Note that the default threshold of 1d-99 results in the absence of regularization (default: ‘1d-99’)
- **plot_always** if set to $\neq$ 0, enable plotting for every iteration. Otherwise, only final results are plotted. (default: '0')
- **plot_vc** if set to $\neq$ 0, enable writing of file with correlation potential (default: '0')
- **plot_vx** if set to $\neq$ 0, enable writing of file with exchange potential (default: '0')
- **plot_vxc** if set to $\neq$ 0,enable writing of file with exchange-correlation potential (default: '0')
- **plot_vref** if set to $\neq$ 0, enable writing of file with reference potential (default: '0')
- **plot_rho_ks** if set to $\neq$ 0, enable writing of file with Kohn-Sham density (default: '0')
- **plot_rho_ref** if set to $\neq$ 0, enable writing of file with reference reference (default: '0')
- **plot_rho_diff** if set to $\neq$ 0, enable writing of file with difference between Kohn-Sham and reference density (default: '0')
- **plot_x** if set to $\neq$ 0, enable writing of file data along x-axis (default: '0')
- **plot_y** if set to $\neq$ 0, enable writing of file data along y-axis (default: '0')
- **plot_z** if set to $\neq$ 0, enable writing of file data along z-axis (default: '0')
- **gridsize**
- **plotrange**
- **verb** determines the level of verbosity in the output file, integer values of 0, 1, 2, and 3 provide different levels of verbosity (default ’0’)

plotting example

[1] J. Erhard, E. Trushin, A. Görling [J. Chem. Phys.](https://aip.scitation.org/doi/full/10.1063/5.0087356) 156, 204124 (2022)  
[2] E. Trushin, A. Görling, [J. Chem. Phys.](https://aip.scitation.org/doi/full/10.1063/5.0056431) 155, 054109 (2021)  
[3] P. Bleiziffer, A. Heßelmann, A. Görling [J. Chem. Phys.](https://doi.org/10.1063/1.4818984) 139, 084113 (2013)