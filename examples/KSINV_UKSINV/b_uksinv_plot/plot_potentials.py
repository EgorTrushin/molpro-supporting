#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def load_potential(filename):
    """Read data from file."""
    coord = list()
    pot = list()
    for line in open(filename):
        aux = line.split()
        coord.append(float(aux[2]))
        pot.append(float(aux[3]))
    coord = np.array(coord)
    pot = np.array(pot)
    return coord, pot

coord, vref_a = load_potential('vref-alpha-final.z')
coord, vxc_a = load_potential('vxc-alpha-final.z')
coord, vref_b = load_potential('vref-beta-final.z')
coord, vxc_b = load_potential('vxc-beta-final.z')

plt.plot(coord, vref_a+vxc_a, color='orangered', label=r'$v_{xc}^{\alpha}$')
plt.plot(coord, vref_b+vxc_b, color='dodgerblue', label=r'$v_{xc}^{\beta}$')

plt.ylabel('Potential (Hartree)', fontsize=16)
plt.xlabel('r (Bohr)', fontsize=16)
plt.xscale("log")
plt.xlim(0, 20)
plt.legend(frameon=False, fontsize=14)
plt.tight_layout()
plt.show()
