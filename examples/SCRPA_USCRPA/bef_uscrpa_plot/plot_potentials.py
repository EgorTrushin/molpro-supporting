#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def load_potential(filename):
    """Read data from file."""
    coord, vref, vrest, v = list(), list(), list(), list()
    for line in open(filename):
        aux = line.split()
        coord.append(float(aux[4]))
        vref.append(float(aux[6]))
        vrest.append(float(aux[7]))
        v.append(float(aux[8]))
    coord = np.array(coord)
    vref, vrest, v = np.array(vref), np.array(vrest), np.array(v)
    return coord, vref, vrest, v

coord, _, _, vxa = load_potential('vxa-final.z')
coord, _, vca, _ = load_potential('vca-final.z')

coord, _, _, vxb = load_potential('vxb-final.z')
coord, _, vcb, _ = load_potential('vcb-final.z')


plt.plot(coord, vxa+vca, color='dodgerblue', label=r'$v_{xc,\alpha}$')
plt.plot(coord, vxb+vcb, color='orangered', label=r'$v_{xc,\beta}$')
#plt.plot(coord, vx+vc, color='orange', label='$v_{xc}$')

plt.ylabel('Potential (Hartree)', fontsize=16)
plt.xlabel('r (Bohr)', fontsize=16)
#plt.ylim(-2, 0.25)
plt.xlim(-5, 5)
plt.legend(frameon=False, fontsize=14)
plt.tight_layout()
plt.show()