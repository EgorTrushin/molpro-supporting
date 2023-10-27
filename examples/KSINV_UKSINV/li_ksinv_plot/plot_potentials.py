#!/usr/bin/env python3

import sys
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

coord, vref = load_potential('vref-total-final.z')
coord, vx = load_potential('vx-total-final.z')
coord, vc = load_potential('vc-total-final.z')
coord, vxc = load_potential('vxc-total-final.z')

plt.plot(coord, vref+vx, color='orangered', label='$v_x$')
plt.plot(coord, vc, color='dodgerblue', label='$v_c$')
plt.plot(coord, vref+vxc, color='lawngreen', label='$v_{xc}$')

plt.ylabel('Potential (Hartree)', fontsize=16)
plt.xlabel('r (Bohr)', fontsize=16)
plt.xscale("log")
plt.xlim(0, 20)
plt.legend(frameon=False, fontsize=14)
plt.tight_layout()
plt.show()
