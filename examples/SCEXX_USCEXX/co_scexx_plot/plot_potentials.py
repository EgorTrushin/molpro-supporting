#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def load_potential(filename):
    """Read data from file."""
    coord, vxref, vxrest, vx = list(), list(), list(), list()
    for line in open(filename):
        aux = line.split()
        coord.append(float(aux[4]))
        vxref.append(float(aux[6]))
        vxrest.append(float(aux[7]))
        vx.append(float(aux[8]))
    coord = np.array(coord)
    vxref, vxrest, vx = np.array(vxref), np.array(vxrest), np.array(vx)
    return coord, vxref, vxrest, vx

coord, vxref, vxrest, vx = load_potential('vx-final.z')

plt.plot(coord, vxref, color='orangered', label='$v_{x,ref}$')
plt.plot(coord, vxrest, color='dodgerblue', label='$v_{x,rest}$')
plt.plot(coord, vx, color='orange', label='$v_x$')

plt.ylabel('Potential (Hartree)', fontsize=16)
plt.xlabel('r (Bohr)', fontsize=16)
plt.xlim(-5, 5)
plt.legend(frameon=False, fontsize=14)
plt.tight_layout()
plt.show()