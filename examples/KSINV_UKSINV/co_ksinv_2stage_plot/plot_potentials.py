#!/usr/bin/env python3

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
plt.tight_layout()
plt.show()
