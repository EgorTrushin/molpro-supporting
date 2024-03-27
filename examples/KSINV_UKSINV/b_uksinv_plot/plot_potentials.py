#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

vref_a = pd.read_csv('vref-alpha-final-z.csv')
vref_b = pd.read_csv('vref-beta-final-z.csv')
vxc_a = pd.read_csv('vxc-alpha-final-z.csv')
vxc_b = pd.read_csv('vxc-beta-final-z.csv')

plt.plot(vxc_a.z, vref_a.value+vxc_a.value, color='orangered', label=r'$v_{xc}^{\alpha}$')
plt.plot(vxc_b.z, vref_b.value+vxc_b.value, color='dodgerblue', label=r'$v_{xc}^{\beta}$')

plt.ylabel('Potential (Hartree)', fontsize=16)
plt.xlabel('r (Bohr)', fontsize=16)
plt.xscale("log")
plt.xlim(0, 20)
plt.legend(frameon=False, fontsize=14)
plt.tight_layout()
plt.show()
