import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt

paths = os.listdir('voids/')

vdisp = np.zeros(len(paths))

i = 0
for path in paths:
    print(i, path)
    vdisp[i] = fits.open('voids/' + path)[2].data['VDISP']
    i+=1
plt.hist(vdisp, bins=100)
plt.show()
