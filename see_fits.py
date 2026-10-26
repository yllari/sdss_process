from astropy.io import fits
import matplotlib.pyplot as plt

d = fits.open('data/clusters/SDSS_spec303993281753147392.fits')
print(d[0].header)
plt.plot(10**d[1].data['loglam'], d[1].data['flux'])
plt.show()