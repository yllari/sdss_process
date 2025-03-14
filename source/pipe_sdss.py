import numpy as np
from astropy.io import fits
from scipy.stats import norm
import matplotlib.pyplot as plt
import requests
import os
from ext_utils import correct_extinction

class pipe_sdss():
    def __vacuum_to_air(self, wave):
        ''' Vacuum to air formula for the wavelength in SDSS DR7
              AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4)
            from https://classic.sdss.org/dr7/products/spectra/vacwavelength.php
        :param wave: wavelength in vacuum
        :return: wavelength in air
        '''
        return wave / (1.0 + 2.735182e-4 + 131.4182 / wave**2 + 2.76249e8 / wave**4)


    def __vacuum_to_air_deriv(self, wave):
        ''' Return the evaulated derivate from the formula __vacum to_air_deriv in the given vacuum wavelength
        :param wave: wavelength in vacuum
        :return: lambda_air derivative
        '''
        deriv = self.__vacuum_to_air(wave)/wave
        deriv += (-2*131.4182/wave**3 - 4*2.76249e8/wave**5)*wave/(1.0 + 2.735182e-4 + 131.4182 / wave**2 + 2.76249e8 / wave**4)
        return deriv


    def filter_wave(self):
        wave_filter = (self.wave_air> self.w_lims[0]) & (self.wave_air < self.w_lims[1])
        self.wave_vac = self.wave_vac[wave_filter]
        self.wave_air = self.wave_air[wave_filter]
        self.logwave = self.logwave[wave_filter]
        self.flux = self.flux[wave_filter]
        self.ivar = self.ivar[wave_filter]
        self.wdisp = self.wdisp[wave_filter]
        self.eflux = self.eflux[wave_filter]
        self.sigma_wave = self.sigma_wave[wave_filter]
        self.max_sigma_wave = np.max(self.sigma_wave)


    def __init__(self,
            filepath : str,
            Av : float,
            eAv : float,
            w_lims : list = [3800, 9000],
            w_margin = 15):
        ''' read spectrum from SDSS fits file
        :param filepath: filepath, default string from which read the SDSS spectrum
        :param Av: mean extinction in V band for extinction correction
        :param eAv: std extinction in V band
        :param w_lims: wavelength  limits for spectrum
        :param w_margin: some margin to leave to avoid boundary problems when processing, in AAº
        '''

        # Documentation for header:
        # https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html
        data_fits = fits.open(filepath)

        self.wave_vac= 10**(data_fits[1].data['loglam'])
        self.wave_air = self.__vacuum_to_air(self.wave_vac)
        self.w_lims = w_lims
        self.logwave = data_fits[1].data['loglam']
        self.flux = data_fits[1].data['flux']
        self.ivar = data_fits[1].data['ivar']
        self.wdisp = data_fits[1].data['wdisp']
        # ivar = 0 means infinite error, effectively bad pixel
        # for our purposes we can just give large values of eflux
        self.eflux = 1./np.sqrt(self.ivar)
        self.eflux[self.ivar == 0] = 999999  # large error
        # leave margin for degrade, by default 15 AA
        self.w_lims = [w_lims[0] - w_margin, w_lims[1] + w_margin]


        # Extinction values
        self.Av = Av
        self.eAv = eAv
        self.sigma_wave = self.sigma_wave_air()

        self.filter_wave()
        self.w_lims = w_lims
        self.max_sigma_wave = np.max(self.sigma_wave)


    def _plot_debug(self):
        fig, ax = plt.subplots(1, 5)
        fig.set_size_inches(16,5)
        ax[0].plot(self.wave_air, self.resolution_power() )
        ax[0].set_xlabel(r'$\lambda$ air')
        ax[0].set_ylabel(r'Resolution power')

        ax[1].plot(self.wave_air, self.sigma_wave_air() )
        ax[1].set_xlabel(r'$\lambda$ air')
        ax[1].set_ylabel(r'$\sigma_\lambda$ air')

        deg = self.degrade()
        ax[2].plot(self.wave_air, self.flux - deg[1] )
        ax[2].set_xlabel(r'$\lambda$ air')
        ax[2].set_ylabel(r'flux - degraded flux')
        ax[2].set_ylim([-10,10])

        # Showing only the most affected range
        filter_wave = (self.wave_air > 4000) & (self.wave_air < 4100)
        ax[3].plot(self.wave_air[filter_wave], self.flux[filter_wave], label='original')
        ax[3].plot(self.wave_air[filter_wave], deg[1][filter_wave], label='degraded')
        ax[3].set_ylim()
        ax[3].set_xlabel(r'$\lambda$ air')
        ax[3].set_ylabel(r'flux')
        ax[3].set_xlim([4000,4100])
        ax[3].legend()

        # Now the less affeted to compare
        filter_wave = (self.wave_air > 6500) & (self.wave_air < 6600)
        ax[4].plot(self.wave_air[filter_wave], self.flux[filter_wave], label='original')
        ax[4].plot(self.wave_air[filter_wave], deg[1][filter_wave], label='degraded')
        ax[4].set_ylim()
        ax[4].set_xlabel(r'$\lambda$ air')
        ax[4].set_ylabel(r'flux')
        ax[4].set_xlim([6500,6600])
        ax[4].legend()
        fig.tight_layout()
        fig.savefig('debug.png')

    def sigma_wave_air(self):
        ''' degrade all wavelength dispersion (i.e, gaussian) to its
        max value

        The sigma in logwave is expressed in terms of 0.001 loglam (1 pixel). From error propagation, we obtain the formula
        sigma_\lambda = \sigma_f \ln(10) \cdot 10^{\loglam}
        where \sigma_f =  1e-4 * wdisp

        We have to correct for the air formula as well, with the derivative
        '''

        sigma_wave_vac = 1e-4 * self.wdisp * np.log(10) * self.wave_vac
        sigma_wave_air = self.__vacuum_to_air_deriv(self.wave_vac)*sigma_wave_vac

        return sigma_wave_air


    def resolution_power(self):
        return self.wave_air/(2.35*self.sigma_wave_air())


    def give_data(self):
        ''' Returns all necessary info from read fits'''
        return self.wave_air, self.flux, self.eflux

    def corr_ext(self):
        self.flux, self.eflux = correct_extinction(self.wave_air, self.flux, self.eflux, self.Av, self.eAv)
        return None

    def degrade(self, pad=10):
        ''' Degrade the spectrum quality to achieve uniform resolution '''

        pad_wave = np.zeros( (int(len(self.wave_air) + 2*pad, ) ) )
        pad_sigma = np.zeros( (int(len(self.wave_air) + 2*pad, ) ) )
        pad_flux = np.zeros( (int(len(self.wave_air) + 2*pad, ) ) )

        low_wave = self.wave_air[0]
        high_wave = self.wave_air[-1]
        low_dx_wave = self.wave_air[1] - low_wave # be careful with nan
        high_dx_wave =  high_wave - self.wave_air[-2] # be careful with nan
        pad_wave[:pad+1] = np.linspace(low_wave - low_dx_wave*pad, low_wave, pad+1) # boundary='extend' behaviour of astropy
        pad_wave[-pad-1:] = np.linspace(high_wave, high_wave + high_dx_wave*pad, pad+1) # boundary='extend' behaviour of astropy
        pad_wave[pad:-pad] = self.wave_air


        sigma_wave = self.sigma_wave
        max_sigma = np.nanmax(sigma_wave) # hotfix to avoid maxsigma zero

        conv_sigma = np.sqrt(max_sigma**2 - sigma_wave**2)
        conv_sigma[conv_sigma < 1e-9] = 1e-9 # dividing by zero fix
        pad_sigma[:pad] = conv_sigma[0]
        pad_sigma[-pad:] = conv_sigma[-1]
        pad_sigma[pad:-pad] = conv_sigma


        pad_flux[:pad] = self.flux[0]
        pad_flux[-pad:] = self.flux[-1]
        pad_flux[pad:-pad] = self.flux

        xy = np.meshgrid(pad_wave, pad_sigma)

        ''' WARNING, IN REALITY WE NEED  DX-DEPENDENT NORMALIZATION'''
        pad_mask = norm.pdf(xy[0].T,  loc=xy[0], scale=xy[1].T)
        pad_mask_norm = np.sum(pad_mask, axis=1)
        pad_mask = pad_mask/pad_mask_norm[None,:] # broadcasting
        gauss_flux = pad_flux@pad_mask

        #from astropy.convolution import Gaussian1DKernel, convolve
        #gaus = Gaussian1DKernel(stddev=5)
        #gauss_flux = convolve(pad_flux,  Gaussian1DKernel(stddev=np.sqrt(max_sigma**2 - np.mean(sigma_wave)**2)),boundary='extend', nan_treatment='interpolate')

        # Update values
        self.wave_air = pad_wave[pad:-pad]
        self.flux = gauss_flux[pad:-pad]
        self.sigma_wave = np.ones(sigma_wave.shape)*max_sigma
        self.max_sigma_wave = max_sigma

        return self.wave_air, self.flux


if __name__ == '__main__':
    filepath = 'sdss_process/data/voids/CAVITY137.fits'
    w_lims = [3800, 7100]

    # For CAVITY 137
    # SDF 0.0492, 0.0009
    # SandF 0.0423, 0.0008
    sdss = pipe_sdss(filepath, 0.0492, 0.0009,w_lims)
    data_0 = sdss.give_data()
    sdss.corr_ext()
    sdss.degrade()

    data_1 = sdss.give_data()
    plt.plot(data_1[0], data_1[1], alpha=.5)
    data_ref = np.genfromtxt('sdss_process/data/CAVITY137_spec_sdss_conv_nonan.txt')
    plt.plot(data_ref[:,1], data_ref[:,2], alpha=0.5)
    lims = [6500, 7000]
    #lims = [4300, 5000]
    lims = [5800, 6200]
    plt.xlim(lims)
    filter_w = (data_ref[:,1] > lims[0]) & (data_ref[:,1] < lims[1])
    plt.ylim([np.min(data_ref[:,2][filter_w]), np.max(data_ref[:,2][filter_w])] )
    plt.savefig('spectrum.png')
