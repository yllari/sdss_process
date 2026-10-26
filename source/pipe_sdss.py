import numpy as np
from astropy.io import fits
from scipy.stats import norm
import matplotlib.pyplot as plt
import requests
import os
import time
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
        wave_filter = (self.wave_air > self.w_lims[0]) & (self.wave_air < self.w_lims[1])
        self.wave_air = self.wave_air[wave_filter]
        self.flux = self.flux[wave_filter]
        self.ivar = self.ivar[wave_filter]
        self.wdisp = self.wdisp[wave_filter]
        self.eflux = self.eflux[wave_filter]
        self.sigma_wave = self.sigma_wave[wave_filter]

        # recalculate max sigma wave for our new wave ranges
        self.max_sigma_wave = np.max(self.sigma_wave)

    def resample(self, spacing):
        min_wave = np.amin(self.wave_air)
        max_wave = np.amax(self.wave_air)
        new_wave = np.linspace(min_wave, max_wave, int( (max_wave-min_wave)/spacing) )

        new_flux = np.interp(new_wave, self.wave_air, self.flux)
        new_sigma_wave = np.interp(new_wave, self.wave_air, self.sigma_wave)
        new_eflux = np.interp(new_wave, self.wave_air, self.eflux)
        new_wdisp = np.interp(new_wave, self.wave_air, self.wdisp)

        self.wave_air = new_wave
        self.flux = new_flux

        self.ivar = 1/new_eflux
        self.eflux = new_eflux

        self.wdisp = new_wdisp
        self.sigma_wave = new_sigma_wave


    def __init__(self,
            filepath : str,
            Av : float,
            eAv : float,
            w_lims : list=[3800, 9000],
            interpolate=False,
            w_margin=15):
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

        self.wave_vac = 10**(data_fits[1].data['loglam'])
        self.wave_air = self.__vacuum_to_air(self.wave_vac)
        self.w_lims = w_lims
        # The flux is in 10^-17 erg/cm/s/A -> 10^{-16}
        self.flux = data_fits[1].data['flux']/10
        self.ivar = data_fits[1].data['ivar']
        self.wdisp = data_fits[1].data['wdisp']
        # ivar = 0 means infinite error, effectively bad pixel
        # for our purposes we can just give large values of eflux
        self.eflux = 1./np.sqrt(self.ivar)/10
        self.eflux[self.ivar == 0] = 999999  # large error
        # leave margin for degrade, by default 15 AA
        self.w_lims = [w_lims[0] - w_margin, w_lims[1] + w_margin]

        # Extinction values
        self.Av = Av
        self.eAv = eAv
        self.sigma_wave = self.sigma_wave_air()

        self.filter_wave()
        self.w_lims = w_lims
        if interpolate:
            self.resample(1)
        self.max_sigma_wave = np.max(self.sigma_wave)


    def _plot_debug(self):
        fig, ax = plt.subplots(1, 5)
        fig.set_size_inches(16,5)
        ax[0].plot(self.wave_air, self.resolution_power() )
        ax[0].set_xlabel(r'$\lambda$ air')
        ax[0].set_ylabel(r'Resolution power')

        ax[1].plot(self.wave_air, self.sigma_wave )
        ax[1].set_xlabel(r'$\lambda$ air')
        ax[1].set_ylabel(r'$\sigma_\lambda$ air')


        # Showing only the most affected range
        nodeg = self.flux 
        print('Integrated flux before degrade', np.sum(self.flux[:-1]*np.diff(self.wave_air)))
        self.degrade()
        print('Integrated flux after degrade', np.sum(self.flux[:-1]*np.diff(self.wave_air)))
        ax[2].plot(self.wave_air, nodeg - self.flux)
        ax[2].plot(self.wave_air, self.eflux)
        ax[2].set_xlabel(r'$\lambda$ air')
        ax[2].set_ylabel(r'flux - degraded flux')
        #ax[2].set_ylim([-10,10])

        filter_wave = (self.wave_air > 4000) & (self.wave_air < 4100)
        ax[3].plot(self.wave_air[filter_wave], nodeg[filter_wave], label='original')
        ax[3].plot(self.wave_air[filter_wave], self.flux[filter_wave], label='degraded')
        ax[3].set_ylim()
        ax[3].set_xlabel(r'$\lambda$ air')
        ax[3].set_ylabel(r'flux')
        ax[3].set_xlim([4000,4100])
        ax[3].legend()

        # Now the less affeted to compare
        filter_wave = (self.wave_air > 6900) & (self.wave_air < 7000)
        ax[4].plot(self.wave_air[filter_wave], self.flux[filter_wave], label='original')
        ax[4].plot(self.wave_air[filter_wave], nodeg[filter_wave], label='degraded')
        ax[4].set_ylim()
        ax[4].set_xlabel(r'$\lambda$ air')
        ax[4].set_ylabel(r'flux')
        ax[4].set_xlim([6900,7000])
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

        ## Alternative formula
        #sigma_wave_vac = 1e-4 * self.wdisp * np.log(10) * self.wave_vac
        #res = self.__vacuum_to_air_deriv(self.wave_vac)*sigma_wave_vac
        
        res = np.zeros(len(self.wdisp))
        res[:-1] = np.diff(self.wave_air)*self.wdisp[:-1]
        res[-1] = res[-2]

        return res


    def resolution_power(self):
        ''' Return the resolution power'''
        return self.wave_air/(2.35*self.sigma_wave)


    def give_data(self):
        ''' Returns all necessary info from read fits'''
        return self.wave_air, self.flux, self.eflux

    def corr_ext(self):
        ''' Correct for Milky Way's extinction '''
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
        max_sigma = np.max(sigma_wave) # hotfix to avoid maxsigma zero

        conv_sigma = np.sqrt(max_sigma**2 - sigma_wave**2)
        conv_sigma[conv_sigma < 1e-9] = 1e-9 # dividing by zero fix
        pad_sigma[:pad] = conv_sigma[0]
        pad_sigma[-pad:] = conv_sigma[-1]
        pad_sigma[pad:-pad] = conv_sigma


        pad_flux[:pad] = self.flux[0]
        pad_flux[-pad:] = self.flux[-1]
        pad_flux[pad:-pad] = self.flux

        xy = np.meshgrid(pad_wave, pad_sigma)

        ''' WARNING, IN REALITY WE NEED  DX-DEPENDENT NORMALIZATION if it is not resampled'''
        pad_mask = norm.pdf(xy[0].T,  loc=xy[0], scale=xy[1].T)
        dx = np.diff(pad_wave)
        dx = np.append(dx, dx[-1])
        pad_mask_norm = np.sum(pad_mask*dx[None, :], axis=1)
        pad_mask = pad_mask/pad_mask_norm[None,:] # broadcasting to normalize kernel
        
        gauss_flux = np.dot(pad_flux, pad_mask*dx[None, :])

        #from astropy.convolution import Gaussian1DKernel, convolve
        #gaus = Gaussian1DKernel(stddev=5)
        #gauss_flux = convolve(pad_flux,  Gaussian1DKernel(stddev=np.sqrt(max_sigma**2 - np.mean(sigma_wave)**2)),boundary='extend', nan_treatment='interpolate')

        # Update values
        self.wave_air = pad_wave[pad:-pad]
        self.flux = gauss_flux[pad:-pad]
        print(print(self.flux))
        self.sigma_wave = np.ones(sigma_wave.shape)*max_sigma
        self.max_sigma_wave = max_sigma

        return self.wave_air, self.flux


if __name__ == '__main__':
    def single_spectrum_save(objid, interpolate):
        import pandas as pd
        from astropy.io import ascii as astro_ascii
        header_path = 'data/voids_header.csv'
        ext_path = 'data/voids_extinction.txt'
        gal_type = 'voids' # void, clusters or filaments
        w_lims = [3800, 7100]
        RV = 3.1

        # reading data and ext_table
        header_df = pd.read_csv(header_path)
        ext_table = astro_ascii.read(ext_path)

        # voids_extinction and voids_header have same ordering of galaxies (indexed)
        gal_id = objid
        i = np.where(header_df['galaxy'] == objid)
        z = header_df['z'].iloc[i]
        filepath = f'data/voids/CAVITY{gal_id}.fits'
        print('Processing galaxy', objid, 'id', gal_id)

        EBV = ext_table['mean_E_B_V_SandF'][i]
        # There is another option, exact position, but without std
        # AV = ext_table['AV_SandF'][i]
        stdEBV = ext_table['stdev_E_B_V_SandF'][i]
        AV = EBV*RV
        eAV = stdEBV*RV
        print('AV', AV)

        sdss = pipe_sdss(filepath, AV, eAV, w_lims, interpolate=interpolate)

        ''' Correct for galactic extinction, degrade the specturm to uniform sigma and filter for wavelenghts
        given by wave_lims'''
        sdss.corr_ext()
        sdss.degrade()
        sdss.filter_wave()

        indx = np.arange(1, len(sdss.wave_air) + 1, 1) # First index necessary for PyPipe3D spectrum formatting
        print('sigma_inst', sdss.max_sigma_wave)
        np.savetxt(f'test_{gal_id}.txt', np.array([indx,sdss.wave_air,sdss.flux, sdss.eflux]).T)
        np.savetxt(f'test_{gal_id}_data.txt', [sdss.max_sigma_wave, z])
    #single_spectrum_save(3534, True)
    filepath = f'data/voids/CAVITY137.fits'
    w_lims = [3800, 7100]

    # For CAVITY 137
    # SDF 0.0492, 0.0009
    # SandF 0.0423, 0.0008
    t_0 = time.time()
    for i in range(1):
        print(i)
        sdss = pipe_sdss(filepath, 0.0492, 0.0009,w_lims, interpolate=False)
        #sdss.degrade()
        sdss._plot_debug()

    print((time.time() - t_0)/100)


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
