import numpy as np
from astropy.io import fits
from scipy.stats import norm
import matplotlib.pyplot as plt

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


    def __filter_wave(self):
        wave_filter = (self.wave_air> self.w_lims[0]) & (self.wave_air < self.w_lims[1])
        self.wave_vac = self.wave_air[wave_filter]
        self.wave_air = self.wave_air[wave_filter]
        self.logwave = data_fits[1].data['loglam'][wave_filter]
        plt.show()
        self.flux = data_fits[1].data['flux'][wave_filter]
        self.ivar = data_fits[1].data['ivar'][wave_filter]
        self.wdisp = data_fits[1].data['wdisp'][wave_filter]
        self.eflux = 1./np.sqrt(self.ivar)


    def __init__(self, filepath, w_lims):

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
        self.eflux = 1./np.sqrt(self.ivar)

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

    def correct_extinction(self):
        ...

    def degrade(self, pad=10):
        result_flux = 3

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


        sigma_wave = self.sigma_wave_air()
        max_sigma = np.nanmax(sigma_wave) # hotfix to avoid maxsigma zero

        conv_sigma = np.sqrt(max_sigma**2 - sigma_wave**2)
        conv_sigma[conv_sigma < 0.0001] = 0.0001
        pad_sigma[:pad] = conv_sigma[0]
        pad_sigma[-pad:] = conv_sigma[-1]
        pad_sigma[pad:-pad] = conv_sigma


        pad_flux[:pad] = self.flux[0]
        pad_flux[-pad:] = self.flux[-1]
        pad_flux[pad:-pad] = self.flux

        dwave = np.zeros((len(pad_flux,)))
        dwave[0:-1] = np.diff(pad_wave)
        dwave[-1] = dwave[-2]
        print(dwave)

        xy = np.meshgrid(pad_wave, pad_sigma)


        #pad_mask = norm.pdf(xy[0].T,  loc=xy[0], scale=xy[1].T)



        ''' WARNING, WIP DX-DEPENDENT NORMALIZATION'''
        pad_mask_norm = np.sum(pad_mask, axis=1)
        #print(pad_mask.shape)
        #min_indx = np.where(pad_mask == np.max(pad_mask))
        #min_indx = min_indx[1][0]
        #plt.imshow(pad_mask > 0)


        from astropy.convolution import Gaussian1DKernel, convolve
        gaus = Gaussian1DKernel(stddev=5)
        print(gaus[0])
        gauss_flux = convolve(self.flux,  Gaussian1DKernel(stddev=np.sqrt(max_sigma**2 - np.mean(sigma_wave)**2)),boundary='extend', nan_treatment='interpolate')


        return self.wave_air, gauss_flux


if __name__ == '__main__':
    filepath = 'data/voids/CAVITY137.fits'
    w_lims = [3800, 7100]
    sdss = pipe_sdss(filepath, w_lims)
    data = sdss.give_data()
    deg  = sdss.degrade()
    plt.plot(*deg)
    plt.plot(sdss.wave_air,sdss.flux)
    #plt.plot(sdss.wave_air, sdss.flux)
    plt.savefig('example.png')
    def plot_resolution(class_data : pipe_sdss):
        fig, ax = plt.subplots()
        ax.set_xlim([3800,9000])
        ax.set_ylim([1000,3000])
        ax.plot(class_data.wave_air, class_data.resolution_power() )
        fig.savefig('resolution.png')
