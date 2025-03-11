import os
import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.interpolate import interp1d
import xml.etree.ElementTree as ET
import extinction
import matplotlib.pyplot as plt
from astropy.convolution import Box1DKernel, Gaussian1DKernel, convolve
import warnings

warnings.filterwarnings("ignore")

table2 = pd.read_table('clusters_list.txt', delim_whitespace=True, names=['names'])

table = pd.read_csv(f'clusters.csv') 

#print(table.galaxy)
#print(f'{table["VGS"][0]:02d}')
#print(table["DEJ2000"][0])
#print(table["DEJ2000"][0].split()[0]+'d'+table["DEJ2000"][0].split()[1]+'m'+table["DEJ2000"][0].split()[2]+'s')

d2 = {'#name': [], 'ra': [], 'dec': [], 'redshift': [], 'sigma_inst': []}
#for index, item in enumerate(table['galaxy']):
for index, item in enumerate(table['specobjid']):
#    if f'CAVITY{str(item)}' not in list(table2['names']) or os.path.isfile(f'spectra/convolved/cavity/CAVITY{str(item)}_spec_sdss_conv.txt') == True: #os.path.isfile(f'spectra/convolved/cavity/CAVITY{str(item)}_spec_sdss_conv.txt'): f'CAVITY_{item:02d}'
#    if f'SDSS_spec{str(item)}' not in list(table2['names']) or os.path.isfile(f'spectra/convolved/filaments/SDSS_spec{str(item)}_spec_sdss_conv.txt') == True: 
    if f'SDSS_spec{str(item)}' not in list(table2['names']) or os.path.isfile(f'spectra/convolved/clusters/SDSS_spec{str(item)}_spec_sdss_conv.txt') == True: 
        continue
    print(f'SDSS_spec{str(item)}')
#    print(f'CAVITY{str(item)}')
    #######
    ### Help: load data from fits tables
    #######
#    t = Table.read(f'./spectra/cavity/CAVITY{str(item)}.fits', hdu=1)
#    t = Table.read(f'./spectra/filaments/SDSS_spec{str(item)}.fits', hdu=1)
    t = Table.read(f'./spectra/clusters/SDSS_spec{str(item)}.fits', hdu=1)
#    if item == 'UGC01087':
#        flux = t['FLUX']
#        logwv = t['LOGLAM']
#        ivar = t['IVAR']
#        wdisp = t['WDISP']
#    else:
    flux = t['flux']
    logwv = t['loglam']
    ivar = t['ivar']
    wdisp = t['wdisp']
        
    if all(x == 0. for x in list(wdisp)):
        print('ERROR!!!:', item, 'has wdisp = 0')
        sigmainst.extend([np.nan])
        continue
    wdisp[wdisp==0.0] = np.nan
    
    #######
    ### Help: derivation of flux errors, conversion between vacuum and air wavelengths, and transformation of all magnitudes to linear wavelength scale
    #######
    eflux = np.sqrt(1./ivar)
    linwv_vac = 10**logwv
    linwv_air = (linwv_vac)/(1.0+2.735182E-4 + 131.4182/(linwv_vac**2)+2.76249E8/(linwv_vac**4))
    wdispAA = np.array([item*(linwv_vac[k+1]-linwv_vac[k]) if k<len(wdisp)-1 else item*(linwv_vac[k]-linwv_vac[k-1]) for k,item in enumerate(wdisp)]) #to have the sigma in AA
    f = interp1d(linwv_air, flux, kind='linear') #to have a regular wavelength step in AA
    wvlin = np.arange(int(np.nanmin(linwv_air)+1.), int(np.nanmax(linwv_air)), 1)
    fluxlin = f(wvlin)
    z0 = interp1d(linwv_air[np.where(np.isfinite(np.array(wdispAA)))], np.array(wdispAA)[np.where(np.isfinite(wdispAA))], kind='linear') # this is to interpolate the nan in the array to cover the whole spectral range (no holes)
    wdispAAint = z0(linwv_air)
    z = convolve(np.array(wdispAAint), Box1DKernel(width=20), boundary='extend', nan_treatment='interpolate') #to smooth the function
    f2 = interp1d(linwv_air, z, kind='linear', fill_value='extrapolate') #to have a regular wavelength step in AA
    wdispg = f2(wvlin)
    
    f3 = interp1d(linwv_air, eflux, kind='linear') #to have a regular wavelength step in AA
    efluxlin = f3(wvlin)
    
    #######
    ### Help: we obtain the corresponding galactic extinction based on the RA and Dec of the galaxy and correct the fluxes from such extinction based on Cardelli curve
    #######
#    ra = table["ra"][index].split()[0]+'h'+table["RAJ2000"][index].split()[1]+'m'+table["RAJ2000"][index].split()[2]+'s'
    ra = str(table["ra"][index])
#    de = table["DEJ2000"][index].split()[0]+'d'+table["DEJ2000"][index].split()[1]+'m'+table["DEJ2000"][index].split()[2]+'s'
    de = str(table["dec"][index])
#    if os.path.isfile(f'spectra/convolved/cavity/ext_CAVITY{str(item)}.xml') == False:
#        os.system("wget -nd 'http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr="+str(ra)+"+"+str(de)+"+equ+j2000' -O spectra/convolved/cavity/ext_CAVITY"+f"{str(item)}"+".xml")
#    if os.path.isfile(f'spectra/convolved/filaments/ext_SDSS_spec{str(item)}.xml') == False:
#        os.system("wget -nd 'http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr="+str(ra)+"+"+str(de)+"+equ+j2000' -O spectra/convolved/filaments/ext_SDSS_spec"+f"{str(item)}"+".xml")
    if os.path.isfile(f'spectra/convolved/clusters/ext_SDSS_spec{str(item)}.xml') == False:
        os.system("wget -nd 'http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr="+str(ra)+"+"+str(de)+"+equ+j2000' -O spectra/convolved/clusters/ext_SDSS_spec"+f"{str(item)}"+".xml")
#    tree = ET.parse(f'spectra/convolved/cavity/ext_CAVITY{str(item)}.xml')
#    tree = ET.parse(f'spectra/convolved/filaments/ext_SDSS_spec{str(item)}.xml')
    tree = ET.parse(f'spectra/convolved/clusters/ext_SDSS_spec{str(item)}.xml')
    root = tree.getroot()
    ext_attr = str(root[1][2][4])
    ext_attr = str(ext_attr.split()[1])[1:-1]
    if ext_attr == 'meanValueSFD':
        Av = root[1][2][4].text
        Av = float(Av.rstrip('\n').split()[0])
        eAv = root[1][2][6].text
        eAv = float(eAv.rstrip('\n').split()[0])
    else: 
        print('Error extracting Av from xml file')
        break
#    wv0 = wvlin/(1+table['redshift'][index])  # we do not need to correct from redshift. THIS HAS TO BE COMMENTED
#    print(table['redshift'][index])
    
    flux_corr = extinction.remove(extinction.ccm89(wvlin.astype('double'), Av, 3.1), fluxlin)
    eflux_corr1 = extinction.remove(extinction.ccm89(wvlin.astype('double'), Av, 3.1), efluxlin)
    eflux_corr = np.sqrt(eflux_corr1**2 + flux_corr**2 * (0.4*np.log(10)*eAv)**2)

    
    #######
    ### Help: this is to plot the instrumental dispersion/resolution power curves and compare them between each other and with Fig.35 of Smee et al. 2013, AJ, 146
    #######
#    plt.scatter(linwv_air, wdisp, marker='.', s=1., zorder=10)
#    plt.scatter(linwv_air, wdispAA, marker='s', s=1., zorder=10)
#    plt.plot(linwv_air, wdispAAint, '--', lw=1., zorder=20)
#    plt.scatter(linwv_air, linwv_air/(2.355*wdispAA), marker='.', s=1., zorder=10)
#    plt.plot(linwv_air, z, '-', lw=1., zorder=20)
#    plt.plot(wvlin, f2(wvlin), '-', lw=1., zorder=20)
    
    #######
    ### Help: this is to degrade the spectrum and have an instrumental dispersion that is the same along the whole spectral range (we need that in order to apply PIPE3D)
    #######
    wmin, wmax = 3800, 7100    ### Help: quick convolve taking the whole spectrum and the mean sigma within the waverange (3800, 7100)
    wrange = (wvlin >= wmin) & (wvlin <= wmax)
    worst_sigma = np.nanmax(wdispg[wrange])
#    sigma_mean = f2((wmax-wmin)/2.)
    sigma_mean = (np.nanmax(wdispg[wrange]) + np.nanmin(wdispg[wrange]))/2.
    sigma = np.sqrt(worst_sigma**2-sigma_mean**2)/(wvlin[1]-wvlin[0])
    flux_final = convolve(flux_corr[wrange], Gaussian1DKernel(stddev=sigma), boundary='extend', nan_treatment='interpolate')
    
    flux_final2 = np.zeros_like(flux_final)  ### Help: better procedure deriving the mean sigma in windows of 50AA and convolving the spectrum by pieces using those small windows
    n = 0
    for k in range(wmin,wmax,50):
        wmin2, wmax2 = k, k+50
        wrange2 = (wvlin >= wmin2) & (wvlin < wmax2)
        if len(wvlin[wrange2])>0:
            sigma_mean2 = (np.nanmax(wdispg[wrange2]) + np.nanmin(wdispg[wrange2]))/2.
            sigma2 = np.sqrt(worst_sigma**2-sigma_mean2**2)/(wvlin[1]-wvlin[0])
            flux_final2[n:n+len(wvlin[wrange2])] = convolve(flux_corr[wrange2], Gaussian1DKernel(stddev=sigma2), boundary='extend', nan_treatment='interpolate')
        n += len(wvlin[wrange2])
    
    #######
    ### Help: this is to plot the different spectra (observed, extinction-corrected, convolved...)
    #######
#    plt.figure(1, figsize=(8,4))
#    plt.grid(b=None, which='major', axis='both', zorder=0)
#    plt.plot(wvlin[wrange], fluxlin[wrange],'k-', lw=1.0)
#    plt.plot(wvlin[wrange], flux_corr[wrange],'r-', lw=1.0)
#    plt.plot(wvlin[wrange], flux_final,'b--', lw=1.0)
#    plt.plot(wvlin[wrange], flux_final2,'g--', lw=1.0)
#    plt.show()
#    
    
    #######
    ### Help: this is to save the convolved (and extinction-corrected) spectrum in a txt file
    #######
    
    d = {'index': np.arange(1,len(wvlin[wrange])+1), 'wavelength': wvlin[wrange], 'flux': flux_final2, 'eflux': np.array(eflux_corr)[wrange]}
    df = pd.DataFrame(data=d)
#    print(df)
    fmt = '%d', '%d', '%1.10f', '%1.10f'
#    np.savetxt(f'spectra/convolved/cavity/CAVITY{str(item)}_spec_sdss_conv.txt', df.values, fmt=fmt)
#    np.savetxt(f'spectra/convolved/filaments/SDSS_spec{str(item)}_spec_sdss_conv.txt', df.values, fmt=fmt)
    np.savetxt(f'spectra/convolved/clusters/SDSS_spec{str(item)}_spec_sdss_conv.txt', df.values, fmt=fmt)    
    
#######
### Help: this is to save the final instrumental dispersion of the spectra
#######
#    d2['#name'].extend([f'VGS_{item:02d}'])
#    d2['#name'].extend([f'CAVITY{str(item)}'])
    d2['#name'].extend([f'SDSS_spec{str(item)}'])
    d2['ra'].extend([float(ra)])
    d2['dec'].extend([float(de)])
    d2['redshift'].extend([table['z'][index]])
    d2['sigma_inst'].extend([worst_sigma])

df2 = pd.DataFrame(data=d2)
fmt2 = '%-15s', '%-1.6f', '%-1.6f', '%-1.6f', '%-1.2f', 
np.savetxt('list_sigma_clusters.txt', df2.values, delimiter='\t\t', header= 'name              RA             Dec           redshift     sigma_inst', fmt=fmt2)
    
##plt.xlim(3800,7000)
##plt.ylim(1000,3000)
#plt.minorticks_on()
#plt.xlabel(r'Wavelength ($\rm \AA$)')
#plt.ylabel(r'Instrumental sigma')
##plt.ylabel(r'Resolving power')
#plt.tight_layout()
##plt.savefig('res_power_sdss.png')
#plt.show()


