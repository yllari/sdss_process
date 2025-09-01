import numpy as np
import requests
import pandas as pd
import xml.etree.ElementTree as ET
import time
import extinction


''' In this module there is all functions needed for extinction-related purposes in the SDSS preparation
pipeline. NOTE: TO OBTAIN EXTINCTIONS IT IS BETTER TO USE https://irsa.ipac.caltech.edu/applications/DUST/ and
upload an IPAC table with all RA and DEC. This is much faster than single http queries and gives A_V calculated.

To convert from CSV to IPAC use utils.csv_to_ipac
'''

def single_redd(objid : str, dataheader):
    ''' objid has to be CAVITY void galaxy or specobjid
    :param objid: CAVITY assigned number for void galaxy or SDSS specobjid
    :return: value of reddening E(B-V) as returned from the IRSA IPAC DUST application
    '''

    header_info  = pd.read_csv(dataheader)
    if 'void' in str(dataheader):
        header_info = header_info[header_info['galaxy'] == objid]
        if len(header_info.index) != 1:
            raise Exception (f'Error retriving info, expected 1 row with objid but instead got {len(header_info.index)}')
        else:
            ra = header_info['ra'].iloc[0]
            dec = header_info['dec'].iloc[0]
    else:
        header_info = header_info[header_info['specobjid'] == objid]
        if len(header_info.index) != 1:
            raise Exception (f'Error retriving info, expected 1 row with objid but instead got {len(header_info.index)}')
        else:
            ra = header_info['ra'].iloc[0]
            dec = header_info['dec'].iloc[0]

    # Now we prepare the query

    # Info about the fetch process https://irsa.ipac.caltech.edu/applications/DUST/docs/dustProgramInterface.html
    query = f"http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={str(ra)}+{str(dec)}"

    t_0 = time.time()
    res_xml = requests.get(query + '+equ+j2000').text
    t_1 = time.time()
    print(res_xml)
    print('query time', t_1 - t_0)
    root = ET.fromstring(res_xml)

    # We are extracting SandF info
    # We trim to avoid the '(mag)' string.
    ''' WARNING: THIS IS REDDENING, NOT EXTINCTION '''
    E_BV = float(next(root.iter('meanValueSandF')).text.rstrip('\n').split()[0])
    sigma_EBV = float(next(root.iter('stdSandF')).text.rstrip('\n').split()[0])

    return E_BV, sigma_EBV

def correct_extinction(wave, flux, eflux, Av, eAv, Rv=3.1):
    ''' Given a spectra, correct for extinction saved in another file
    and return
    :param wave: wavelength for the measured flux
    :param flux: value fo flux in weach wavelength
    :param eflux: uncertainty associated to flux for each wavelength
    :param a_v: mean extinction in V band
    :param r_v: Ratio of total to selective extinctio, A_V/E(B-V)
    :return: corrected flux and associated eflux
    '''


    ext_model = extinction.ccm89(wave.astype(float), Av, r_v=Rv)
    flux_corr= extinction.remove(ext_model, flux)
    #eflux_corr = extinction.remove(ext_model, eflux)

    eflux_corr = np.sqrt((10**(0.4*ext_model)*eflux)**2 + (flux_corr*0.4*np.log(10)*eAv*ext_model/Av)**2)

    return flux_corr, eflux_corr


if __name__ == '__main__':
    t_0 = time.time()
    print(*single_redd(54708, 'data/voids_header.csv'))
    t_1 = time.time()
    print('Elapsed time', t_1-t_0)
