import pandas as pd
import numpy as np
import astropy.io.ascii as astro_ascii
from pipe_sdss import pipe_sdss
import matplotlib.pyplot as plt

'''Try processing voids data'''

if __name__ == '__main__':
    header_path = 'data/voids_header.csv'
    ext_path = 'data/voids_extinction.txt'
    gal_type = 'voids' # void, clusters or filaments
    w_lims = [3800, 7100]
    RV = 3.1

    # reading data and ext_table
    header_df = pd.read_csv(header_path)
    ext_table = astro_ascii.read(ext_path)

    if gal_type == 'voids':
        ...
    else:
        gal_id = header_df['specobjid'].iloc[i]
        filepath = f'data/{gal_type}/SDSS_spec{gal_id}.fits'
    for i in header_df.index:
        # voids_extinction and voids_header have same ordering of galaxies (indexed)
        gal_id = header_df['galaxy'].iloc[i]
        z = header_df['z'].iloc[i]
        filepath = f'data/voids/CAVITY{gal_id}.fits'
        print('Processing galaxy', i, 'id', gal_id)

        EBV = ext_table['mean_E_B_V_SandF'][i]
        # There is another option, exact position, but without std
        # AV = ext_table['AV_SandF'][i]
        stdEBV = ext_table['stdev_E_B_V_SandF'][i]
        AV = EBV*RV
        eAV = stdEBV*RV
        print('AV', AV)

        sdss = pipe_sdss(filepath, AV, eAV, w_lims, interpolate=True)

        ''' Correct for galactic extinction, degrade the specturm to uniform sigma and filter for wavelenghts
        given by wave_lims'''
        sdss.corr_ext()
        sdss.degrade()
        sdss.filter_wave()

        indx = np.arange(1, len(sdss.wave_air) + 1, 1) # First index necessary for PyPipe3D spectrum formatting
        print('sigma_inst', sdss.max_sigma_wave)
        np.savetxt(f'processed/{gal_type}/{gal_id}.txt', np.array([indx,sdss.wave_air,sdss.flux, sdss.eflux]).T)
        np.savetxt(f'processed/{gal_type}/{gal_id}_data.txt', [sdss.max_sigma_wave, z])
