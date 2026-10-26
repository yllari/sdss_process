import pandas as pd
import numpy as np
import astropy.io.ascii as astro_ascii
from pipe_sdss import pipe_sdss
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib.pyplot as plt
import time
import os

'''Try processing voids data'''

def call_save_pipesdss(filepath, out_path, AV, eAV, w_lims, gal_id, z, interpolate=False):
    ''' Correct for galactic extinction, degrade the specturm to uniform sigma and filter for wavelenghts
    given by wave_lims'''
    sdss = pipe_sdss(filepath, AV, eAV, w_lims, interpolate=interpolate)

    # Comment or uncomment depending on what you want to do: correct extinction, degrade spectrum and filter wavelengths
    sdss.corr_ext()
    #sdss.degrade()
    sdss.filter_wave()

    indx = np.arange(1, len(sdss.wave_air) + 1, 1) # First index necessary for PyPipe3D spectrum formatting
    # Air wavelength is used here
    np.savetxt(f'{out_path}/{gal_type}/{gal_id}.txt', np.array([indx, sdss.wave_air, sdss.flux, sdss.eflux]).T)
    np.savetxt(f'{out_path}/{gal_type}/{gal_id}_data.txt', [sdss.max_sigma_wave, z])

    return None

def multiprocess_sdss(in_path, out_path, header_path, ext_path, gal_type, w_lims, RV, ncores, interpolate=False):

    os.makedirs(f'{out_path}/{gal_type}', exist_ok=True) 
    procs = []
    # reading data and ext_table
    header_df = pd.read_csv(header_path)
    ext_table = astro_ascii.read(ext_path)

    if gal_type == 'voids':
        primary_id = 'galaxy'
        file_root = 'CAVITY'
    else:
        primary_id = 'specobjid'
        file_root = 'SDSS_spec'
    with ProcessPoolExecutor(ncores) as executor:        
        for i in header_df.index:
            # voids_extinction and voids_header have same ordering of galaxies (indexed)
            gal_id = header_df[primary_id].iloc[i]
            z = header_df['z'].iloc[i]
            filepath = f'{in_path}/{gal_type}/{file_root}{gal_id}.fits'

            EBV = ext_table['mean_E_B_V_SandF'][i]
            # There is another option, exact position, but without std
            # AV = ext_table['AV_SandF'][i]
            stdEBV = ext_table['stdev_E_B_V_SandF'][i]
            AV = EBV*RV
            eAV = stdEBV*RV
            proc = executor.submit(call_save_pipesdss, filepath, out_path, AV, eAV, w_lims, gal_id, z,  interpolate=interpolate)
            procs.append(proc)
        total = len(procs)
        current = 0
        for proc in as_completed(procs):
            _ = proc.result()
            current += 1
            print(f'Progress: {current}/{total}')



def check_all_files(gal_type, primary_id):
    header_df = pd.read_csv(header_path)
    failed = []
    for gal_id in header_df[primary_id]:
        try:
            spectra = np.genfromtxt(f'processed/{gal_type}/{gal_id}.txt')
        except:
            failed.append(gal_id)
            print('Failed on', gal_id)
    print('No errors')
    return failed


if __name__ == '__main__':
    gal_type = 'clusters' # voids, clusters or filaments
    gal_id = 'specobjid'
    ext_path = f'data/{gal_type}_extinction.txt'
    header_path = f'data/{gal_type}_header.csv'
    in_path = 'data'
    w_lims = [3800, 7100]
    out_path = 'processed_mastar'
    RV = 3.1
    ncores = 4
    t_0 = time.time()
    multiprocess_sdss(in_path, out_path, header_path, ext_path, gal_type, w_lims, RV, ncores, interpolate=False)
    print('Elapsed time', time.time() - t_0)
    #print(check_all_files(gal_type, gal_id))


