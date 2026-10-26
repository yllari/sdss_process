SDSSS data processing for pyPipe3D

Run: 
```
   python main.py

   Change:
    # Comment or uncomment depending on what you want to do: correct extinction, degrade spectrum and filter wavelengths
    
    sdss.corr_ext()
    #sdss.degrade()
    sdss.filter_wave()
```

Needs:
```
    gal_type = 'clusters' # voids, clusters or filaments
    gal_id = 'specobjid' # identifier

    ext_path = f'data/{gal_type}_extinction.txt' # extinction file from IPAC with mean_E_B_V_SandF, std_E_B_V_SandF
    header_path = f'data/{gal_type}_header.csv' # file with redshift, z, and list of galaxies with same order as in ext_path
    in_path = 'data' # folder wher edata is 
    w_lims = [3800, 7100] # wavelength limits
    out_path = 'processed_mastar' # folder where data is out 
    RV = 3.1 # R_V for extinction law, Cardelli
    ncores = 4 # Number of cores
