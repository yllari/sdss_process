import numpy as np
import pandas as pd
from astropy.table import Table

def csv_to_ipac(filepath, columns, overwrite=True):
    ''' This is meant to trim and transform to IPAC to RA DEC, only few coolumns'''
    df = pd.read_csv(filepath)
    df = df[columns]

    col_str = ''.join("_" + str(col) for col in columns)

    table = Table.from_pandas(df)

    table.write(filepath.replace('.csv','') + '_ipac' + col_str + '.txt', format='ipac', overwrite=overwrite)

if __name__ == '__main__':
    csv_to_ipac('data/filaments_header.csv', ['ra', 'dec'])
    


