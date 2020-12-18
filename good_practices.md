- By default, np.std(x) returns the population standard deviation of an array—to get the sample standard deviation (N-1), set the optional ```ddof``` equal to the integer 1. 
- Read data columns into arrays from ipac format .tbl files:
```
from astropy.io import ascii
import numpy as np
path = '\home\example_file.tbl'
r = ascii.read(path,format='ipac',delimiter='|')
red_dates = np.array(r['hjd'])
red_mags=np.array(r['mags'])
red_mag_errors=np.array(r['magerr'])
``` 
- Good routine for pre-treating ZTF data prior to analysis: 
  1. Remove known "clump(s)" from data with the ```removeInterval()``` function.
  2. Sort data with the ```sortData()``` function.
