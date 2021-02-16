# ysopy 
This is the documentation for all the files and functions in our ysopy package. **Note: documentation is not caught up with some of the code yet**

## Guidelines
#### Naming and Style Conventions 
- Module (file) names should be in snake_case, and functions in camelCase. 
- Example: 
```
from ysopy.plotting_funcs import plotLightCurve
```
#### Requirements 
- numpy
- matplotlib
- scipy
- sklearn
- seaborn
- progress
- astropy
- pandas
#### General notes
- All functions automatically import the modules they require.  
## handy_scripts.py
- **Overview:** this file contains miscellaneous helpful routines. 
### calculatePeakLocs()
- **Summary:** this function uses the ```scipy.find_peaks``` function to return the indices of all relative extrema in an array of magnitudes. 
- **Parameters:**
  - x: the array of magnitudes (which has been previously sorted to a corresponding array of dates)
  - w_val: the required minimal horizontal distance between neighbouring peaks (scipy). This function is different than the ```scipy.find_peaks()``` it's based on because it finds relative minimums as well as relative maximums (it does this by find peaks in the array multiplied by negative one). 
#### Example: 
```python
from ysopy.handy_scipts import calculatePeakLocs
peak_indices = calculatePeakLocks(x=example_array_of_mags,w_val=3)
```
### queryCoordSimbad()
- **Summary:** this function returns the string identifier from Simbad for a given input coordinate (it takes the closest result from the Simbad result table).
- **Parameters:**
  - raw_coord: the coordinate for the object in question should be in this form: ```raw_coord='20:50:32.32+44:26:17.4'```.
  - search_radius: the search radius (in arcseconds), should be an integer value. ```search_radius=5``` is the recommended default.
##### Example: 
```python
from ysopy.handy_scipts import queryCoordSimbad
obsid = queryCoordSimbad(raw_coord='20:54:24.41+44:48:17.3',search_radius=5)
```
### removeInterval()
- **Summary:** this function removes data that fall between a parameter defined interval of dates. ```removeInterval(...)[0]``` returns the modified dates array, and ```removeInterval(...)[1]``` returns a list of the modified input arrays, so to call the modified input arrays individually, use an additional index, e.g. ```removeInterval(...)[1][0]``` returns the first modified array that you provided in ```y```. 
- **Parameters:**
  - x: the sorted date array. 
  - y: a list of predefined array names (e.g. ```y=[mags,magerrs]```).   
  - intervals: this should be a string in which which should be formated like this: ```'lower_date_bound:upper_date_bound'```. All data in the provided open interval (lower_date_bound,upper_date_bound) will not be present in the return date array and the ```y``` arrays. For example, if one wanted to remove all observations between HJD 2456202 and HJD 2456205 in a data set, set ```interval``` to ```'2456202:2456205'```.
#### Example: 
```python
#Import(s)
from ysopy.handy_scripts import removeInterval
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

#Action
r = ascii.read(data_file,format='ipac',delimiter='|')
r_dates = r['mag']
r_mags=r['hjd']
r_magerrs = r['magerr']

f = removeInterval(x=r_dates,y=[r_mags,r_magerrs],interval='2458455.6:2458455.8')
r_dates2=f[0]
r_mags2=f[1][0] 
r_magerrs2=f[1][1]

f=removeInterval(x=r_dates2,y=[r_mags2,r_magerrs2],interval='2458437:2458438')
r_dates3=f[0]
r_mags3=f[1][0]
r_magerrs3 = f[1][1]

plt.errorbar(r_dates3,r_mags3,yerr=r_magerrs3,color='tomato',lw=0,marker='o',ms=2,elinewidth=1)
plt.ylabel('mag')
plt.xlabel('HJD')
plt.gca().invert_yaxis()
plt.show()
```
<img src="https://github.com/noahbraunf/ysopy/blob/main/images/no_clumps.png" width="350" height="230">
- Note how the clumps in the below zoomed-in image (between HJD 2458455.6-2458455.8 and HJD 2458437-2458438) are not present in the above plot. 
<img src="https://github.com/noahbraunf/ysopy/blob/main/images/clumps_zoomin.png" width="350" height="230">

### returnDistances()
- **Summary:** this function returns an array of all the distances between date values in an array. It can really be replaced by the one line np.diff(a) function that it's based off (it used to be a more extensive function, but then we realized we can just simplify it to np.diff). 
- **Parameters:** 
  - x: the array of dates in which you wish to find the differences between consecutive elements. 
#### Example:
```python
from ysopy.handy_scripts import returnDistances
differences=returnDistances(example_date_array)
```
### sortData()
- **Summary:** this function sorts an array of dates and then sorts its corresponding, provided array of mags to the same order (so corresponding dates and mags are matched). The function returns the sorted dates array first and then returns the sorted mags array (as usual call these individually by indexing the function's output). 
- **Parameters:**
  - x: the date array. 
  - y: the corresponding array of magnitudes. 
#### Example: 
```python
from ysopy.handy_scripts import sortData
sorted_dates=sortData(example_dates_array,example_mags_array)[0]
sorted_mags=sortDate(example_dates_array,example_mags_array)[1]
```
### returnSI()
- **Summary:** this function, when used with ```sortFI()```, provides an alternative way to sort complementary lists according to one sorted list. ```returnSI``` stands for ```returnSortedIndices``` because this function takes in an input list, sorts it, and returns the indices of the sorted elements.
- **Parameters:** 
  - data: the list/array of dates that you want to sort and return the sorted indices from (for sorting other arrays to that sorted dates list/array). 
### sortFI()
- **Summary:** this function sorts and returns a list/array based on a given list of indices to sort it to. 
- **Parameters:**
  - data: the list/array to sort based on the following provided list of indices 
  - indices: the list of indices to sort data on. 
### clean_clusters()
- **Summary:** this function removes the weird "clusters" in ZTF lightcurve data. The function does not remove clusters based on *a priori* knowledge of where the clusters occur in time. Rather, it bins the data into one day bins from ```min(dates)``` to ```max(dates)``` and calculates the mean "bin_count" as well as the standard deviation of the "bin_counts". Then it removes data that occur in bins which have more data points that the cutoff, which is calculated as ```np.mean(bin_counts)+tolerance*(np.std(bin_counts,ddof=1))```. ```tolerance``` is parameter defined number of standard deviations from the mean that the user decides on. For example, a tolerance set to ```1``` means that all bins with counts > mean counts + one standard deviation are deleted from the returned lists. 
- **Parameters:**
  - dates: the list of observation dates (MJD, HJD, JD, etc.)
  - paired_lists: a list of the lists containing data values paired to the date values in ```dates```, e.g. it may have lists of mags and magerrs. 
  - tolerance: integer/float for use in the cutoff calculation mentioned previously. 
#### Example: 
```python
from ysopy.handy_scripts import clean_clusters
r_dates = [] #Pretend this list is populated
r_mags = [] # ""
r_magerrs = [] # ""

cf = clean_clusters(dates=r_dates,paired_lists=[r_mags,r_magerrs],tolerance=3)

cleaned_r_dates = cf[0]
cleaned_r_mags = cf[1][0]
cleaned_r_magerrs = cf[1][1]
```
## interpolation.py
- **Overview:** this file contains routines that may be a great help when interpolating data.
### returnGoodIntervals()
- **Summary:** this function returns sub-interval arrays that meet certain criteria (maximum distance between elemental points, minimum cardinality of sub-interval) from an input array of dates. This function was originally written to improve interpolation quality by getting rid of regions where the data are sparsely sampled. There are three outputs of the function, which can be called by index. The first is the number of good sub-intervals returned, the second is a list of the good sub-interval arrays of dates, and the third is a list of the good sub-interval mags. 
- **Parameters:**
  - x: the sorted array of dates that you wish to split into good sub-intervals. 
  - y: the sorted array of magnitudes that correspond to the dates in x (as this function goes ahead and modifies the mag array while it's modifying the date array).
  - max_sep: the maximum separation in days allowed between elements in a good sub-interval.  
  - min_card: the minimum cardinality for good intervals. Intervals of data with lesser cardinality will be removed. 
#### Example: 
```python
#Import(s)
import matplotlib.pyplot as plt
from ysopy.interpolation import returnGoodIntervals

#Action
num_arrays = returnGoodRegions(x=sgd,y=sgm,max_sep=20,min_card=3)[0] 
x_arrays = returnGoodRegions(x=sgd,y=sgm,max_sep=20,min_card=3)[1]
y_arrays = returnGoodRegions(x=sgd,y=sgm,max_sep=20,min_card=3)[2]
i = 0
colors = ['red','blue','green','yellow','purple','black','brown','pink','orange']
for elem in x_arrays:
    plt.scatter(elem,y_arrays[i],s=4,c=colors[i])
    i += 1
plt.gca().invert_yaxis()
plt.xlabel('HJD')
plt.ylabel('Mag')
plt.show()
```
<img src="https://github.com/thissop/YSOs/blob/main/ysospy/images/mincardinaction.png" width="350" height="230">
- Note that in this example, sgd was a predefined array of dates, and that sgm was a predefined array of mags. Also note that isolated data present in the graph below have been removed in the plot of the returnGoodIntervals output above.
<img src="https://github.com/thissop/YSOs/blob/main/ysospy/images/betterplotofprev.png" width="350" height="230">

## plotting_funcs.py
- **Overview:** this file contains general plotting routines. 
### plotLightCurve()
- **Summary:** you can pass multiple date and mag arrays for plotting. You can also set whether the output is a plot opened in a new window or a saved file (because if we set it to open new plots for hundreds of files...forget *:(){ :|:& };:*, opening hundreds of plots would probably dos your computer just as easily xD).  
- **Parameters:** 
  - x: should be a list of all the names of the date arrays you want to plot. 
  - y: should be a list of all the names of the magnitude arrays you want to plot, with each array in the same indice position as its corresponding date array in x. 
  - colors: a list of colors for each of the lines you want to plot, e.g. colors=['green','red'] will set the first line green and the second red. 
  - x_label: the string label for the x-axis. 
  - y_label: the string label for the y-axis. 
  - plot_title: the string label for the plot title.
  - line_labels: this should be a list of string labels for each of the lines being plotted. They should be in the same list order as the items in x in y. 
  - plot_type: there are four plot_type options:
    1. 'scatter': simple scatter plot
    2. 'plot': simple line plot
    3. 'Scatter_error': scatter plot with error bars 
    4. 'plot_error': line plot with error bars
  - out_type: should either be 'show', which results in plt.show(), or a string filepath (which will save the plot as the given path).
  - error_arrays: the names of the arrays of error values. Set it to 'N/A' when plot_type is not set to either 'scatter_error' or 'plot_error'
#### Example
```python
from ysopy.plotting_funcs import plotLightCurve
plotLightCurve(x=[sgd,srd],y=[sgm,srm],colors=['green','red'],x_label='HJD',y_label='Mag',plot_title='Example Plot',line_labels=['Green Band','Red Band'],plot_type='scatter',out_type='show',error_arrays='N/A')
```
<img src="https://github.com/thissop/YSOs/blob/main/ysospy/images/example_plotLightCurve.png" width="350" height="230">

*Note that in this example, sgm and srm represent previously defined arrays of stellar magnitudes sorted by their corresponding date arrays, sgd (sorted green dates) and srd (sorted red dates)* 
## variability.py
- **Overview:** this file contains variability characterization routines. 
### codyM()
- **Summary**: this function calculates and returns the asymmetry metric M that first appeared in Cody et al. 2014. 
- **Parameters**
  - x: array of mags
### lomgScargle()
- **Summary:** this functions runs off Astropy's lomb-scargle class. Creates a plot that can be opened or saved and returns results. Note: we calculate periodograms as a function of frequency but for ease of interpretation, we display periodograms as a function of period. 
- **Parameters:**
  - x: array of observation dates
  - y: array of observation mags
  - yerr: array of magerrs 
  - fap_levels: list of false alarm probability levels to calculate. Example: if you wanted to calculate the 95% and 99% confidence false alarm probabilities, set ```fap_levels=[0.05,0.01]```.
  - min_period: minimum period to calculate (days)
  - max_period: maximum period to calculate (days)
  - out_type: if set to 'show' the plot is shown (i.e. ```plt.show()```), otherwise the function will save the plot to the path set to ```out_type```.
- **Returns:**
  - ```[0]```: best period, in days
  - ```[1]```: power of best period
  - ```[2-n]```: each calculated false alarm level in the same order as provided in ```fap_levels```. Example: if three separate false alarm probabilities were calculated, the function would return a total of five items (the last of which having the index ```[4]```, of course). 
#### Example: 
```python
import pandas as pd
from ysopy.variability import lombScargle
r_file = '/home/2MASS_J20470481+4349114_r.csv'
df = pd.read_csv(r_file)
r_mags = list(df['mag'])
r_mjds = list(df['mjd'])
r_magerrs = list(df['magerr'])

lombScargle(x=r_mjds,y=r_mags,yerr=r_magerrs,fap_levels=[0.1, 0.05, 0.01, 0.001],min_period=1.25,max_period=399.99,out_type='/home/Figure_1.svg')
``` 
### sergisonDistribution()
- **Summary:** this is a function for normalizing magnitudes as in Sergison et al. 2019 §4. It returns two outputs, which can be returned individually by calling  the first or second index of the function, which returns either the AH68 metric or an array of the normalized magnitudes, respectively. 
- **Parameters:**
  - x: an array or list of magnitudes that you wish to normalize and or find the AH68 metric for. 
  - percentiles: the AH68 metric is calculated with the 16th and 84th percentiles (as in the paper), but in case you want to calculate a similiar measure of the amplitude of variability in the array, you can enter the percentiles as integer items in a list and set that list equal to the percentiles argument, e.g. percentiles=[5,95] (the range that confines 90% of the data). 
#### Example: 
```python
from ysopy.plotting_funcs import sergisonDistribution()
import seaborn as sns
import matplotlib.pyplot as plt
norm_mags = sergisonDistribution(x=srm,percentiles=[16,84])[1] #percentiles for AH68 metric
sns.histplot(norm_mags,kde=True)
plt.xlabel('Normalized Magnitudes')
plt.ylabel('Counts')
plt.show()
```
<img src="https://github.com/thissop/YSOs/blob/main/ysospy/images/sergisonDist.png" width="350" height="230">
- Note that in this example, srm is a predefined magnitudes array. 

