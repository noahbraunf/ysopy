def calculatePeakLocs(x,w_val):
    #Import(s)
    from scipy.signal import find_peaks
    import numpy as np
    
    #Action
    x = np.array(x)
    fakepeaks, _ = find_peaks(x,width=w_val)
    realpeaks, _ = find_peaks((x*-1),width=w_val)
    outarray = np.concatenate((fakepeaks,realpeaks))
    return(outarray)

def queryCoordSimbad(raw_coord,search_radius):
    #Import(s)
    import numpy as np
    from astropy import coordinates as coord
    from astropy import units as u
    from astroquery.simbad import Simbad
    from astropy.coordinates.sky_coordinate import SkyCoord

    #Action
    c = SkyCoord(raw_coord,unit=(u.hourangle,u.deg))
    c = c.to_string('hmsdms')
    result_table = Simbad.query_region(coord.SkyCoord(c, frame='icrs'),radius=('0d0m'+str(search_radius)+'s'))
    names_col = result_table['MAIN_ID']
    id = str(names_col[0])[1:]
    return id

def returnDistances(x):
    #Import(s)
    import numpy as np
    
    #Action
    input_array = np.array(x)
    out_array = np.diff(input_array)
    return out_array
    
def removeIntervals(x,y,intervals):
    #Import(s)
    import numpy as np 

    #Action
    x = list(x)
    y = list(y)
    
    for item in intervals:
        lower_bound = float(item.split(':')[0])
        upper_bound = float(item.split(':')[1])
        for elem in x:
            if elem < upper_bound:
                if elem > lower_bound:
                    elem_index = x.index(elem)
                    x.remove(elem)
                    y.remove(y[elem_index])
    
    return np.array(x),np.array(y)

def sortData(x,y):
    #Import(s)
    import numpy as np
    
    #Action
    
    unsorted_dates = list(x)
    
    sorted_dates = list(np.sort(unsorted_dates))
    sorted_arrays = []
    for item in y:
        unsorted_list = list(item)
        sorted_list = []
        for elem in sorted_dates:
            #newIndex = sorted_dates.index(elem)
            oldIndex = unsorted_dates.index(elem)
            sorted_list.append(unsorted_list[oldIndex])
        sorted_array = np.array(sorted_list)
        sorted_arrays.append(sorted_array)

    sorted_dates = np.array(sorted_dates)
    
    return sorted_dates, sorted_arrays
