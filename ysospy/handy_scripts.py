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
    
def removeInterval(x,y,interval):
    #Import(s)
    import numpy as np 

    #Action
    x = list(x)

    new_x = []
    indices = []
    new_ys = []
    lower_bound = float(interval.split(':')[0])
    upper_bound = float(interval.split(':')[1])    

    for item in x:
        if item > upper_bound or item < lower_bound:
            new_x.append(item)
            indices.append(x.index(item))

    for item in y:
        item = list(item)
        new_y = []
        for elem in indices:
            new_y.append(item[elem])
        new_ys.append(np.array(new_y))

    return np.array(new_x),new_ys

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

def returnSI(data):
    n = 0
    tempTup = []
    data = list(data)
    for elem in data:
        tempTup.append((elem, n))
        n += 1
    tempTup.sort()
    new_indices = []
    for elem in tempTup:
        new_indices.append(elem[1])
    return new_indices

def sortFI(data, indices): #indices must be less than or equal to data
    data=list(data)
    indices=list(indices)
    rearranged = []
    for elem in indices:
        rearranged.append(data[elem])
    return rearranged   
