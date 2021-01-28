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

def clean_clusters(dates,paired_lists,tolerance):
    #Import(s)
    import numpy as np
    import math

    #Action

    min_date = math.floor(np.min(dates))
    max_date = math.ceil(np.max(dates))
    date_bins = list(range(min_date,max_date,1))
    binned_dates = np.digitize(x=dates,bins=date_bins)
    bin_counts = np.bincount(binned_dates)

    nz_bin_counts = bin_counts[np.nonzero(bin_counts)] #Get rid of bin_counts = 0 for calculations

    mean_counts = np.mean(nz_bin_counts)
    std_counts = np.std(nz_bin_counts,ddof=1)

    cutoff = mean_counts+float(tolerance)*(std_counts)

    intervals_to_remove = []
    
    bin_counts = list(bin_counts)

    for count in bin_counts:
        if count > cutoff:
            bin_index = bin_counts.index(count)
            lower_edge = date_bins[bin_index-1]
            upper_edge = date_bins[bin_index]
            edge_string = str(lower_edge)+':'+str(upper_edge)
            intervals_to_remove.append(edge_string)

    indices_to_remove = []

    for interval in intervals_to_remove:
        interval_list = interval.split(':')
        lower = float(interval_list[0])
        upper = float(interval_list[1])
        for date in dates:
            if date > lower and date < upper:
                date_index = dates.index(date)
                indices_to_remove.append(date_index)

    for i in sorted(indices_to_remove, reverse=True):
        del dates[i]
        for paired_list in paired_lists:
            del paired_list[i]

    return dates, paired_lists
