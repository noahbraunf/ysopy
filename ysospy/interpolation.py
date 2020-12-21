def returnGoodRegions(x,y,max_sep,min_card):
    #Import(s)
    import numpy as np

    #Action
    dates_list = list(x)
    dates_array = np.array(x)
    differences = list(np.diff(dates_array))
    break_indices = []

    for item in differences:
        if item >= max_sep:
            break_indices.append(differences.index(item)+1)
    break_indices = np.array(break_indices)

    fx = np.split(ary=dates_array,indices_or_sections=break_indices)
    fy = np.split(ary=np.array(y),indices_or_sections=break_indices)
   
    #This section gets rid of intervals that don't meat the min_card limit
    i = 0
    arrays_to_remove = []
    for item in fx:
        if len(item) < min_card:
            arrays_to_remove.append(i)
        i += 1    
    
    for i in sorted(arrays_to_remove, reverse=True): #this is legendary
        del fx[i]
        del fy[i]

    return len(fx), fx, fy
