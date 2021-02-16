def codyM(data):
    #Import(s)
    from scipy import stats 
    import numpy as np
    
    #Action
    data = np.array(data)
    m_metric = (np.mean([stats.mstats.mquantiles(data,prob=0.9),stats.mstats.mquantiles(data,prob=0.1)])-np.median(data))/np.sqrt(((data-data.mean())**2).sum()/len(data))
    return m_metric

def lombScargle(x, y, yerr, fap_levels, min_period, max_period, out_type):
    # Import(s)
    from astropy.timeseries import LombScargle
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib

    # Action
    dates = np.array(x)
    mags = np.array(y)
    # Calculate periodogram
    frequencies, powers = LombScargle(dates, mags, yerr).autopower(method='fastchi2', minimum_frequency=(1 / max_period),
                                                             maximum_frequency=(1 / min_period))
    # Some statistics
    false_alarm_levels = list(LombScargle(dates, mags, yerr).false_alarm_level(fap_levels))
    # Phase the data
    best_frequency = frequencies[np.argmax(powers)]
    best_power = powers[np.argmax(powers)]
    T = 1 / (float(best_frequency))
    phased_dates = np.mod(dates, T) / T  # grabbed this from feets package documentation
    phased_dates_cycle_2 = phased_dates + 1

    # Generate plots
    plt.rcParams['font.family'] = 'serif'
    ax1 = plt.subplot(222)
    periods = 1 / frequencies
    ax1.plot(periods, powers)
    
    #Plot FAP Levels
    colors = ['lightgrey','silver','darkgray','gray','dimgray']

    for i in range(len(fap_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
        confidence_label = str(100*(1-fap_levels[i]))+'% FAP'
        ax1.hlines(y=(false_alarm_levels[i]), xmin=min_period, xmax=max_period, color = colors[i],linestyles='--',
               label=confidence_label)

    ax1.set_xscale('log')
    ax1.set_xlabel('Period d')
    ax1.set_xscale('log')
    ax1.set_ylabel('Power')
    #ax1.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax1.set_title('Periodogram', fontsize=10)
    box = ax1.get_position()
    ax1.set_position([box.x0,box.y0,box.width*0.5,box.height])
    ax1.legend(bbox_to_anchor=(1.15,0.5),loc='center',fontsize=4)

    ax2 = plt.subplot(221)
    xlabel = 'Phase (P = '+str(round((1 / best_frequency), 3))+' d)'
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel('Mag')
    ax2.scatter(phased_dates, mags, s=2)
    ax2.scatter(phased_dates_cycle_2, mags, s=2, c='C0')
    ax2.invert_yaxis()
    ax2.set_title("Folded Light Curve", fontsize=10)
    ax2.locator_params(axis='y', nbins=5)

    ax3 = plt.subplot(212)
    ax3.errorbar(dates, mags, yerr=yerr, lw=0,elinewidth=0.5)
    ax3.scatter(dates,mags,s=2)
    ax3.invert_yaxis()
    ax3.set_ylabel('Mag')
    ax3.set_xlabel('MJD')
    ax3.set_title('Light Curve', fontsize=10)
    ax3.locator_params(axis='y', nbins=5)

    plt.subplots_adjust(wspace=0.4, hspace=0.45)

    if out_type == 'show':
        plt.show()
        plt.clf()
    else:
        plt.savefig(out_type, dpi=200, format='svg')
        plt.clf()
    
    out_list = [(1/best_frequency),best_power]
    for fap in false_alarm_levels:
        out_list.append(float(fap))

    return out_list
        
def sergisonDistribution(x,percentiles):
    #Description: Returns variability amplitude metric and normalized magnitudes array as in Sergison et al. 2019 
    
    #Import(s)
    import numpy as np
    
    #Action
    mag_list = list(x)
    mean_mag = np.mean(np.array(mag_list))
    lower_percentile,upper_percentile = np.percentile(np.array(mag_list),[percentiles[0], percentiles[1]])
    amp_metric = lower_percentile-upper_percentile
    

    normalized_mags = []

    for item in mag_list:
        normalized_mag = (item-mean_mag)/amp_metric
        normalized_mags.append(normalized_mag)

    return amp_metric, np.array(normalized_mags) # Return results

def sokolovskyNu(x,xerr):
    #Import(s)
    import numpy as np

    #Action
    x = np.array(x)
    xerr = np.array(xerr)
    
    x_minus_xerr = x-xerr
    x_plus_xerr = x+xerr
    
    min_val = np.amin(x_minus_xerr)
    max_val = np.amax(x_plus_xerr)
    
    nu = (min_val-max_val)/(min_val+max_val)
    return nu
