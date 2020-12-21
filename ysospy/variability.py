def calculateFluxAsymmetry(x):
    #Description: Returns the flux asymmetry metric M from Bredall et al. 2019
    
    #Import(s)
    import numpy as np
    from scipy import stats

    #Action
    x = list(x)
    
    bottom10 = [] #bottom 10% of magnitude measurements
    top10 = [] #top 10% of magnitude measurements

    for item in x:
        if stats.percentileofscore(x,item) <= 10:
            bottom10.append(item)
        elif stats.percentileofscore(x,item) >= 90:
            top10.append(item)
    
    topandbottom10 = np.array(bottom10+top10)
    m10 = np.mean(topandbottom10)
    x_median = np.median(x) 
    sd = np.std(x)
    M = (m10-x_median)/(sd) #Flux Asymmetry metric, M
    
    return M

def lombScargle(id,x,y,min_period,max_period,out_type):
    #Import(s)
    from astropy.timeseries import LombScargle
    import matplotlib.pyplot as plt
    from astropy.io import ascii
    import numpy as np 

    #Action
    dates = np.array(x)
    mags = np.array(y)
    #Calculate periodogram
    frequencies, powers = LombScargle(dates,mags).autopower(method='fastchi2',minimum_frequency=(1/max_period),maximum_frequency=(1/min_period))
    #Some statistics
    false_alarm_levels = list(LombScargle(dates,mags).false_alarm_level([0.1,0.05,0.01]))
    #Phase the data
    best_frequency = frequencies[np.argmax(powers)]
    T = 1/(float(best_frequency))
    phased_dates = np.mod(dates, T)/T #grabbed this from feets package documentation
    phased_dates_cycle_2 = phased_dates+1

    #Generate plots
    ax1 = plt.subplot(222)
    periods=1/frequencies
    ax1.plot(periods,powers)
    ax1.hlines(y=(false_alarm_levels[0]),xmin=min_period,xmax=max_period,colors=['silver'],linestyles='dashed',label='90% FAP')
    ax1.hlines(y=(false_alarm_levels[1]),xmin=min_period,xmax=max_period,colors=['darkgrey'],linestyles='dashed',label='95% FAP')
    ax1.hlines(y=(false_alarm_levels[2]),xmin=min_period,xmax=max_period,colors=['gray'],linestyles='dashed',label='99% FAP')
    ax1.set_xscale('log')
    ax1.set_xlabel('Period d')
    ax1.set_xscale('log')
    ax1.set_ylabel('Power')
    ax1.set_yticks([0.1,0.3,0.5,0.7,0.9])
    ax1.set_title('Periodogram',fontsize=10)
    ax1.legend(prop={'size':5})

    ax2 = plt.subplot(221)
    ax2.set_xlabel('Phase')
    ax2.set_ylabel('Mag')
    ax2.scatter(phased_dates,mags,s=2)
    ax2.scatter(phased_dates_cycle_2,mags,s=2,c='C0')
    ax2.invert_yaxis()
    ax2.set_title("Folded Light Curve, P: "+str(round((1/best_frequency),3))+' d',fontsize=10)
    ax2.locator_params(axis='y',nbins=5)
    
    ax3 = plt.subplot(212)
    ax3.scatter(dates,mags,s=2)
    ax3.invert_yaxis()
    ax3.set_ylabel('Mag')
    ax3.set_xlabel('HJD')
    ax3.set_title('Light Curve',fontsize=10)
    ax3.locator_params(axis='y',nbins=5)
    

    plt.subplots_adjust(wspace=0.3,hspace=0.45)

    if out_type == 'show':
        plt.show()
        plt.clf()
    else: 
        plt.savefig(out_type,dpi=200)
        plt.clf()
        
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
