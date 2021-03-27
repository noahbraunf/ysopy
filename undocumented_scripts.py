def q_diagnostic(id,mjd,mag,magerr,dv,min_per,max_per,false_alarm_levels,out_type):
    ''' 
    Create diagnostic plot for lomb-scargle periodogram as well as q-analysis. 
    Incorporates the best working version of our period search and quasi-periodicty
    routines. 
    '''
    # Import(s)
    from astropy.timeseries import LombScargle
    from pandas import Series, concat, date_range, to_datetime
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator

    #Action

    JD = mjd + 2400000.5

    # Generate a Dirac Comb, our window function
    time = to_datetime(JD, unit="D", origin="julian")
    time_not_obs = date_range(time.min(), time.max(), periods=1000)
    base = Series(np.zeros(len(time_not_obs)), index=time_not_obs)
    teeth = Series(np.ones(len(time)), index=time)
    dirac_comb = concat([base, teeth]).sort_index()
    
    minf = 1/max_per
    maxf = 1/min_per

    # Periodogram of the window function
    JD_W = dirac_comb.index.to_julian_date()
    mag_W = dirac_comb.values
    periodogram_W = LombScargle(JD_W, mag_W)
    freq_W, power_W = periodogram_W.autopower(method='fastchi2',minimum_frequency=minf, maximum_frequency=maxf)

    # Periodogram of original light curve
    periodogram = LombScargle(JD, mag, magerr)
    ls_freqs, ls_powers = periodogram.autopower(method='fastchi2',minimum_frequency=minf, maximum_frequency=maxf)

    # Mask out peak window-function frequencies from the data with a notch
    # width of dv (default should be 0.03) Hz on either side.
    high_power_W = power_W.mean() + 2 * power_W.std()
    
    pwff = freq_W[np.where(power_W>high_power_W)] # pwff = Peak Window Function Frequencies
    
    wffitr = np.array([]) # wfftr = Window Function Frequency Indices To Remove

    ls_pers = 1/ls_freqs

    for f in pwff:
        per = 1/f
        #print((per-0.02),per,(per+0.02))
        #good_idx = np.invert(np.logical_and((ls_freqs<(f+dv)), (ls_freqs>(f-dv)))) # Phew that is clean
        good_idx = np.invert(np.logical_and(ls_pers<per+dv,ls_pers>per-dv))
        bap = np.where(good_idx==False)
        wffitr = np.append(wffitr,bap)#.astype(int)
        wffitr = np.unique(wffitr).astype(int)    
    
    cleaned_powers = np.delete(ls_powers,wffitr)
    cleaned_freqs = np.delete(ls_freqs,wffitr)

    '''
    cleaned_powers = ls_powers[good_idx]
    cleaned_freqs = ls_freqs[good_idx]
    '''

    # Calculate FAPs
    faps = periodogram.false_alarm_level(false_alarm_levels)
    
    # Mask known aliase ranges (Lunar, etc.)
    cleaned_periods = 1/cleaned_freqs
    mask_lunar = np.invert(np.logical_and(cleaned_periods>26, cleaned_periods<30))
    cleaned_freqs = np.array(cleaned_freqs)[mask_lunar]
    cleaned_powers = np.array(cleaned_powers)[mask_lunar]

    # Find best results
    best_index = np.argmax(cleaned_powers)
    best_power = cleaned_powers[best_index]
    best_freq = cleaned_freqs[best_index]
    best_per = 1/best_freq

    # Fold the light curve
    T = 1 / (float(best_freq))
    phased_dates = np.mod(mjd, T) / T  # Got this from the feets package documentation
    phased_dates_cycle_2 = phased_dates + 1

    # Calculate Q and get residuals plot
    qp_results = quas_per(mjd=mjd,mag=mag,magerr=magerr,per=best_per)
    q = qp_results[0]
    residuals = qp_results[1]

    # Calculate m
    m = codyM(x=mag)

    def create_plot():
        # Plot
        plt.rcParams['font.family'] = 'serif'
        fig, axs = plt.subplots(2,2)
        fig.suptitle('Q: '+str(round(q,3))+'; M: '+str(round(m,3)),fontsize='medium')

        # Periodogram plot
        periods = 1 / cleaned_freqs
        axs[1,1].plot(periods[np.argwhere(periods<26)], cleaned_powers[np.argwhere(periods<26)])
        axs[1,1].plot(periods[np.argwhere(periods>30)], cleaned_powers[np.argwhere(periods>30)],color='C0')
        
        #Plot FAP Levels
        colors = ['lightgrey','silver','darkgray','gray','dimgray']

        for i in range(len(false_alarm_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
            confidence_label = str(100*(1-false_alarm_levels[i]))+'% FAP'
            axs[1,1].hlines(y=(faps[i]), xmin=min_per, xmax=max_per, color = colors[i],linestyles='--',
                label=confidence_label)

        axs[1,1].set_xscale('log')
        axs[1,1].set_xlabel('Period d')
        axs[1,1].set_xscale('log')
        axs[1,1].set_ylabel('Power')
        axs[1,1].set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
        axs[1,1].yaxis.set_minor_locator(AutoMinorLocator())
        axs[1,1].set_ylim(0,1)
        axs[1,1].set_xlim(0.5,250)
        axs[1,1].set_title('Periodogram', fontsize=10)
        box = axs[1,1].get_position()
        axs[1,1].set_position([box.x0,box.y0,box.width*0.5,box.height])
        axs[1,1].legend(bbox_to_anchor=(1.15,0.5),loc='center',fontsize=4)

        # Folded light curve plot
        xlabel = 'Phase (P = '+str(round((1 / best_freq), 3))+' d)'
        axs[0,1].set_xlabel(xlabel)
        axs[0,1].set_ylabel('Mag')
        axs[0,1].scatter(phased_dates, mag, s=2)
        axs[0,1].scatter(phased_dates_cycle_2, mag, s=2, c='C0')
        axs[0,1].xaxis.set_minor_locator(AutoMinorLocator())
        axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())
        axs[0,1].invert_yaxis()
        axs[0,1].set_title("Folded Light Curve", fontsize=10)
        axs[0,1].locator_params(axis='y', nbins=5)

        # Unfolded light curve plot
        axs[0,0].errorbar(mjd, mag, yerr=magerr, lw=0,elinewidth=0.5)
        axs[0,0].scatter(mjd,mag,s=2)
        axs[0,0].invert_yaxis()
        axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
        axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())
        axs[0,0].set_ylabel('Mag')
        axs[0,0].set_xlabel('MJD')
        axs[0,0].set_title('Light Curve', fontsize=10)
        axs[0,0].locator_params(axis='y', nbins=5)

        # Residuals plot
        axs[1,0].scatter(mjd,residuals,s=2)
        axs[1,0].axhline(y=0,xmin=0,xmax=1,color='black',lw=0.75)
        axs[1,0].set_title('Residual Plot')
        axs[1,0].xaxis.set_minor_locator(AutoMinorLocator())
        axs[1,0].yaxis.set_minor_locator(AutoMinorLocator())
        axs[1,0].set_xlabel('MJD')
        axs[1,0].set_ylabel('Residual Mag',fontsize=10)

        plt.subplots_adjust(wspace=0.4, hspace=0.50)

    if out_type == 'show' or out_type == 'Show':
        create_plot()
        plt.show()
        plt.clf()
    
    elif out_type == 'None' or out_type == 'none':
        plt.clf()
    
    else:
        create_plot()
        q_str = str(q) 
        q_str_split = q_str.split('.')
        decimal_portion = q_str_split[1]
        decimal_portion = decimal_portion[0:3]
        q_str = q_str_split[0]+'.'+decimal_portion

        actual_path = out_type.replace('***',q_str)
        actual_path = actual_path.replace('+++',id)
        plt.savefig(actual_path, dpi=200, format=out_type[-3:]) # Flexible save type (svg, png, etc.)
        plt.clf()
        plt.close()
    
    out_list = [best_power, 1/best_freq,q,m]
    for fap in faps: 
        out_list.append(fap)
    
    return out_list
 
def fig_plot_2(data_file):
    '''
    Create figure two plot (from Bredall et al. 2020) from data file which has all the IDs, Qs, and Ms 
    '''
    
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.ticker import MultipleLocator

    # Action

    df = pd.read_csv(data_file)
    ids = list(df['ID'])
    q_raw = np.array(df['Q'])
    q_mask = np.logical_and(q_raw>0,q_raw<1)
    m = np.array(df['M'])[q_mask]
    q_raw = q_raw[q_mask]
    
    # Transform Q
    q = (q_raw-min(q_raw))*(1/(max(q_raw)-min(q_raw)))
    
    # Plot
    plt.rcParams['font.family'] = 'serif'

    plt.scatter(q,m,marker='.',color='orange',s=2) #slategray

    plt.figtext(0.18,0.9,'Periodic',ha='center',fontsize='small')
    plt.figtext(0.5,0.9,'Quasi-Periodic',ha='center',fontsize='small')
    plt.figtext(0.845,0.9,'Aperiodic',ha='center',fontsize='small')

    plt.figtext(0.91,0.685,'Bursting',rotation=270,fontsize='small')
    plt.figtext(0.91,0.43,'Symmetric',rotation=270,fontsize='small')
    plt.figtext(0.91,0.22,'Dipping',rotation=270,fontsize='small')

    plt.axhline(y=0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axhline(y=-0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=0.15,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=0.85,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.xlim(0,1)
    plt.ylim(-1,1)
    plt.gca().invert_yaxis()
    plt.xlabel('Quasi-Periodicity (Q)')
    plt.ylabel('Flux Asymmetry (M)')
    plt.axes().yaxis.set_minor_locator(MultipleLocator(0.05))
    plt.axes().xaxis.set_minor_locator(MultipleLocator(0.05))
    plt.show()
