def plotLightCurve(x,y,colors,x_label,y_label,plot_title,line_labels,plot_type,out_type,error_arrays):
    #Import(s)
    import matplotlib.pyplot as plt 
    import numpy as np
    import seaborn as sns
    
    #Action
    sns.set_style('darkgrid')
    plt.rcParams['font.family'] = 'Nimbus Roman'
    plt.title(plot_title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    

    item_index = 0
    for item in x:
        if plot_type == 'scatter':
            plt.scatter(item,y[item_index],s=2,c=(colors[item_index]),label=line_labels[item_index])
        elif plot_type == 'plot':
            plt.plot(item,y[item_index],markersize=2,linewidth=1,color=colors[item_index],label=line_labels[item_index])
        elif plot_type == 'Scatter_error':
            plt.error(item,y[item_index],yerr=error_arrays[item_index],marker='o',color=colors[item_index],ms=2,linewidth=0,label=line_labels[item_index])
        elif plot_type == 'plot_error':
            plt.error(item,y[item_index],yerr=error_arrays[item_index],marker='o',color=colors[item_index],ms=2,linewidth=0,label=line_labels[item_index])
        item_index = item_index+1
    
    plt.gca().invert_yaxis()
    plt.legend()
    
    if out_type == 'show':
        plt.show()
        plt.clf()
    else: 
        plt.savefig(out_type)
        plt.clf()
