import numpy as np
import matplotlib.pyplot as plt

def plotxy(x, y, aspect = 'auto', xlabel = '$t$', ylabel = '$f(t)$', ylim = (None, None)):
    
    fig = plt.figure(figsize=(7,7))
    
    ax = plt.subplot(111)
    ax.set_aspect(aspect)
    ax.set(xlabel = xlabel, ylabel = ylabel)
    
#    ax.set_ylim(ylim[0],ylim[1])
    ax.plot(x, y)
    
    
    
    plt.show()



def readXYData(fileName, delimiter = '\t', skip_header = 1, usecols = None, dtype = None):
    
    data = np.genfromtxt(fileName, delimiter = delimiter, skip_header = skip_header, usecols = usecols, dtype = dtype)
    
    return data  
    
def plotxyFile(fileName, delimiter = '\t', skip_header = 1, usecols = None, dtype = None, aspect = 'auto', xlabel = '$t$', ylabel = '$f(t)$', ylim = (None, None)):
    
    data = readXYData(fileName, delimiter = delimiter, skip_header = skip_header, usecols = usecols, dtype = dtype)
    
    X = data[0:, 0]
    Y = data[0:, 1]
    
    plotxy(X, Y, aspect = aspect, xlabel = xlabel, ylabel = ylabel, ylim = ylim)
    
    
    
