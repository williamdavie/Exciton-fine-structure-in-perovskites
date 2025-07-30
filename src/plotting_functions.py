import matplotlib.pyplot as plt
import matplotlib.axes as Axes
import numpy as np 

#--------------------------------------Formatting--------------------------------------

def fetchColour(colourname: str,) -> str:
    '''
    Fetch colour from main report colour palette. 
    
    options: 'blue','lightred','lime','red','lilac','purple','yellow','darkblue'
    '''

    colournames = ['green','blue','lightred','lime','red','lilac','purple','yellow','darkblue']
    
    # Bands = blue
    
    colourHex = ['#BEE4A8','#00BFC4','#F7766D','#3FF772','#F7473B','#8F42F4','#9302B8','#FDE725','#2A788E']
    
    return colourHex[colournames.index(colourname)]
    
    
def runFormat() -> None:
    
    plt.rcParams['figure.constrained_layout.use'] = True
    plt.rcParams["font.family"] = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['axes.linewidth'] = 1
    plt.rc('axes', titlesize=20)     
    plt.rc('axes', labelsize=22)    
    plt.rc('xtick', labelsize=16)    
    plt.rc('ytick', labelsize=16)    
    plt.rc('legend', fontsize=20)
    plt.rc('figure', titlesize=16)
    
    
#-----------------------------------Angle Landscape------------------------------------


def imshowLandscape(ax: Axes, values: np.ndarray, minXY: float, maxXY: float, Nsteps: int,
                       title: str='',colorbar: bool=True,text: bool=False):
    '''
    Given a set of values for a x-y landscape spanned by min-max, add an imshow plot to the axis 
    '''

    assert values.shape == (Nsteps,Nsteps)
    
    angles = np.linspace(minXY,maxXY,Nsteps)
    print(angles)
    
    labels = []
    for i in range(angles.shape[0]):
        if i % 2 == 0:
            labels.append(f'{angles[i]}')
        else:
            labels.append('')
            
    ax.set(xlabel=r'$\beta$ ($^o$)',ylabel=r'$\delta$ ($^o$)')
    ax.set_title(title,pad=10)
    ax.set_xticks(np.arange(0,Nsteps))
    ax.set_yticks(np.arange(0,Nsteps))
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    imshow = ax.imshow(values, origin='lower',cmap='plasma')
    
    if text:
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
        
                if i <= j:
                    text = f'{values[i, j]:.2f}'
                    ax.text(j, i, text, ha='center', va='center', color='black')
    
    if colorbar:
        # Use the axis' figure to add the colorbar
        cbar = ax.figure.colorbar(imshow, ax=ax, pad=0.04)
        
        
#-----------------------------------Band Structure------------------------------------

def fetchHighSymLabels(kpoints: np.ndarray, highsymmetrypoints: dict):
    

    highsymKeys = list(highsymmetrypoints.keys())
    label_indices = [0]
    labels = [highsymKeys[0]]
    for i in range(len(kpoints)):
            
        for j in range(len(highsymKeys)):
                
            a,b,c = highsymmetrypoints[highsymKeys[j]]
            a2,b2,c2 =  kpoints[i]
                
            if (a, b, c) == (a2, b2, c2):
                    
                label_indices.append(i)
                labels.append(highsymKeys[j])
                
    return label_indices, labels
            

def addBandStructure(ax: Axes, kpoints: np.ndarray, bandData: np.ndarray, efermi: float, startband: int=0, endband: int=0,
                     col: str='blue',alpha: float=1,label: str='',linestyle: str='-',linewidth: float=1.5):
    '''
    adds a band structure to an predifined ax.
    '''
    
    numKpoints, numBands = bandData.shape
    
    if endband == 0:
        endband = numBands

    labeled = False
    for band in range(startband,endband):
        
        if labeled == False:
            ax.plot(bandData[:,band]-efermi,c=col,alpha=alpha,label=label,linestyle=linestyle,linewidth=linewidth)
            labeled = True
        else:
            ax.plot(bandData[:,band]-efermi,c=col,alpha=alpha,linestyle=linestyle,linewidth=linewidth)
    
    return ax