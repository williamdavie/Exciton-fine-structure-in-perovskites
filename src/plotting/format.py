'''
File to store figure format.
'''

import matplotlib.pyplot as plt

colourPalette = ['#F24405','#FA7F08','#9EF8EE','#22BABB','#348888']

def runFormat() -> None:
    
    plt.rcParams['figure.constrained_layout.use'] = True
    plt.rcParams["font.family"] = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rc('axes', titlesize=20)     
    plt.rc('axes', labelsize=20)    
    plt.rc('xtick', labelsize=14)    
    plt.rc('ytick', labelsize=14)    
    plt.rc('legend', fontsize=18)
    plt.rc('figure', titlesize=16)