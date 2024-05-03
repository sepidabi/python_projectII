# import modules
import sparsetools as sp
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
import numpy as np
from sepid import *
import scipy as sci
from Tkinter import *
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

# Font properties
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }

# Directories
datadir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
outdir = '/home/seki2695/OUTPUT/project2/'
invdir = '/home/seki2695/INV/stic/II/slabs_tseries/'

# Difinitions
res = 0.0375 # CHROMIS pixel size in arcse
cad = 8
fr_select = [29,75,135,167] # frames to display
w_pos = 13 # wavelength position to display (K2V, K3, K2R) = (10,13,16)


cut_file = restore(datadir+'crispex_3950_2018-07-22T08:23:58_scans=0-212_time-corrected_2021Sep08_120534.csav')
#cut_file = restore(datadir+'crispex_3950_2018-07-22T08:23:58_scans=0-212_time-corrected_2021Sep08_120700.csav')

cut = np.transpose(cut_file.loop_slab[w_pos, :,:])

plt.close('all')

fig = plt.figure(figsize = (8,5.5))

plt.imshow(cut, cmap = 'gray', origin = 'lower')

plt.tight_layout()
plt.show()
