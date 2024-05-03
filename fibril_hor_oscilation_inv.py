'''
fibril_hor_oscilation_inv.py

program to fast plot
the inversion results
of cuts across a fibril
'''

from sepid import *
import sparsetools as sp
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
import numpy as np
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

# getting the fibril and cuts information
fibdir = datadir+'fr29/'
fib_file = (file_search(fibdir,'crispex*3950*.csav'))[0]
fib = restore(fibdir+fib_file)
cut0_file = (file_search(fibdir,'crispex*3950*.csav'))[1]
cut1_file = (file_search(fibdir,'crispex*3950*.csav'))[2]
cut2_file = (file_search(fibdir,'crispex*3950*.csav'))[3]

xmin = int(np.min([np.min([cut1.x_coords, cut2.x_coords,
                           cut0.x_coords]),np.min(fib.x_coords)]) - 5)
ymin = int(np.min([np.min([cut1.y_coords, cut2.y_coords,
                           cut0.y_coords]),np.min(fib.y_coords)]) - 5)
xmax = int(np.max([np.max([cut1.x_coords, cut2.x_coords,
                           cut0.x_coords]),np.max(fib.x_coords)]) + 5)
ymax = int(np.max([np.max([cut1.y_coords, cut2.y_coords,
                           cut0.y_coords]),np.max(fib.y_coords)]) + 5)


# Reading the inv results
cut0 = sp.model('slab_nostokes_atmosout_cycle1_31_172118.nc')
cut1 = sp.model('slab_nostokes_atmosout_cycle1_31_172059.nc')
cut2 = sp.model('slab_nostokes_atmosout_cycle1_31_172143.nc')


fig, ax = plt.subplots(1,3)
ax0 = ax[0]
ax1 = ax[1]
ax2 = ax[2]

# axis info
rr = 40/1.5
ytick_pos = np.arange(0,(np.round((ymax-ymin)/rr)+1)*rr,rr)
ytick_lab = np.intc(np.round(ytick_pos*res, decimals = 0))
xtick_pos = np.arange(0,(np.round((xmax-xmin)/rr)+1)*rr,rr)
xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))
xtick_lab1 = np.intc(np.round(xtick_pos*res*(np.sqrt((cut1.x_coords[0]-cut1.x_coords[1])**2 +
                                                     (cut1.y_coords[0]-cut1.y_coords[1])**2)/cut1.loop_size),
                              decimals = 0))
xtick_lab2 = np.intc(np.round(xtick_pos*res*(np.sqrt((cut2.x_coords[0]-cut2.x_coords[1])**2 +
                                                     (cut2.y_coords[0]-cut2.y_coords[1])**2)/cut2.loop_size),
                              decimals = 0))
xtick_lab4 = np.intc(np.round(xtick_pos*res*(np.sqrt((cut0.x_coords[0]-cut0.x_coords[1])**2 +
                                                     (cut0.y_coords[0]-cut0.y_coords[1])**2)/cut0.loop_size),
                              decimals = 0))

tt = 30
ttick_pos = np.arange(0,(np.round((213)/tt)+1)*tt,tt)
ttick_lab = ttick_pos*cad/60

temp_min = 4.2e3

ax0.imshow(cut0.temp[0,:,:,tau], cmap = 'inferno', origin = 'lower', vmin = temp_min)
ax1.imshow(cut1.temp[0,:,:,tau], cmap = 'inferno', origin = 'lower', vmin = temp_min)
ax2.imshow(cut2.temp[0,:,:,tau], cmap = 'inferno', origin = 'lower', vmin = temp_min)
