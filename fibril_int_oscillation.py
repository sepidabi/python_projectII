# fibril_hor_oscillation.py

'''Meant for plotting a sub field of view
with a fibril and its horizontal oscillations
(the cuts along the fibrils are already defined via crispex)
'''

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

# Reading the observed data
pref = ['6173','8542','3950']
file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
file_ca8 =file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak =file_search(datadir,'crispex*'+pref[2]+'*.fcube')

cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])
cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])

nx = cube_cak.shape[4]
ny = cube_cak.shape[3]
nt = cube_cak.shape[0]
nw = cube_cak.shape[2]

# getting the fibril and cuts information
ind = 4
fibdir = datadir+'fr29/'
files_n = len(file_search(fibdir,'crispex*3950*.csav'))/4

for ff in range(files_n):

    fib_file = (file_search(fibdir,'crispex*3950*.csav'))[0+4*ff]
    fib = restore(fibdir+fib_file)
    cut1_file = (file_search(fibdir,'crispex*3950*.csav'))[1+4*ff]
    cut2_file = (file_search(fibdir,'crispex*3950*.csav'))[2+4*ff]
    cut3_file = (file_search(fibdir,'crispex*3950*.csav'))[3+4*ff]
    #cut4_file = (file_search(fibdir,'crispex*3950*.csav'))[4]
    
    cut1 = restore(fibdir+cut1_file)
    cut2 = restore(fibdir+cut2_file)
    cut3 = restore(fibdir+cut3_file)
    
    xmin = int(np.min([np.min([cut1.x_coords, cut2.x_coords,
                               cut3.x_coords]),np.min(fib.x_coords)]) - 5)
    ymin = int(np.min([np.min([cut1.y_coords, cut2.y_coords,
                               cut3.y_coords]),np.min(fib.y_coords)]) - 5)
    xmax = int(np.max([np.max([cut1.x_coords, cut2.x_coords,
                               cut3.x_coords]),np.max(fib.x_coords)]) + 5)
    ymax = int(np.max([np.max([cut1.y_coords, cut2.y_coords,
                               cut3.y_coords]),np.max(fib.y_coords)]) + 5)
    
    #######################
    ## plot oscilations
    #######################
    
    xfig = 8
    yfig = 4.97
    plt.close('all')
    f = plt.figure(figsize = (xfig, yfig))
    
# plot settings    
    gs1 = gridspec.GridSpec(2, 2) # grid scale of the left col.
    gs1.update(left=0.06,
               right=0.46,
               wspace=0.,
               bottom=0.1,
               top=0.98,
               hspace = 0.
    )
    
    
    gs = gridspec.GridSpec(1, 3) # grid scale of the right col.
    gs.update(left=0.53,
              right=0.99,
              wspace=0.02,
              bottom=0.1,
              top=0.98,
              hspace = 0.
    )
    
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
    xtick_lab4 = np.intc(np.round(xtick_pos*res*(np.sqrt((cut3.x_coords[0]-cut3.x_coords[1])**2 +
                                                         (cut3.y_coords[0]-cut3.y_coords[1])**2)/cut3.loop_size),
                                  decimals = 0))
    
    tt = 30
    ttick_pos = np.arange(0,(np.round((213)/tt)+1)*tt,tt)
    ttick_lab = ttick_pos*cad/60
    
    
    ax13 = plt.subplot(gs1[1,0],adjustable = 'box')#,aspect = 'equal')
    ax13.set_xlabel(r'x [arcsec]', fontdict = font)
    ax13.set_ylabel(r'y [arcsec]', fontdict = font)
    ax13.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax13.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax13.tick_params(which='minor', length=2)
    ax13.set_xticks(xtick_pos)
    ax13.set_xticklabels(xtick_lab, fontdict = font)
    ax13.set_yticks(ytick_pos)
    ax13.set_yticklabels(ytick_lab, fontdict = font)
    
    ax14 = plt.subplot(gs1[1,1],adjustable = 'box')#,aspect = 'equal')
    ax14.set_xlabel(r'x [arcsec]', fontdict = font)
    #ax14.set_ylabel(r'y [arcsec]', fontdict = font)
    ax14.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax14.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax14.tick_params(which='minor', length=2)
    ax14.set_xticks(xtick_pos)
    ax14.set_xticklabels(xtick_lab, fontdict = font)
    ax14.set_yticks(ytick_pos)
    ax14.set_yticklabels([])
    
    ax11 = plt.subplot(gs1[0,0],adjustable = 'box')#,aspect = 'equal')
    #ax11.set_xlabel(r'x [arcsec]', fontdict = font)
    ax11.set_ylabel(r'y [arcsec]', fontdict = font)
    ax11.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax11.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax11.tick_params(which='minor', length=2)
    ax11.set_xticks(xtick_pos)
    ax11.set_xticklabels([])
    ax11.set_yticks(ytick_pos)
    ax11.set_yticklabels(ytick_lab, fontdict = font)
    
    ax12 = plt.subplot(gs1[0,1],adjustable = 'box')#,aspect = 'equal')
    #ax12.set_xlabel(r'x [arcsec]', fontdict = font)
    #ax12.set_ylabel(r'y [arcsec]', fontdict = font)
    ax12.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax12.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax12.tick_params(which='minor', length=2)
    ax12.set_xticks(xtick_pos)
    ax12.set_xticklabels([])
    ax12.set_yticks(ytick_pos)
    ax12.set_yticklabels([])
    
    
    ax2 = plt.subplot(gs[0,0],adjustable = 'box')#,aspect = 'equal')
    ax2.set_xlabel(r'length [arcsec]', fontdict = font)
    ax2.set_yticks(ttick_pos)
    ax2.set_yticklabels(ttick_lab, fontdict = font)
    ax2.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.set_ylabel(r'$t$ [m]', fontdict = font)
    ax2.set_xticks(xtick_pos)
    ax2.set_xticklabels(xtick_lab4, fontdict = font)
    ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.set_ylim(0,nt-1)
    
    ax3 = plt.subplot(gs[0,1],adjustable = 'box')#,aspect = 'equal')
    ax3.set_yticks([])
    ax3.set_xlabel(r'length [arcsec]', fontdict = font)
    ax3.set_yticks(ttick_pos)
    ax3.set_yticklabels([])
    ax3.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax3.set_xticklabels(xtick_lab2, fontdict = font)
    ax3.set_xticks(xtick_pos)
    ax3.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax3.set_ylim(0,nt-1)
    
    ax4 = plt.subplot(gs[0,2],adjustable = 'box')#,aspect = 'equal')
    ax4.set_yticks([])
    ax4.set_xlabel(r'length [arcsec]', fontdict = font)
    ax4.set_yticks(ttick_pos)
    ax4.set_yticklabels([])
    ax4.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax4.set_xticklabels(xtick_lab1, fontdict = font)
    ax4.set_xticks(xtick_pos)
    ax4.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax4.set_ylim(0,nt-1)
    
    #######################
    ## PLOTS roi
    #######################
    vmin = None#2.e-10 #1.9180524e-10
    vmax = None#1.e-8 #1.09262785e-08
    
    ax11.imshow(cube_cak[fr_select[0],0,w_pos,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower',
                vmin = vmin, vmax=vmax)
    ax11.text(1, 1, str(fr_select[0]*8/60)+' m', color = 'red')
    
    ax12.imshow(cube_cak[fr_select[1],0,w_pos,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower',
                vmin = vmin, vmax=vmax)
    ax12.text(1, 1, str(fr_select[1]*8/60)+' m', color = 'red')
    
    ax13.imshow(cube_cak[fr_select[2],0,w_pos,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower',
                vmin = vmin, vmax=vmax)
    ax13.text(1, 1, str(fr_select[2]*8/60)+' m', color = 'red')
     
    ax14.imshow(cube_cak[fr_select[3],0,w_pos,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower',
                vmin = vmin, vmax=vmax)
    ax14.text(1, 1, str(fr_select[3]*8/60)+' m', color = 'red')
    #ax14.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'red')
    ax14.plot(cut1.x_coords-xmin, cut1.y_coords-ymin, color = 'red')
    ax14.plot(cut2.x_coords-xmin, cut2.y_coords-ymin, color = 'red')
    ax14.plot(cut3.x_coords-xmin, cut3.y_coords-ymin, color = 'red')
    ax14.text(cut1.x_coords[0]-xmin-4, cut1.y_coords[0]-ymin-4, '(3)', fontdict=font, color = 'red')
    ax14.text(cut2.x_coords[0]-xmin-4, cut2.y_coords[0]-ymin-4, '(2)', fontdict=font, color = 'red')
    ax14.text(cut3.x_coords[0]-xmin-4, cut3.y_coords[0]-ymin-4, '(1)', fontdict=font, color = 'red')
    
    
    ratio =0.8 # DO NOT CHANGE THIS RATIO otherwise the output image goes bananas!
    ax2.imshow(cut3.loop_slab[w_pos,:,:], cmap = 'gray', origin = 'lower', aspect=ratio,
               vmin = vmin, vmax=vmax)
    ax2.plot(fr_select[0], alpha = 0.5, linestyle = '--')
    ax2.text(2,205,'(1)', color = 'red')
    
    ax3.imshow(cut2.loop_slab[w_pos,:,:], cmap = 'gray', origin = 'lower', aspect=ratio,
               vmin = vmin, vmax=vmax)
    ax3.text(2,205,'(2)', color = 'red')
    
    ax4.imshow(cut1.loop_slab[w_pos,:,:], cmap = 'gray', origin = 'lower', aspect=ratio,
               vmin = vmin, vmax=vmax)
    ax4.text(2,205,'(3)', color = 'red')
    
    plt.show()
    
    # wpos
    if (w_pos == 13): wpos='K3'
    elif (w_pos==10): wpos = 'K2V'
    else: wpos = 'K2R'
    
    filename = outdir+'oscilation'+wpos+'_'+str(ff)+'.pdf'
    f.savefig(filename, quality = 100)
    f.tight_layout()
    print 'file saved to: ' + filename


    
