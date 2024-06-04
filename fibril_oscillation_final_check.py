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
from scipy.optimize import curve_fit
from datetime import date
from datetime import datetime
import scipy.ndimage as spnd
import scipy.signal as sgnl
from mpl_toolkits import mplot3d
from scipy.signal import argrelextrema as extrem
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import os

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
objdir = '/home/seki2695/OUTPUT/project2/class/'
objresdir = '/home/seki2695/OUTPUT/project2/result_object/'
removedir = outdir+'deleted_oscillations'

# Difinitions
res = 0.0375 # CHROMIS pixel size in arcsec
cad = 8
#fr_select = [29,75,135,167] # frames to display
w_pos = 13 # wavelength position to display (K2V, K3, K2R) = (10,13,16)

plt.close("all")

class oscillation:
    def __init__(self, time, vpos_f, vlos_f, pos_xmax, pos_ymax, pos_xmin, pos_ymin, pos_ncc_xmax, pos_ncc_ymax, los_xmax, los_ymax, los_xmin, los_ymin, los_ncc_xmax, los_ncc_ymax):
        self.time = time
        self.vpos_f = vpos_f
        self.vlos_f = vlos_f
        self.pos_xmax = pos_xmax
        self.pos_ymax = pos_ymax
        self.pos_xmin = pos_xmin
        self.pos_ymin = pos_ymin
        self.pos_ncc_xmax = pos_ncc_xmax
        self.pos_ncc_ymax = pos_ncc_ymax
        self.los_xmax = los_xmax
        self.los_ymax = los_ymax
        self.los_xmin = los_xmin
        self.los_ymin = los_ymin
        self.los_ncc_xmax = los_ncc_xmax
        self.los_ncc_ymax = los_ncc_ymax


class oscillation_result:
    def __init__(self,
                 per_los,
                 per_los_m,
                 per_los_dom,
                 per_pos,
                 per_pos_m,
                 per_pos_dom,
                 A_los_m,
                 A_pos_m,
                 dphi,
                 ncc_min,
                 ncc_max,
                 mid_coord,
                 cor_rate,
                 theta,
                 per_m,
                 per_dom,
                 ratio,
                 dist,
                 ):
        
        self.per_los = per_los
        self.per_los_m = per_los_m
        self.per_los_dom = per_los_dom
        self.per_pos = per_pos
        self.per_pos_m = per_pos_m
        self.per_pos_dom = per_pos_dom
        self.A_los_m = A_los_m
        self.A_pos_m = A_pos_m
        self.dphi = dphi
        self.ncc_min = ncc_min
        self.ncc_max =  ncc_max
        self.mid_coord = mid_coord
        self.cor_rate = cor_rate
        self.theta = theta
        self.per_m = per_m
        self.per_dom = per_dom
        self.ratio = ratio
        self.dist = dist


# Reading the observed data
pref = ['6173','8542','3950']
file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
file_ca8 =file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak =file_search(datadir,'crispex*'+pref[2]+'*.fcube')

cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])
cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])

fibdir = datadir+'big_cut/'

# plotting elements
tt = 30
rr = 40/1.5

# the cross-correlation threshold 
corr_thre = 0.5 # above which the oscillations are assumed to be (anti)correlated

# plotting settings
fig = plt.figure(figsize = (16,9))
gs = gridspec.GridSpec(3,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])
ax3 = plt.subplot(gs[2,0])
ax4 = plt.subplot(gs[0,1])
ax5 = plt.subplot(gs[1,1])
ax6 = plt.subplot(gs[2,1])

files = file_search(fibdir,'crispex*3950*.csav')
for ii in range(len(files)):
    print(ii, files[ii])
i0 = input("Enter the number of the cut: ")
for i in range(i0,len(files)):

    cut_file = files[i]
    cut = restore(fibdir+cut_file)
    cube_raw = np.mean(cut.loop_slab[w_pos-3:w_pos+3,:,:],axis = 0).squeeze()*1e9
    oscils = file_search(outdir, "cut"+cut_file[-11:-5]+"*txt")

    # cleaning and sharpenning the intensity map
    wg = 0.15 # weight of the background intensity from the wings
    cube_bg = wg*np.mean(cut.loop_slab[:w_pos-6,:,:]+cut.loop_slab[w_pos+8:,:,:],axis=0).squeeze()*1e9
    cube_trunc = cube_raw - cube_bg # removing the background intensity
    # contrast correction
    cube_med = exposure.rescale_intensity(sgnl.medfilt2d(cube_trunc,kernel_size = [3,1]), out_range=(0, 1.))
    # sharpenning the maps alond the time axis to eliminate the jitters within frames
    cube_sharp = exposure.rescale_intensity(sharpen(cube_trunc, sigma =[3,1], unsigma = [1,3]), out_range=(0, 1.))
    # gamma correcting to enhance the intensity
    cube = exposure.rescale_intensity(exposure.adjust_gamma(cube_sharp,0.1), out_range=(0, 1.))
    wstd = 3 # weight of the standard deviation
    i_range = wstd*np.abs(np.std(cube)) # intensity range
    i_range2 = wstd*np.abs(np.std(cube))/2
    vmin, vmax = np.mean(cube)-i_range, np.mean(cube)+i_range#0.95, 0.99 # imshow intensity cropping range
    vmin2,vmax2 = np.mean(cube)-i_range2, np.mean(cube)+i_range2
    xx = cut.loop_slab.shape[2]

    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    ax5.clear()
    ax6.clear()

    # Intensity
    #------------
    ax1.imshow(np.transpose(cube),
               cmap = 'gray', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               vmin=vmin,
               vmax = vmax,
               aspect = 0.6,
               )
    ax2.imshow(np.transpose(cube),
               cmap = 'gray', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               vmin=vmin,
               vmax = vmax,
               aspect = 0.6,
               )
    ax3.imshow(np.transpose(cube),
               cmap = 'gray', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               vmin=vmin2,
               vmax = vmax2,
               aspect = 0.6,
               )

    oo = 0
    for oo in range(0,len(oscils)):
        print ([oo,(oscils[oo])[-10:-4]])
    o0 = raw_input(str(i)+"-cut"+cut_file[-11:-5]+": start from oscillation no: (default=0) ").strip() or '0'
    o0 = int(o0)
    for o in range(int(o0),len(oscils)):
        
        per_los = np.zeros(0)
        per_pos = np.zeros(0)
        osc_fname = outdir+oscils[o]
        print("==========================================")
        print(" OSCILLATION: " + osc_fname[31:-4] + " ")
        print("==========================================")
        coord = np.loadtxt(osc_fname)
        
        # Extract the oscillation objects
        obj_name = file_search(objdir, osc_fname[-30:-4]+"*.obj")[0]
        obj = load_obj(objdir + obj_name)
        objres_name = file_search(objresdir, 'result' + osc_fname[-30:-4]+"*.obj")[0]
        objres = load_obj(objresdir + objres_name)
        
        #amin, amax = int(np.min(coord[:,0])), int(np.max(coord[:,0]))+1
        tmin, tmax = int(np.min(coord[:,1])), int(np.max(coord[:,1]))+1
        trange = np.linspace(tmin, tmax, tmax-tmin+1)
        osc_fname_smth = outdir+"smth"+osc_fname[31:]
        smth_y = np.loadtxt(osc_fname_smth)
        
        # plots for sanity check
        #------------------------------
        plt.title(osc_fname[31:-4])
        # clearing the subplots
        ax4.clear()
        ax5.clear()
        ax6.clear()
        if o>o0:
            (over3.pop(0)).remove()
            (over2.pop(0)).remove()
        
        # intensity plots
        ax1.set_title(osc_fname[-10:-4])
        ax1.plot(coord[:,1], coord[:,0],color = 'wheat', alpha = 0.5, label = 'Manual fit', linestyle = '--', linewidth = 1)
        ax1.plot(trange, smth_y, color = 'orangered', alpha = 0.5, label = r'smthd I$_{\rm{max}}$', linewidth = 2)
        over2 = ax2.plot([trange[0], trange[-1]], [smth_y[0], smth_y[-1]], 'r.', alpha = 0.5)
        over3 = ax3.plot([trange[0], trange[-1]], [smth_y[0], smth_y[-1]], 'r.', alpha = 0.5)
        [y,x] = coord[0,:]
        ax1.text(x,y,osc_fname[-10:-4], color = 'orangered')
        
        
        # plot LOS and POS
        ax4.axhline(0, linestyle = ":", color = "gray", linewidth = 1, alpha = 0.5)
        ax4.plot(obj.vlos_f, label = r"LOS", color = 'orangered') # smoothed vLOS
        ax4.plot(obj.vpos_f, label = r"POS", color = 'darkgray') # smoothed vLOS
        ax4.legend()
        
        # Auto-correlation plots
        ax5.axhline(0, linestyle = ":", color = "gray", linewidth = 1, alpha = 0.5)
        ax5.plot(obj.time, obj.ncc_pos, color = "gray", alpha = 0.75)
        ax5.plot(obj.time, obj.ncc_los, color = "orangered", alpha = 0.75)    

        # Cross-correlation plots
        if objres.dphi>=-5:
            ax6.axhline(0, linestyle = ":", color = "gray", linewidth = 1, alpha = 0.5)
            ax6.axvline(0, linestyle = ":", color = "gray", linewidth = 1, alpha = 0.5)
            ax6.plot(objres.ncc, label = r"NCC", color = "black", alpha = 0.75)
            ax6.legend()
            
        plt.tight_layout()
        plt.show()
        
        ucheck1 = raw_input("Accept oscillation path ([y] / n)? ").strip() or "y"
        if ucheck1 =='y':
            print('The oscillation properties are accepted as it was.')
            #stop()
        if ucheck1=='n':
            ucheck2 = raw_input('Delete the oscillation ([NO!] / yes)? ').strip() or "NO!"
            if ucheck2=='NO!':
                ucheck3 = raw_input('ok then ;) Shall we redefine the dphi? (y / [n])').strip() or "n"
                if ucheck3=='y':
                    print('Changing the dphi value from '+str(objres.dphi) + ' to:')
                    objres.dphi = input('Enter the new value: ')
                    save_obj(objresdir + objres_name, objres)
                if ucheck3=='n':
                    print('So nothing was changed...')
                    stop()
                    
            if ucheck2=='yes':
                os.system('mv -v '+objresdir+'*'+osc_fname[-30:-4]+'* '+removedir)
                os.system('mv -v '+outdir+'*'+osc_fname[-30:-4]+'* '+removedir)
                os.system('mv -v '+objdir+'*'+osc_fname[-30:-4]+'* '+removedir)
                #stop()

        # intensity plots
        ax1.plot(trange, smth_y, color = 'white', alpha = 0.5, label = r'smthd I$_{\rm{max}}$', linewidth = 2)
        ax1.text(x,y,osc_fname[-10:-4], color = 'white')
        
