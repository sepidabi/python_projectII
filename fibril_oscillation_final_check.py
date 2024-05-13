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

for ii in range(len(file_search(fibdir,'crispex*3950*.csav'))):
    print(ii, file_search(fibdir,'crispex*3950*.csav')[ii])
cc = input("Enter the number of the cut: ")

cut_file = (file_search(fibdir,'crispex*3950*.csav'))[cc]
cut = restore(fibdir+cut_file)
cube_raw = np.mean(cut.loop_slab[w_pos-3:w_pos+3,:,:],axis = 0).squeeze()*1e9

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
i_range = np.abs(wstd*np.std(cube)) # intensity range
vmin, vmax = np.mean(cube)-i_range, np.mean(cube)+i_range#0.95, 0.99 # imshow intensity cropping range

xx = cut.loop_slab.shape[2]

# plotting settings
fig = plt.figure(figsize = (8,9))
gs = gridspec.GridSpec(2,1)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])

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

for oo in range(len(file_search(outdir, "cut"+cut_file[-11:-5]+"*txt"))):

    per_los = np.zeros(0)
    per_pos = np.zeros(0)
    osc_fname = outdir+file_search(outdir, "cut"+cut_file[-11:-5]+"*txt")[oo]
    coord = np.loadtxt(osc_fname)
    amin, amax = int(np.min(coord[:,0])), int(np.max(coord[:,0]))+1
    tmin, tmax = int(np.min(coord[:,1])), int(np.max(coord[:,1]))+1
    trange = np.linspace(tmin, tmax, tmax-tmin+1)
    osc_fname_smth = outdir+"smth"+osc_fname[31:]
    smth_y = np.loadtxt(osc_fname_smth)
    
    ax1.plot(coord[:,1], coord[:,0],color = 'wheat', alpha = 0.5, label = 'Manual fit', linestyle = '--', linewidth = 1)
    ax1.plot(trange, smth_y, color = 'orange', alpha = 0.5, label = r'smthd I$_{\rm{max}}$', linewidth = 2)
    [y,x] = coord[0,:]
    ax1.text(x,y,osc_fname[-10:-4], color = 'orange', fontdict = font)
    plt.show()
    ucheck = input("Accept oscillation path? ([y] / n)", default="y")
    if ucheck =='y':
        # cut geometry
        x1 = cut.x_coords[0]
        x2 = cut.x_coords[1]
        if (cut_indx=='120700'):
            y2 = cut.y_coords[0]
            y1 = cut.y_coords[1]
        else:
            y1 = cut.y_coords[0]
            y2 = cut.y_coords[1]
        r = np.sqrt((x1-x2)**2+(y1-y2)**2)
        sinth = np.abs(y2-y1)/r
        costh = np.abs(x2-x1)/r
        theta = np.rad2deg(np.arcsin(sinth))
        
        # on the FOV map
        d = np.mean(coord[:,0])
        yi = d*sinth + y1
        xi = d*costh + x1
        mid_coord = np.array([xi,yi])
        osc_indx = str(jjj)
        osc_obj = (file_search(objdir, "cut*"+cut_indx+"*-no"+osc_indx+".obj"))[0]
        osc = load_obj(objdir+osc_obj)
        # categorization
        vpos_f = osc.vpos_f
        vlos_f = osc.vlos_f
        vpos_f_norm = (vpos_f-np.mean(vpos_f)) / (np.std(vpos_f)*len(vpos_f))
        vpos_ff_norm = (vpos_f-np.mean(vpos_f)) / (np.std(vpos_f))
        vlos_f_norm = (vlos_f-np.mean(vlos_f)) / (np.std(vlos_f))
        vlos_ff_norm = (vlos_f-np.mean(vlos_f)) / (np.std(vlos_f)*len(vlos_f))
        ncc = sgnl.correlate(vpos_f_norm, vlos_f_norm, mode = 'full')
        ncc_min = np.min(ncc)
        ncc_max = np.max(ncc)

        # OSCILLATION VALIDITY CHECK
        if ((len(osc.los_ymin)!=0 and len(osc.los_ymax)!=0)
            and (len(osc.pos_ymin)!=0 and len(osc.pos_ymax)!=0)):

            # average calculations
            #---------------------------

            # LOS Amplitude
            A_los_m = np.mean([np.mean(np.abs(osc.los_ymax)), np.mean(np.abs(osc.los_ymin))])
            # POS Amplitude
            A_pos_m = np.mean([np.mean(np.abs(osc.pos_ymax)), np.mean(np.abs(osc.pos_ymin))])
                    
            # LOS Period
            pe_l = np.sort(np.append(osc.los_xmax, osc.los_xmin))
            
            if (len(pe_l)==1):
                per_los_one = pe_l*2
            if (len(pe_l)==2):
                per_los_one = 2*np.abs(osc.los_xmax-osc.los_xmin)
            if (len(pe_l)==3):
                per_los_one = (pe_l[-1] - pe_l[0])
            if (len(pe_l)>3):
                for l in range (len(pe_l)-3):
                    per_los_one = (pe_l[l+3]-pe_l[l+1]+pe_l[l+2]-pe_l[l])/2.
            per_los = np.append(per_los, per_los_one)
            per_los_m =  np.mean(per_los)
            per_los_dom = np.max(per_los)    
            #per_los_dom = np.append(per_los_dom, per_los_one[np.argmin(np.abs(per_los_one-osc.los_ncc_xmax[0]))])
            #A_los_m = np.append(A_los_m, np.mean([np.mean(np.abs(osc.los_ymax)), np.mean(np.abs(osc.los_ymin))]))
            
            
            # POS period
            pe_p = np.sort(np.append(osc.pos_xmax, osc.pos_xmin))
            
            if (len(pe_p)==1):
                per_pos_one = pe_p*2
            if (len(pe_p)==2):
                per_pos_one = 2*np.abs(osc.pos_xmax-osc.pos_xmin)
            if (len(pe_p)==3):
                per_pos_one = (pe_p[-1] - pe_p[0])
            if (len(pe_p)>3):
                for m in range (len(pe_p)-3):
                    per_pos_one = (pe_p[m+3]-pe_p[m+1]+pe_p[m+2]-pe_p[m])/2.
                    
            per_pos = np.append(per_pos, per_pos_one) 
            per_pos_dom = np.max(per_pos)    
            per_pos_m = np.mean(per_pos)
            #per_pos_dom = np.append(per_pos_dom, per_pos_one[np.argmin(np.abs(per_pos_one-osc.pos_ncc_xmax[0]))])
            #A_pos_m = np.append(A_pos_m, np.mean([np.mean(np.abs(osc.pos_ymax)), np.mean(np.abs(osc.pos_ymin))]))
            
            
            ssize = 30 # plotting marker size
            
            per_dom = np.max([per_pos_dom, per_los_dom])
            per_m = np.mean([per_pos_m, per_los_m])
            lag = np.arange(-len(vlos_f_norm)+1, len(vpos_f_norm))*8.

            Ppos_m = per_pos_m
            Ppos_dom = per_pos_dom
            Plos_m = per_los_m
            Plos_dom = per_los_dom
            dp = np.round(np.abs(per_pos_m - per_los_m), decimals = 1)
            #dp_dom = np.round(np.abs(per_pos_dom - per_los_dom), decimals = 1)

            if (np.max(ncc)>corr_thre
                and np.abs(lag[np.argmax(ncc)])<=per_dom
                and np.max(ncc)>np.abs(np.min(ncc))):
                
                P_cor = np.abs(lag[np.argmax(ncc)])
                ind = ind+1
                '''print(
                    str(ind)+')'+
                    osc_fname[31:-4]+ '          ' +
                    str(dp) + '          ' +
                    #str(dp_dom) + '          ' +
                    str(np.round(Ppos_m, decimals = 1)) + '          ' +
                    str(np.round(Ppos_dom, decimals = 1)) + '          ' +
                    str(np.round(Plos_m, decimals = 1))+ '          ' +
                    str(np.round(Plos_dom, decimals = 1)) + '          ' +
                    str(np.round(P_cor
                                 , decimals = 1))
                )'''

                ucheck2 = input('Still yes? ([1] / 0)', default = 1)
                if ucheck2==1:
                    dphi = lag[np.argmax(ncc)]/per_dom
                if ucheck2==0:
                    dphi = -5
                    
            elif (np.min(ncc)<-corr_thre
                  and np.abs(lag[np.argmin(ncc)])<=per_dom
                  and np.max(ncc)<np.abs(np.min(ncc))):

                ind = ind+1

                P_cor = np.abs(lag[np.argmin(ncc)])
                
                ucheck2 = input('Still yes? ([1] / 0)', default = 1)
                if ucheck2==1:
                    dphi = lag[np.argmin(ncc)]/per_dom
                    #if (dphi_one<0):
                    #   dphi_one = 2+dphi_one
                else:
                    dphi = -5

            else:
                dphi = -5
        else:
            dphi = -10
            per_m = 0

        ################
        # filling in the result object
        ################
        ratio = 0
        dist = 0
        cor_rate = 0
        
        result = oscillation_result(per_los,
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
                                    )
        objresname = objresdir + 'result' + osc_fname[-30:-4] + "-no"+ str(jjj)+ ".obj"
        save_obj(objresname, result)
        
    else:
        ucheck3 = input('delete the oscillation? ([NO!] / yes)', default = 'NO!')
        if ucheck3=='NO!':
            print('ok ;)')
        if ucheck=='yes':
            
    stop()

