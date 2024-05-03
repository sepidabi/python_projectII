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

# Difinitions
res = 0.0375 # CHROMIS pixel size in arcsec
cad = 8
#fr_select = [29,75,135,167] # frames to display
w_pos = 13 # wavelength position to display (K2V, K3, K2R) = (10,13,16)

plt.close("all")

#functions
def test_func(x, a, b, c, d):
    return a * np.sin(b * x - c)+d

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


# boxcar smoothing

def smooth(x,window_len=11,window='flat'):

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2):-(window_len/2)]


# Reading the observed data
pref = ['6173','8542','3950']
file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
file_ca8 =file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak =file_search(datadir,'crispex*'+pref[2]+'*.fcube')

cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])
cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])

# CaK Core
cak_wc = 13
fr = 59
cak_core = unsharp(cube_cak[fr,0, cak_wc, :, :], alpha = 0.5, sigma = 5.)
# photosphere CPtot
fe_cont = cube_fe[fr,0,-1,:,:]
CPtot = unsharp(0.5*(np.mean(cube_fe[fr,3,0:7,:,:], axis=0) - np.mean(cube_fe[fr,3,7:-1,:,:], axis=0))/fe_cont, alpha = 0.5, sigma = 2.)

# plotting settings
fig = plt.figure(figsize = (8,3))
gs = gridspec.GridSpec(ncols = 1,nrows = 1, figure = fig)
ax0 = plt.subplot(gs[0,0])

nx = cube_cak.shape[4]
ny = cube_cak.shape[3]
nt = cube_cak.shape[0]
nw = cube_cak.shape[2]

fibdir = datadir+'big_cut/'

# to test the curve fitting
ti = np.arange(0, nt, 1.)

# plotting elements
tt = 30
rr = 40/1.5

tau_i = 23

per_los = np.zeros(0)
per_los_m = np.zeros(0)
per_los_dom = np.zeros(0)
per_pos = np.zeros(0)
per_pos_m = np.zeros(0)
per_pos_dom = np.zeros(0)
amp_los = np.zeros(0)
amp_los_m = np.zeros(0)
amp_pos = np.zeros(0)
amp_pos_m = np.zeros(0)
A_los_one = np.zeros(0)
A_los = np.zeros(0)
A_los_m = np.zeros(0)
A_pos_one = np.zeros(0)
A_pos = np.zeros(0)
A_pos_m = np.zeros(0)
dphi = np.zeros(0)

ncc_min = np.zeros(0)
ncc_max =  np.zeros(0)

dphi = np.zeros(0)
mid_coord_x = np.zeros(0)
mid_coord_y = np.zeros(0)
cor_rate = np.zeros(0)
theta = np.zeros(0)
per_dom = np.zeros(0)

# Extract the oscillation objects

# User contribution to the specify the oscillation
for iii in range(len(file_search(objdir, "cut*osc0.obj"))):
    cut_indx = file_search(objdir, "cut*osc0.obj")[iii][3:12]
    
    # coords of the centere of oscillation
    cut_file = (file_search(fibdir,'crispex*3950*'+cut_indx+'*.csav'))[0]
    cut = restore(fibdir+cut_file)
    cube = cut.loop_slab[w_pos,:,:]

    per_los_one = np.zeros(0)
    per_pos_one = np.zeros(0)
    
    for jjj in range(len(file_search(objdir, 'cut'+cut_indx+'*osc*.obj'))):
        
        osc_fname = outdir+file_search(outdir, "cut"+cut_file[-11:-5]+"*txt")[jjj]
        coord = np.loadtxt(osc_fname)

        # cut geometry
        x1 = cut.x_coords[0]
        x2 = cut.x_coords[1]
        if (cut_indx=='08_120700'):
            y2 = cut.y_coords[0]
            y1 = cut.y_coords[1]
        else:
            y1 = cut.y_coords[0]
            y2 = cut.y_coords[1]
        r = np.sqrt((x1-x2)**2+(y1-y2)**2)
        sinth = np.abs(y2-y1)/r
        costh = np.abs(x2-x1)/r
        theta = np.append(theta, np.rad2deg(np.arcsin(sinth)))
        
        # on the FOV map
        # sin = (y2-y1)/r = (yi-y1)/d
        d = np.mean(coord[:,0])
        yi = d*sinth + y1
        xi = d*costh + x1
        osc_mid_coord = np.array([xi,yi])
        osc_indx = str(jjj)
        osc_obj = (file_search(objdir, "cut*"+cut_indx+"_osc"+osc_indx+".obj"))[0]
        osc = load_obj(objdir+osc_obj)
        # categorization
        vpos_f = osc.vpos_f
        vlos_f = osc.vlos_f
        vpos_f_norm = (vpos_f-np.mean(vpos_f)) / (np.std(vpos_f)*len(vpos_f))
        vpos_ff_norm = (vpos_f-np.mean(vpos_f)) / (np.std(vpos_f))
        vlos_f_norm = (vlos_f-np.mean(vlos_f)) / (np.std(vlos_f))
        vlos_ff_norm = (vlos_f-np.mean(vlos_f)) / (np.std(vlos_f)*len(vlos_f))
        ncc = sgnl.correlate(vpos_f_norm, vlos_f_norm, mode = 'full')
        ncc_min = np.append(ncc_min, np.min(ncc))
        ncc_max = np.append(ncc_max, np.max(ncc))
        
        # average calculations
        if ((len(osc.los_ymin)!=0 and len(osc.los_ymax)!=0)
            and (len(osc.pos_ymin)!=0 and len(osc.pos_ymax)!=0)):
            # LOS Amplitude
            #A_los_m = np.append(A_los_m, np.mean([np.mean(np.abs(osc.los_ymax)), np.mean(np.abs(osc.los_ymin))]))
            # POS Amplitude
            #A_pos_m = np.append(A_pos_m, np.mean([np.mean(np.abs(osc.pos_ymax)), np.mean(np.abs(osc.pos_ymin))]))
            
        # OSCILLATION VALIDITY CHECK
        
        #if (len(osc.los_ncc_xmax)!=0 and len(osc.pos_ncc_xmax)!=0):
            #print(cut_indx, osc_fname)
            # LOS Period
            pe_l = np.sort(np.append(osc.los_xmax, osc.los_xmin))
            
            if (len(pe_l)==2):
                per_los_one = np.append(per_los_one, 2*2*np.abs(osc.los_xmax-osc.los_xmin))
            if (len(pe_l)==3):
                per_los_one = np.append(per_los_one, 2*(pe_l[-1] - pe_l[0]))
            if (len(pe_l)>3):
                for l in range (len(pe_l)-3):
                    per_los_one = np.append(per_los_one, (pe_l[l+3]-pe_l[l+1]+pe_l[l+2]-pe_l[l])/2.)
            #if (len(pe_l)==1):
             #   per_los_one = pe_l*2
                
            per_los = np.append(per_los, per_los_one)
            #per_los_dom = np.append(per_los_dom, per_los_one[np.argmin(np.abs(per_los_one-osc.los_ncc_xmax[0]))])
            per_los_m = np.append(per_los_m, np.mean(per_los_one))
            A_los_m = np.append(A_los_m, np.mean([np.mean(np.abs(osc.los_ymax)), np.mean(np.abs(osc.los_ymin))]))
            
            
            # POS period
            pe_p = np.sort(np.append(osc.pos_xmax, osc.pos_xmin))
            
            if (len(pe_p)==2):
                per_pos_one = np.append(per_pos_one, 2*2*np.abs(osc.pos_xmax-osc.pos_xmin))
            if (len(pe_p)==3):
                per_pos_one = np.append(per_pos_one, 2*(pe_p[-1] - pe_p[0]))
            if (len(pe_p)>3):
                for m in range (len(pe_p)-3):
                    per_pos_one = np.append(per_pos_one, (pe_p[m+3]-pe_p[m+1]+pe_p[m+2]-pe_p[m])/2.)
                    
            per_pos = np.append(per_pos, per_pos_one) 
            #per_pos_dom = np.append(per_pos_dom, per_pos_one[np.argmin(np.abs(per_pos_one-osc.pos_ncc_xmax[0]))])
            per_pos_m = np.append(per_pos_m, np.mean(per_pos_one))
            A_pos_m = np.append(A_pos_m, np.mean([np.mean(np.abs(osc.pos_ymax)), np.mean(np.abs(osc.pos_ymin))]))
            
            
            ssize = 30 # plotting marker size
            
            per_one_dom = np.max([np.max(per_pos_one), np.max(per_los_one)])
            per_dom = np.append(per_dom, np.mean([np.mean(per_pos_one), np.mean(per_los_one)]))
            lag = np.arange(-len(vlos_f_norm)+1, len(vpos_f_norm))*8.
            
            if (np.max(ncc)>0.5 and np.abs(2*lag[np.argmax(ncc)])<per_one_dom):
                dphi_one = 2*lag[np.argmax(ncc)]/per_one_dom
                #if (dphi_one<0):
                 #   dphi_one = 2+dphi_one
                dphi = np.append(dphi, dphi_one)
                #print('max', str(iii)+'-'+str(jjj), 2*lag[np.argmax(ncc)], per_one_dom, dphi_one)
            elif (np.min(ncc)<-0.5 and np.abs(2*lag[np.argmin(ncc)])<per_one_dom):
                dphi_one = 2*lag[np.argmin(ncc)]/per_one_dom
                #if (dphi_one<0):
                 #   dphi_one = 2+dphi_one
                dphi = np.append(dphi, dphi_one)
                #print('min', str(iii)+'-'+str(jjj), 2*lag[np.argmin(ncc)], per_one_dom, dphi_one)
            else:
                dphi = np.append(dphi, -5)
        else:
            dphi = np.append(dphi,-10)
            per_dom = np.append(0,per_dom)

        mid_coord_x = np.append(mid_coord_x, osc_mid_coord[0])
        mid_coord_y = np.append(mid_coord_y, osc_mid_coord[1])

        # color coding the scatter plot
        #if (np.abs(np.min(ncc))>np.max(ncc)):
         #   cor_rate = np.append(cor_rate, np.min(ncc))
        
        #if (np.abs(np.min(ncc))<=np.max(ncc)):
         #   mid_coord_x = np.append(mid_coord_x, osc_mid_coord[0])
          #  mid_coord_y = np.append(mid_coord_y, osc_mid_coord[1])
           # cor_rate = np.append(cor_rate, np.max(ncc))
            
#stop()            
            
##########
# JOINT PLOTS
##########

# Period PLOT
# ========
fontsize = 8.
xy = np.linspace(-0,1000,100)

plt.close('all')
fig = plt.figure(figsize = (4,4))
bins = 50
gs = gridspec.GridSpec(ncols = 5,nrows = 5)

ax_joint = fig.add_subplot(gs[1:5,0:4])
ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_y = fig.add_subplot(gs[1:5,4])

xmin, ymin = 110, 110
xmax, ymax = 550, 550

# Period
x = per_pos_m
y = per_los_m


angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
ax_joint.text(xmax - 60,ymax - 55,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

ax_joint.scatter(x, y, alpha = 0.25, color = 'orangered')
ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = bins)
ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = bins)

# Set ax limits on marginals
ax_joint.set_xlim([xmin,xmax])
ax_joint.set_ylim([ymin,ymax])
ax_marg_y.set_ylim([ymin,ymax])
ax_marg_x.set_xlim([xmin,xmax])

pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
ax_joint.annotate('pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')

# Turn off ticks and labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
plt.setp(ax_marg_y.get_xticklabels(), visible=False)
plt.setp(ax_marg_x.get_yticklabels(), visible=False)
ax_marg_x.get_yaxis().set_visible(False)
ax_marg_y.get_xaxis().set_visible(False)

# horizontal hist
ax_marg_x.spines['top'].set_visible(False)
ax_marg_x.spines['right'].set_visible(False)
ax_marg_x.spines['left'].set_visible(False)

# vertical hist
ax_marg_y.spines['top'].set_visible(False)
ax_marg_y.spines['right'].set_visible(False)
ax_marg_y.spines['bottom'].set_visible(False)

# joint plot
ax_joint.spines['top'].set_visible(False)
ax_joint.spines['right'].set_visible(False)

# Set labels on joint
ax_joint.set_xlabel(r' $\overline{P}_{\rm{POS}}$ [s]', fontdict = font)
ax_joint.set_ylabel(r' $\overline{P}_{\rm{LOS}}$ [s]', fontdict = font)
ax_joint.tick_params(labelsize = fontsize)
ax_marg_x.tick_params(labelsize = fontsize)
ax_marg_y.tick_params(labelsize = fontsize)
ax_joint.tick_params(labelsize = fontsize)

plt.tight_layout()
plt.subplots_adjust(left = 0.14,
                    bottom = 0.11,
                    right = 1.,
                    top = 1.
)

filename = outdir + 'analys_pp.pdf'
fig.savefig(filename, dpi = 1000)

# Amplitude PLOT
# ==========
fontsize = 8.
xy = np.linspace(-0,1000,100)

plt.close('all')
fig = plt.figure(figsize = (4,4))
bins = 50
gs = gridspec.GridSpec(ncols = 5,nrows = 5)

ax_joint = fig.add_subplot(gs[1:5,0:4])
ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_y = fig.add_subplot(gs[1:5,4])

xmin, ymin = 0,0
xmax, ymax = 5.5,5.5

# Amplitude
x = A_pos_m
y = A_los_m


angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
ax_joint.text(xmax - 0.6,ymax - 0.55,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

ax_joint.scatter(x, y, alpha = 0.25, color = 'orangered')
ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = bins)
ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = bins/3)

# Set ax limits on marginals
ax_joint.set_xlim([xmin,xmax])
ax_joint.set_ylim([ymin,ymax])
ax_marg_y.set_ylim([ymin,ymax])
ax_marg_x.set_xlim([xmin,xmax])

pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
ax_joint.annotate('pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')

# Turn off ticks and labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
plt.setp(ax_marg_y.get_xticklabels(), visible=False)
plt.setp(ax_marg_x.get_yticklabels(), visible=False)
ax_marg_x.get_yaxis().set_visible(False)
ax_marg_y.get_xaxis().set_visible(False)

# horizontal hist
ax_marg_x.spines['top'].set_visible(False)
ax_marg_x.spines['right'].set_visible(False)
ax_marg_x.spines['left'].set_visible(False)

# vertical hist
ax_marg_y.spines['top'].set_visible(False)
ax_marg_y.spines['right'].set_visible(False)
ax_marg_y.spines['bottom'].set_visible(False)

# joint plot
ax_joint.spines['top'].set_visible(False)
ax_joint.spines['right'].set_visible(False)

# Set labels on joint
ax_joint.set_xlabel(r' $\overline{A}_{\rm{POS}}$ [km s$^{-1} $]', fontdict = font)
ax_joint.set_ylabel(r' $\overline{A}_{\rm{LOS}}$ [km s$^{-1} $]', fontdict = font)
ax_joint.tick_params(labelsize = fontsize)
ax_marg_x.tick_params(labelsize = fontsize)
ax_marg_y.tick_params(labelsize = fontsize)
ax_joint.tick_params(labelsize = fontsize)

plt.tight_layout()
plt.subplots_adjust(left = 0.14,
                    bottom = 0.11,
                    right = 1.,
                    top = 1.
)

filename = outdir + 'analys_AA.pdf'
fig.savefig(filename, dpi = 1000)


# POS PLOT
# ======
fontsize = 8.
xy = np.linspace(-0,1000,100)

plt.close('all')
fig = plt.figure(figsize = (4,4))
bins = 50
gs = gridspec.GridSpec(ncols = 5,nrows = 5)

ax_joint = fig.add_subplot(gs[1:5,0:4])
ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_y = fig.add_subplot(gs[1:5,4])

xmin, ymin = 110,0
xmax, ymax = 550,5.5

# POS
x = per_pos_m
y = A_pos_m


#angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
#ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
#ax_joint.text(xmax - 60,ymax - 0.55,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

ax_joint.scatter(x, y, alpha = 0.25, color = 'orangered')
ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = bins)
ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = bins)

# Set ax limits on marginals
ax_joint.set_xlim([xmin,xmax])
ax_joint.set_ylim([ymin,ymax])
ax_marg_y.set_ylim([ymin,ymax])
ax_marg_x.set_xlim([xmin,xmax])

#pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
#pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
#ax_joint.annotate('pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')

# Turn off ticks and labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
plt.setp(ax_marg_y.get_xticklabels(), visible=False)
plt.setp(ax_marg_x.get_yticklabels(), visible=False)
ax_marg_x.get_yaxis().set_visible(False)
ax_marg_y.get_xaxis().set_visible(False)

# horizontal hist
ax_marg_x.spines['top'].set_visible(False)
ax_marg_x.spines['right'].set_visible(False)
ax_marg_x.spines['left'].set_visible(False)

# vertical hist
ax_marg_y.spines['top'].set_visible(False)
ax_marg_y.spines['right'].set_visible(False)
ax_marg_y.spines['bottom'].set_visible(False)

# joint plot
ax_joint.spines['top'].set_visible(False)
ax_joint.spines['right'].set_visible(False)

# Set labels on joint
ax_joint.set_xlabel(r' $\overline{P}_{\rm{POS}}$ [s]', fontdict = font)
ax_joint.set_ylabel(r' $\overline{A}_{\rm{POS}}$ [km s$^{-1} $]', fontdict = font)
ax_joint.tick_params(labelsize = fontsize)
ax_marg_x.tick_params(labelsize = fontsize)
ax_marg_y.tick_params(labelsize = fontsize)
ax_joint.tick_params(labelsize = fontsize)

plt.tight_layout()
plt.subplots_adjust(left = 0.14,
                    bottom = 0.11,
                    right = 1.,
                    top = 1.
)

filename = outdir + 'analys_POS.pdf'
fig.savefig(filename, dpi = 1000)


# LOS PLOT
# ======
fontsize = 8.
xy = np.linspace(-0,1000,100)

plt.close('all')
fig = plt.figure(figsize = (4,4))
bins = 50
gs = gridspec.GridSpec(ncols = 5,nrows = 5)

ax_joint = fig.add_subplot(gs[1:5,0:4])
ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_y = fig.add_subplot(gs[1:5,4])

xmin, ymin = 110,0
xmax, ymax = 550,4

# LOS
x = per_los_m
y = A_los_m


#angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
#ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
#ax_joint.text(xmax - 60,ymax - 0.55,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

ax_joint.scatter(x, y, alpha = 0.25, color = 'orangered')
ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = bins)
ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = bins/2)

# Set ax limits on marginals
ax_joint.set_xlim([xmin,xmax])
ax_joint.set_ylim([ymin,ymax])
ax_marg_y.set_ylim([ymin,ymax])
ax_marg_x.set_xlim([xmin,xmax])

#pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
#pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
#ax_joint.annotate('pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')

# Turn off ticks and labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
plt.setp(ax_marg_y.get_xticklabels(), visible=False)
plt.setp(ax_marg_x.get_yticklabels(), visible=False)
ax_marg_x.get_yaxis().set_visible(False)
ax_marg_y.get_xaxis().set_visible(False)

# horizontal hist
ax_marg_x.spines['top'].set_visible(False)
ax_marg_x.spines['right'].set_visible(False)
ax_marg_x.spines['left'].set_visible(False)

# vertical hist
ax_marg_y.spines['top'].set_visible(False)
ax_marg_y.spines['right'].set_visible(False)
ax_marg_y.spines['bottom'].set_visible(False)

# joint plot
ax_joint.spines['top'].set_visible(False)
ax_joint.spines['right'].set_visible(False)

# Set labels on joint
ax_joint.set_xlabel(r' $\overline{P}_{\rm{LOS}}$ [s]', fontdict = font)
ax_joint.set_ylabel(r' $\overline{A}_{\rm{LOS}}$ [km s$^{-1} $]', fontdict = font)
ax_joint.tick_params(labelsize = fontsize)
ax_marg_x.tick_params(labelsize = fontsize)
ax_marg_y.tick_params(labelsize = fontsize)
ax_joint.tick_params(labelsize = fontsize)

plt.tight_layout()
plt.subplots_adjust(left = 0.14,
                    bottom = 0.11,
                    right = 1.,
                    top = 1.
)

filename = outdir + 'analys_LOS.pdf'
fig.savefig(filename, dpi = 1000)

plt.show()
stop()

