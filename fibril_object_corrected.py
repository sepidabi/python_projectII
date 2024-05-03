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
dphi = np.load('dphi.npy') # np.zeros(0)

ncc_min = np.zeros(0)
ncc_max =  np.zeros(0)

mid_coord_x = np.zeros(0)
mid_coord_y = np.zeros(0)
cor_rate = np.zeros(0)
theta = np.zeros(0)
per_dom = np.zeros(0)
ratio = np.zeros(0)
dist = np.zeros(0)

# Extract the oscillation objects

# User contribution to the specify the oscillation
print('-----------------------------------------------------------------------------------------------------------------------')
print('index                                  |dP|           Ppos_m       Ppos_dom       Plos_m       Plos_dom       P_corr')
print('-----------------------------------------------------------------------------------------------------------------------')
ind = 0

for iii in range(len(file_search(objdir, "cut*osc0.obj"))):
    cut_indx = file_search(objdir, "cut*osc0.obj")[iii][3:12]
    
    # coords of the centere of oscillation
    cut_file = (file_search(fibdir,'crispex*3950*'+cut_indx+'*.csav'))[0]
    cut = restore(fibdir+cut_file)
    cube = cut.loop_slab[w_pos,:,:]

    for jjj in range(len(file_search(objdir, 'cut'+cut_indx+'*osc*.obj'))):
        
        per_los_one = np.zeros(0)
        per_pos_one = np.zeros(0)
    
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
                per_los_one = np.append(per_los_one, 2*np.abs(osc.los_xmax-osc.los_xmin))
            if (len(pe_l)==3):
                per_los_one = np.append(per_los_one, (pe_l[-1] - pe_l[0]))
            if (len(pe_l)>3):
                for l in range (len(pe_l)-3):
                    per_los_one = np.append(per_los_one, (pe_l[l+3]-pe_l[l+1]+pe_l[l+2]-pe_l[l])/2.)
            #if (len(pe_l)==1):
             #   per_los_one = pe_l*2
            per_los_dom = np.max(per_los_one)    
            per_los = np.append(per_los, per_los_one)
            #per_los_dom = np.append(per_los_dom, per_los_one[np.argmin(np.abs(per_los_one-osc.los_ncc_xmax[0]))])
            per_los_m = np.append(per_los_m, np.mean(per_los_one))
            A_los_m = np.append(A_los_m, np.mean([np.mean(np.abs(osc.los_ymax)), np.mean(np.abs(osc.los_ymin))]))
            
            
            # POS period
            pe_p = np.sort(np.append(osc.pos_xmax, osc.pos_xmin))
            
            if (len(pe_p)==2):
                per_pos_one = np.append(per_pos_one, 2*np.abs(osc.pos_xmax-osc.pos_xmin))
            if (len(pe_p)==3):
                per_pos_one = np.append(per_pos_one, (pe_p[-1] - pe_p[0]))
            if (len(pe_p)>3):
                for m in range (len(pe_p)-3):
                    per_pos_one = np.append(per_pos_one, (pe_p[m+3]-pe_p[m+1]+pe_p[m+2]-pe_p[m])/2.)
                    
            per_pos_dom = np.max(per_pos_one)    
            per_pos = np.append(per_pos, per_pos_one) 
            #per_pos_dom = np.append(per_pos_dom, per_pos_one[np.argmin(np.abs(per_pos_one-osc.pos_ncc_xmax[0]))])
            per_pos_m = np.append(per_pos_m, np.mean(per_pos_one))
            A_pos_m = np.append(A_pos_m, np.mean([np.mean(np.abs(osc.pos_ymax)), np.mean(np.abs(osc.pos_ymin))]))
            
            
            ssize = 30 # plotting marker size
            
            per_one_dom = np.max([np.max(per_pos_one), np.max(per_los_one)])
            per_dom = np.append(per_dom, np.mean([np.mean(per_pos_one), np.mean(per_los_one)]))
            lag = np.arange(-len(vlos_f_norm)+1, len(vpos_f_norm))*8.

            Ppos_m = np.mean(per_pos_one)
            Ppos_dom = np.max(per_pos_one)
            Plos_m = np.mean(per_los_one)
            Plos_dom = np.max(per_los_one)
            dp = np.round(np.abs(Ppos_m - Plos_m), decimals = 1)
            dp_dom = np.round(np.abs(per_pos_dom - per_los_dom), decimals = 1)

            corr_thre = 0.5 # the cross-correlation threshold above which the oscillations are assumed to be (anti)correlated

            if (np.max(ncc)>corr_thre and np.abs(lag[np.argmax(ncc)])<=per_one_dom and np.max(ncc)>np.abs(np.min(ncc))):
                #test = np.abs(np.diff(lag[extrem(ncc, np.greater)]))[np.argmin(np.abs(np.diff(lag[extrem(ncc, np.greater)])) - per_one_dom)],#]),

                      #np.abs(np.diff(lag[extrem(ncc, np.greater)])),
                      #np.abs(lag[
                      #np.where(np.diff(lag[extrem(ncc, np.greater)]) - per_one_dom<0)
                      #)
                P_cor = np.abs(lag[np.argmax(ncc)])
                ind = ind+1
                print(
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
                )
                #print('Plos = ', per_pos_one)
                #print('Ppos = ', per_los_one)
                #str(np.round(dp*100/np.mean([Ppos_m, Plos_m]))))# + '%Ppos' + ' = ' + str(np.round(dp*100/Plos_m)) + '%Plos')
                #', Pm = ' + str(np.round(np.mean([np.mean(per_los_one),np.mean(per_pos_one)]))))                
                #print('max', str(iii)+'-'+str(jjj), 2*lag[np.argmax(ncc)], per_one_dom, dphi_one)
                #ucheck = input('1 or 0? ')
                #if ucheck==1:
                 #   dphi_one = lag[np.argmax(ncc)]/per_one_dom
                  #  #if (dphi_one>1):
                   #     #dphi_one = 2 - dphi_one
                   # dphi = np.append(dphi, dphi_one)
                   # #if (dphi_one<0):
                   # #   dphi_one = 2+dphi_one
                #if ucheck==0:
                 #   dphi = np.append(dphi, -5)
                    
            elif (np.min(ncc)<-corr_thre and np.abs(lag[np.argmin(ncc)])<=per_one_dom and np.max(ncc)<np.abs(np.min(ncc))):
                #test = np.abs(np.diff(lag[extrem(ncc, np.less)]))[np.argmin(np.abs(np.diff(lag[extrem(ncc, np.less)])) - per_one_dom)],#]),
                ind = ind+1

                P_cor = np.abs(lag[np.argmin(ncc)])
                print(
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
                )
                #str(np.round(dp*100/np.mean([Ppos_m, Plos_m]))))# + '%Ppos' + ' = ' + str(np.round(dp*100/Plos_m)) + '%Plos')
                #print('Plos = ', per_pos_one)
                #print('Ppos = ', per_los_one)

                #ucheck = input('1 or 0? ')
                #if ucheck==1:
                 #   dphi_one = lag[np.argmin(ncc)]/per_one_dom
                  #  #if (dphi_one<0):
                   # #   dphi_one = 2+dphi_one
                   # dphi = np.append(dphi, dphi_one)
                #if ucheck==0:
                 #   dphi = np.append(dphi, -5)

            #else:
             #   dphi = np.append(dphi, -5)
        else:
            #dphi = np.append(dphi,-10)
            per_dom = np.append(0, per_dom)

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
fontsize = 10
xy = np.linspace(-0,1000,100)

plt.close('all')
fig = plt.figure(figsize = (8,8))
bins = 50
gs = gridspec.GridSpec(ncols = 9,nrows = 9)

ax_joint = fig.add_subplot(gs[1:5,0:4])
ax_joint2 = fig.add_subplot(gs[1:5,4:8])
ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_x2 = fig.add_subplot(gs[0,4:8])
ax_joint3 = fig.add_subplot(gs[5:9,0:4])
ax_joint4 = fig.add_subplot(gs[5:9,4:8])
ax_marg_y = fig.add_subplot(gs[1:5,8])
ax_marg_y2 = fig.add_subplot(gs[5:9,8])

xmin, ymin = 0, 0
xmax, ymax = 5.5, 5.5

# Period
x = per_pos_m.squeeze()*1e-2
y = per_los_m.squeeze()*1e-2

#xmin, xmax = 1.19, 4.50
#ymin, ymax = 1.19, 4.50


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
ax_joint.annotate('a) pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')


# Amplitudes
x = A_los_m
y = A_pos_m

#xmin, xmax = 0, 5.5
#ymin, ymax = 0, 5.5


angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
ax_joint4.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
ax_joint4.text(xmax - 1.25,ymax - 0.5,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

ax_joint4.scatter(x, y, alpha = 0.25, color = 'orangered')
ax_marg_x2.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = 20)
ax_marg_y2.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = 60)

# Set ax limits on marginals
ax_joint4.set_xlim([xmin,xmax])
ax_joint4.set_ylim([ymin,ymax])
ax_marg_y2.set_ylim([ymin,ymax])
ax_marg_x2.set_xlim([xmin,xmax])

pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
ax_joint4.annotate('d) pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')

# LOS
x = per_pos_m.squeeze()*1e-2
y = A_pos_m

#xmin, xmax = 1.19, 4.50
#ymin, ymax = 0, 5.5


angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
ax_joint3.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
ax_joint3.text(xmax - 1.25,ymax - 0.5,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

ax_joint3.scatter(x, y, alpha = 0.25, color = 'orangered')
ax_joint3.set_xlim([xmin,xmax])
ax_joint3.set_ylim([ymin,ymax])

pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
ax_joint3.annotate('c) pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')


# POS
x = A_los_m
y = per_los_m.squeeze()*1e-2

#xmin, xmax = 0, 5.5
#ymin, ymax = 1.19, 4.50


angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
ax_joint2.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
ax_joint2.text(xmax - 1.25,ymax - 0.5,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

ax_joint2.scatter(x, y, alpha = 0.25, color = 'orangered')
ax_joint2.set_xlim([xmin,xmax])
ax_joint2.set_ylim([ymin,ymax])

pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
ax_joint2.annotate('b) pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')

# Turn off ticks and labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
plt.setp(ax_marg_y.get_xticklabels(), visible=False)
plt.setp(ax_marg_x.get_yticklabels(), visible=False)
ax_marg_x.get_yaxis().set_visible(False)
ax_marg_y.get_xaxis().set_visible(False)
plt.setp(ax_marg_x2.get_xticklabels(), visible=False)
plt.setp(ax_marg_y2.get_yticklabels(), visible=False)
plt.setp(ax_marg_y2.get_xticklabels(), visible=False)
plt.setp(ax_marg_x2.get_yticklabels(), visible=False)
ax_marg_x2.get_yaxis().set_visible(False)
ax_marg_y2.get_xaxis().set_visible(False)

# horizontal hist
ax_marg_x.spines['top'].set_visible(False)
ax_marg_x.spines['right'].set_visible(False)
ax_marg_x.spines['left'].set_visible(False)
ax_marg_x2.spines['top'].set_visible(False)
ax_marg_x2.spines['right'].set_visible(False)
ax_marg_x2.spines['left'].set_visible(False)
# vertical hist
ax_marg_y.spines['top'].set_visible(False)
ax_marg_y.spines['right'].set_visible(False)
ax_marg_y.spines['bottom'].set_visible(False)
ax_marg_y2.spines['top'].set_visible(False)
ax_marg_y2.spines['right'].set_visible(False)
ax_marg_y2.spines['bottom'].set_visible(False)

# joint plot
ax_joint.spines['top'].set_visible(False)
ax_joint2.spines['top'].set_visible(False)
ax_joint2.spines['right'].set_visible(False)
ax_joint4.spines['right'].set_visible(False)


# Set labels on joint
#ax_joint2.set_xlabel(r' $\overline{P}_{\rm{POS}}$ [s]', fontdict = font)
ax_joint.set_ylabel(r' $\overline{P}_{\rm{LOS}}$ [hs]', fontdict = font)
ax_joint.tick_params(labelsize = fontsize)
ax_marg_x.tick_params(labelsize = fontsize)
ax_marg_y.tick_params(labelsize = fontsize)
ax_joint4.set_xlabel(r' $\overline{A}_{\rm{LOS}}$ [km s$^{-1}$]', fontdict = font)
ax_joint3.set_ylabel(r' $\overline{A}_{\rm{POS}}$ [km s$^{-1}$]', fontdict = font)
ax_joint3.set_xlabel(r' $\overline{P}_{\rm{POS}}$ [hs]', fontdict = font)
ax_joint4.set_yticklabels([])
ax_joint2.set_yticklabels([])
ax_joint2.set_xticklabels([])
ax_joint.set_xticklabels([])
#ax_joint3.set_xticklabels([])
ax_joint.tick_params(labelsize = fontsize)
ax_joint2.tick_params(labelsize = fontsize)
ax_joint3.tick_params(labelsize = fontsize)
ax_joint4.tick_params(labelsize = fontsize)
ax_marg_x2.tick_params(labelsize = fontsize)
ax_marg_y2.tick_params(labelsize = fontsize)

plt.tight_layout()
plt.subplots_adjust(left = 0.06,
                    bottom = 0.06,
                    right = 1.,
                    top = 1.,
                    wspace = 0.2,
                    hspace = 0.2
)
plt.show()
#stop()
filename = outdir + 'analys_tot.pdf'
#fig.savefig(filename, dpi = 1000)
plt.close(fig)


##############
# COLORMAP creation
##############
red = (178/255.,34/255.,34/255.)
orange = (1, 128/255., 0)
yellow = (1, 1, 0)
colors = [red, orange, yellow, orange, red]#, orange, yellow, orange, red]
cmap_name = 'phi'
nbin = 100
cmap_phi = LinearSegmentedColormap.from_list(cmap_name, colors, N=nbin)


# D-Phi histogram
#==========
plt.close('all')
fifi = plt.figure(figsize = (4,3))
#gs = gridspec.GridSpec(ncols = 1,nrows = 1, figure = fifi)
#axfi = fig.add_subplot(gs[0,0])
axfi = plt.subplot(111)#, projection='polar')
n, bins, patches  = axfi.hist(dphi[(dphi>=-1)], bins = 14, edgecolor='white', linewidth=1., color = 'orangered', alpha = 0.5, range = [-1,1])
bin_centers = 0.5 * (bins[:-1] + bins[1:])
axfi.tick_params(axis = 'both', labelsize = 8.)
axfi.set_xlabel(r'$\Delta \phi$ [rad]', fontdict = font)
axfi.set_ylabel(r'Number of oscillation cases', fontdict = font)
axfi.set_xticks([-1,-.75, -.5, -.25, 0,0.25,0.5,0.75,1])
axfi.set_xticklabels([r'-$\pi$', r'-3$\pi$/4', r'-$\pi$/2', r'-$\pi$/4', r'0', r'$\pi$/4', r'$\pi$/2', r'3$\pi$/4', r'$\pi$'], fontdict = font)
axfi.set_yticks([0,3,6,9,12,15])
axfi.set_ylim(0,16)
axfi.yaxis.set_minor_locator(AutoMinorLocator(3))

axins = inset_axes(axfi, width=1.2, height=.91,
                    bbox_to_anchor=(0.62, 0.57),
                    bbox_transform=axfi.transAxes, loc=3, borderpad=0)

nn, bbins, ppatches = axins.hist(np.abs(dphi[(dphi>=-1)]), bins = 7, range = [0,1], edgecolor = 'white', linewidth = 1., color = 'orangered', alpha = 0.5)
axins.tick_params(axis = 'both', labelsize = 8.)
axins.set_xlabel(r'|$\Delta \phi$| [rad]', fontdict = font)
#axins.set_ylabel(r'Number of oscillation cases', fontdict = font)
axins.set_yticks([0,5,10,15])
axins.set_xticks([0,0.5,1])
axins.set_xticklabels([r'0', r'$\pi$/2', r'$\pi$'], fontdict = font)
#axins.yaxis.set_minor_locator(AutoMinorLocator(2))

plt.tight_layout()
filename = outdir + 'analys_dphi.pdf'
plt.show()
#fifi.savefig(filename, dpi = 1000)
#stop()
plt.close(fifi)

#axfi.set_xlim([-bin_centers[0]/2,1])

# This is  the colormap I'd like to use.
#cm = plt.cm.get_cmap(cmap_phi)

#col = bin_centers - min(bin_centers)
#col /= max(col)

#for c, p in zip(col, patches):
  #  plt.setp(p, 'facecolor', cm(c))


# SCATTER PLOT
#==========
cor_rate = dphi
ssize = 300
# RoI coords
x1, x2, y1, y2 = 425, 900, 400, 1110
xx, yy = x2-x1, y2-y1
plt.close("all")
map_cak = plt.figure(figsize = (4,6))
gs_map = gridspec.GridSpec(ncols = 1,nrows = 1, figure = map_cak)
ax_map =  plt.subplot(gs_map[0,0])
ax_map.imshow(cak_core[y1:y2, x1:x2], cmap = 'gray', origin = 'lower')
ax_map.contourf(np.abs(CPtot[y1:y2, x1:x2]).squeeze()*1e2, cmap = 'Blues', alpha = 0.5, levels = [2,8])
contour = ax_map.contour(np.abs(CPtot[y1:y2, x1:x2]).squeeze()*1e2, colors = ['white', 'white'], alpha = 0.75, levels = [2,8], linewidths = [0.5,0.5])
for i in range (len(mid_coord_x)):
    if (dphi[i]>=-1):
        sc = ax_map.scatter(x = mid_coord_x[i]-x1, y = mid_coord_y[i]-y1, c = np.abs(cor_rate[i]), marker=(3, 0, theta[i]), cmap = cmap_phi, s = np.power(per_dom[i]/np.max(per_dom), 2)*ssize, alpha = 0.5, vmin =0, vmax = 1)
        dist_i = np.sqrt(contour.find_nearest_contour(mid_coord_x[i]-x1, mid_coord_y[i]-y1)[5])
        dist = np.append(dist, dist_i)
        ratio = np.append(ratio,dist_i/per_dom[i])
    elif (dphi[i]==-5):
        ax_map.scatter(x = mid_coord_x[i]-x1, y = mid_coord_y[i]-y1, marker=(3, 0, theta[i]) , s = np.power(per_dom[i]/np.max(per_dom), 2)*ssize, alpha = 0.5, color = 'gray')
        #dist_i = np.sqrt(contour.find_nearest_contour(mid_coord_x[i]-x1, mid_coord_y[i]-y1)[5])
        #dist = np.append(dist, dist_i)
        #ratio = np.append(ratio,dist_i/per_dom[i])
    else:
        ax_map.scatter(x = mid_coord_x[i]-x1, y = mid_coord_y[i]-y1, marker='x' , s = 25, alpha = 0.5, color = 'white')
#ax_map.scatter(x=100,  y = 100, marker = (3,0,0), s = np.power(140/np.max(per_dom), 2)*ssize)
hands =np.array( [140,240,340,440.])
sizes = np.sqrt(np.power(hands/np.max(per_dom), 2)*ssize)
#sizes = np.power(sizes,0.5)
legend_elements = [Line2D([0], [0], marker='^', color='w', label='140',markerfacecolor='gray', markersize=sizes[0], linestyle='None', markeredgecolor = 'darkgray'), #linewidth = 1.),
                   Line2D([], [], marker='^', color='w', label='240',markerfacecolor='gray', markersize=sizes[1], linestyle='None', markeredgecolor = 'darkgray'),
                   Line2D([], [], marker='^', color='w', label='340',markerfacecolor='gray', markersize=sizes[2], linestyle='None', markeredgecolor = 'darkgray'),
                   Line2D([], [], marker='^', color='w', label='440',markerfacecolor='gray', markersize=sizes[3], linestyle='None', markeredgecolor = 'darkgray'),
                   Line2D([], [], marker='^', color='w', label='',markerfacecolor='gray', markersize=0, linestyle='None'),
                   Line2D([], [], marker='^', color='w', label='',markerfacecolor='gray', markersize=0, linestyle='None'),
]

lghd = ax_map.legend(handles=legend_elements, loc='lower right', framealpha = 0., title = r'P$_{\rm{mean}}$ [s]',borderpad = 1.2, prop={'size' : 8.}, labelspacing = 1.2)#, title_fontsize = 8.,markerscale = 0.5,borderpad = 0.15)
#lghd.legendHandles[0]._legmarker.set_markersize(15)
# colorbar
axin = inset_axes(ax_map,
                  width="25%",  # width = 10% of parent_bbox width
                  height="1.5%",  # height : 50%
                  loc='lower right',
)

plt.setp(lghd.get_texts(), color='w')
plt.setp(lghd.get_title(), color='w')

cb = plt.colorbar(sc, cax=axin, orientation="horizontal", ticks = [-1, -.5, 0, 0.5,  1])
#cb.ax.set_xticklabels(cb_tick_labels, fontdict = font)
cb.ax.set_xticklabels([ r'0',  r'$\pi / 2$', r'$\pi$'], color = 'white')
cb.ax.xaxis.tick_top()
cb.ax.xaxis.set_tick_params(size = 8)
cb.ax.tick_params(axis = 'x', labelsize = 8., length = 2, color = 'white')
cb.ax.set_title(r'$|\Delta \phi$| [rad]', fontdict = font, pad = 18, color = 'white')
cb.outline.set_edgecolor('white')

sc_fact = 266.7 # axis scaling factor
ytick_pos = np.arange(0,(np.round(yy/sc_fact))*sc_fact,sc_fact/2.)
ytick_lab = np.round(ytick_pos*res).astype(int)
xtick_pos = np.arange(0,(np.round(xx/sc_fact))*sc_fact,sc_fact/2.)
xtick_lab = np.round(xtick_pos*res).astype(int)

ax_map.set_xticks(xtick_pos)
ax_map.set_xlim(0,xx)
ax_map.set_ylim(0,yy)
ax_map.set_yticks(ytick_pos)
ax_map.set_ylabel(r'y [arcsec]', fontdict=font)
ax_map.set_yticklabels(ytick_lab, fontdict = font)
ax_map.xaxis.set_minor_locator(AutoMinorLocator(5))
ax_map.yaxis.set_minor_locator(AutoMinorLocator(5))
ax_map.tick_params(which='minor', length=2)
ax_map.set_xticks(xtick_pos)
ax_map.set_xticklabels(xtick_lab, fontdict = font)
ax_map.set_xlabel(r'x [arcsec]', fontdict=font)

gs_map.update(left=0.115,
              right=0.975,
              wspace=0.1,
              bottom=0.02,
              top=0.98,
              hspace = 0.0
)
filename = outdir + 'analys_fibril_coords.pdf'
#map_cak.savefig(filename, dpi = 1000)
print('file saved to: '+ filename)
plt.show()

# fibril idstance and period correlation
plt.figure()
plt.plot(dist/ratio, dist, 'r.')

z = np.polyfit(dist/ratio, dist, 1)
p = np.poly1d(z)
plt.plot(dist/ratio, p(dist/ratio))

plt.xlabel('P [s]')
plt.ylabel('D [pixel]')
plt.title('corr. fibril d. from mag. patch \nvs. avg. osc. period')

slope = (p(dist/ratio)[-1]-p(dist/ratio)[0])/(dist[-1]/ratio[-1]-dist[0]/ratio[0])
plt.text(270,250,'D = '+str(np.round(slope, decimals = 1))+'*P')

plt.show()

stop()
plt.close(map_cak)




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

xmin, ymin = 80, 80
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
ax_joint.set_xlabel(r' $\overline{P}_{\rm{POS}}$ s]', fontdict = font)
ax_joint.set_ylabel(r' $\overline{P}_{\rm{LOS}}$ s]', fontdict = font)
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
#fig.savefig(filename, dpi = 1000)
stop()

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
#fig.savefig(filename, dpi = 1000)
stop()

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

xmin, ymin = 100,0
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
#fig.savefig(filename, dpi = 1000)
stop()

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

xmin, ymin = 80,0
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
#fig.savefig(filename, dpi = 1000)

plt.show()
        
