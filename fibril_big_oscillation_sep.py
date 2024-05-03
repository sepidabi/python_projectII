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
from scipy.optimize import curve_fit
from datetime import date
from datetime import datetime
import scipy.ndimage as spnd
import scipy.signal as sgnl
from mpl_toolkits import mplot3d
from skimage import exposure
from scipy.fftpack import fft2, ifft2


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

#functions
def test_func(x, a, b, c, d):
    return a * np.sin(b * x - c)+d

#boxcar smoothing
'''def smooth(y, w):
    N = y.shape[0]
    r = np.zeros(N)
    for i in range(N):
        if(i==0 or i==N-1):
            r[i] = y[i]
        elif(i>(w-1.)/2. and i<=N-(w+1.)/2.):
            r[i] = np.average(y[int(np.rint(i-w/2.)) : int(np.rint(i+w/2.-1))])
        else:
            r[i] = (y[i-1]+y[i]+y[i+1])/3.
    return r
'''
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
    return y[(window_len/2-1):-(window_len/2)]

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

fibdir = datadir+'big_cut/'

cc = 6

cut_file = (file_search(fibdir,'crispex*3950*.csav'))[cc]
cut = restore(fibdir+cut_file)
cube = cut.loop_slab[w_pos,:,:]

xx = cut.loop_slab.shape[2]

# to test the curve fitting
ti = np.arange(0, nt, 1.)

# plotting elements
tt = 30
ttick_pos = np.arange(0,(np.round((213)/tt)+1)*tt,tt)
ttick_lab = ttick_pos*cad/60

rr = 40/1.5
xtick_pos = np.arange(0,(np.round((xx)/rr)+1)*rr,rr)
xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))

# Extracting the inversion results
inv_res = sp.model(invdir+
                   file_search(invdir, "*atmosout*"+cut_file[-11:-5]+"*nc")[0])

temp_cube = inv_res.temp[0]
vlos_cube = inv_res.vlos[0]
vturb = inv_res.vturb[0]
ltau = inv_res.ltau[0,0,0,:]
ntau = ltau.size

for l in range(ntau):
    vlos_cube[:,:,l] = (vlos_cube[:,:,l]
                        - spnd.filters.gaussian_filter(inv_res.vlos[0,:,:,l], [65,65], mode = 'constant')).squeeze()*1e-5
    
tau_i = 23
    

# intensity oscillation PLOTs
plt.close('all')
fig = plt.figure(figsize = (4,4.6))
gs = gridspec.GridSpec(2,1, figure = fig)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])

fig2 = plt.figure(figsize = (8,4.45))
gs2 = gridspec.GridSpec(2,2, figure = fig2)
axx1 = plt.subplot(gs2[0,0])
axx2 = plt.subplot(gs2[1,0])
axx3 = plt.subplot(gs2[0,1])
axx4 = plt.subplot(gs2[1,1])

for oo in range(len(file_search(outdir, "cut"+cut_file[-11:-5]+"*txt"))):

    osc_fname = outdir+file_search(outdir, "cut"+cut_file[-11:-5]+"*txt")[oo]
    coord = np.loadtxt(osc_fname)
    amin, amax = int(np.min(coord[:,0])), int(np.max(coord[:,0]))+1
    tmin, tmax = int(np.min(coord[:,1])), int(np.max(coord[:,1]))+1
    trange = np.linspace(tmin, tmax, tmax-tmin+1)
    
    # interpolating the points in between the clicked fit
    y_interp = np.interp(trange, coord[:,1], coord[:,0])
    osc_fname_interp = outdir+"interp"+osc_fname[31:]
    np.savetxt(osc_fname_interp, y_interp, fmt='%3.8f', encoding=None)
    
    # track max intensity based on the interp click_fit
    y_imax = np.zeros(tmax-tmin+1)
    dist = 1. # width of the interval in which the max intensity is taken
    for i in range(tmax-tmin+1):
        y_imax[i] = np.argmax(cube[i, int(y_interp[i]-dist):int(y_interp[i]+dist)])+int(y_interp[i])-dist

    smth_y = smooth(y_imax, 8)
    osc_fname_smth = outdir+"smth"+osc_fname[31:]
    np.savetxt(osc_fname_smth, smth_y, fmt='%3.8f', encoding=None)

    # Calculation of the velocity of the oscillation in the plane of sky
    xp = trange*cad
    fxp = smth_y*res*725
    der = derivative(xp, fxp, dx = 1e-5)
    x, fx, f1x = der[0]-np.min(der[0]), der[1], der[2]

    # Intensity
    #------------
    im_cak = ax1.imshow(np.transpose(cube).squeeze()*1e9,
               cmap = 'gray', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               #vmin=2.e-9,
               #vmax = 3.5e-9,
               aspect = 0.6,
    )

    im_cak2 = ax2.imshow(gamma(np.transpose(cube).squeeze()*1e9,0.1)-9.,
               cmap = 'gray', origin = 'lower', #aspect=ratio
               interpolation = 'hermite',
               vmin=1.2,
               vmax = 1.9,
               aspect = 0.6,
    )
    ax2.plot(coord[:,1], coord[:,0],color = 'wheat', alpha = 0.5, label = 'Manual fit', linestyle = '--', linewidth = 1)
    if (oo ==7 or oo==5 or oo==12 or oo==6):
        ax2.plot(trange, smth_y, color = 'red', alpha = 0.5, label = r'smthd I$_{\rm{max}}$', linewidth = 2)
        if (oo==6):
            ax2.text(trange[0]+1, smth_y[0]+1.5, "1", color = 'red')
        if (oo==12):
            ax2.text(trange[0]+1, smth_y[0]+1.5, "2", color = 'red')
        if (oo==5):
            ax2.text(trange[0]+1, smth_y[0]+1.5, "3", color = 'red')
        if (oo==7):
            ax2.text(trange[0]+1, smth_y[0]+1.5, "4", color = 'red')

    elif (oo==6):
        print("nothing!")
    else:
        ax2.plot(coord[:,1], coord[:,0],color = 'wheat', alpha = 0.5, label = 'Manual fit', linestyle = '--', linewidth = 1)
        ax2.plot(trange, smth_y, color = 'orange', alpha = 0.5, label = r'smthd I$_{\rm{max}}$', linewidth = 2)


    # temp
    #--------
    tau = 23
    temp = np.transpose(temp_cube[:,:,tau]).squeeze()*1e-3
    temp_sharpened = exposure.rescale_intensity(sharpen(temp, sigma =[1,3], unsigma = [3,1]), out_range=(np.min(temp), np.max(temp)))
    im_temp = axx3.imshow(temp_sharpened,
               vmin = 3.7, vmax = 4.8,
               cmap = 'inferno', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               #vmin=2.e-9,
               #vmax = 3.5e-9,
               aspect = 0.6,
    )
    
    if (oo==6):
        print("nothing!2")
    else:
        axx3.plot(trange, smth_y, '--', color = 'white', alpha = 0.75, label = r'smthd I$_{\rm{max}}$', linewidth = 2.)


    # v_los
    #--------
    vlos = np.transpose(vlos_cube[:,:,tau])
    vlos_lim0 = np.round(np.max(np.abs(vlos)), decimals = 1)
    vlos_sharpened  = exposure.rescale_intensity(sharpen(vlos, sigma =[1,3], unsigma = [3,1]), out_range=(np.min(vlos), np.max(vlos)))
    im_vlos = axx4.imshow(vlos_sharpened,
               vmin = -vlos_lim0, vmax = vlos_lim0,
               cmap = 'bwr', origin = 'lower', #aspect=ratio
               # interpolation = 'bilinear',
               #vmin=2.e-9,
               #vmax = 3.5e-9,
               aspect = 0.6,
    )
    
    # temp
    #--------
    tau = 28
    temp = np.transpose(temp_cube[:,:,tau]).squeeze()*1e-3
    temp_sharpened = exposure.rescale_intensity(sharpen(temp, sigma =[1,3], unsigma = [3,1]), out_range=(np.min(temp), np.max(temp)))
    im_temp0 = axx1.imshow(temp_sharpened,
               vmin = 4., vmax =5,
               cmap = 'inferno', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               #vmin=2.e-9,
               #vmax = 3.5e-9,
               aspect = 0.6,
    )
    
    # v_los
    #--------
    vlos = np.transpose(vlos_cube[:,:,tau])
    vlos_lim0 = np.round(np.max(np.abs(vlos)), decimals = 1)
    vlos_sharpened  = exposure.rescale_intensity(sharpen(vlos, sigma =[1,3], unsigma = [3,1]), out_range=(np.min(vlos), np.max(vlos)))
    im_vlos0 = axx2.imshow(vlos_sharpened,
               vmin = -vlos_lim0, vmax = vlos_lim0,
               cmap = 'bwr', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               #vmin=2.e-9,
               #vmax = 3.5e-9,
               aspect = 0.6,
    )
    
    if (oo==6):
        print("nothing!2")
    else:
        axx4.plot(trange, smth_y, color = 'black', alpha = 0.5, label = r'smthd I$_{\rm{max}}$', linewidth = 1.)


    # smooth fit to the oscillations
    Imax_smth = np.zeros(tmax-tmin+1)
    temp_smth = np.zeros((tmax-tmin+1,ntau))
    vlos_smth = np.zeros((tmax-tmin+1,ntau))

    for i in range(tmax-tmin+1):
        Imax_smth[i] = cube[i,int(np.rint(smth_y[i]))]
        for j in range(ntau):
            temp_smth[i,j] = temp_cube[i,int(np.rint(smth_y[i])),j]
            vlos_smth[i,j] = vlos_cube[i,int(np.rint(smth_y[i])),j]


ax1.set_ylabel(r'length [arcsec]', fontdict = font)
ax1.set_xticks(ttick_pos)
ax1.set_xticklabels([])
ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
ax1.set_xlabel('')
ax1.set_yticks(xtick_pos)
ax1.set_yticklabels(xtick_lab, fontdict = font)
ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
ax1.set_xlim(0,nt-1)    
ax1.set_ylim(0,201)

ax2.set_ylabel(r'length [arcsec]', fontdict = font)
ax2.set_xticks(ttick_pos)
ax2.set_xticklabels(ttick_lab, fontdict = font)
ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
ax2.set_xlabel('t [min]', fontdict = font)
ax2.set_yticks(xtick_pos)
ax2.set_yticklabels(xtick_lab, fontdict = font)
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.set_xlim(0,nt-1)
ax2.set_ylim(0,201)

axx1.set_ylabel(r'length [arcsec]', fontdict = font)
axx1.set_xticks(ttick_pos)
axx1.set_xticklabels([])
axx1.xaxis.set_minor_locator(AutoMinorLocator(4))
axx1.set_xlabel('')
axx1.set_yticks(xtick_pos)
axx1.set_yticklabels(xtick_lab, fontdict = font)
axx1.yaxis.set_minor_locator(AutoMinorLocator(5))
axx1.set_xlim(0,nt-1)    
axx1.set_ylim(0,201)
axx1.set_title(r'log($\tau_{500}$) = -3.5', size = 9)

axx2.set_ylabel(r'length [arcsec]', fontdict = font)
axx2.set_xticks(ttick_pos)
axx2.set_xticklabels(ttick_lab, fontdict = font)
axx2.xaxis.set_minor_locator(AutoMinorLocator(4))
axx2.set_xlabel('t [min]', fontdict = font)
axx2.set_yticks(xtick_pos)
axx2.set_yticklabels(xtick_lab, fontdict = font)
axx2.yaxis.set_minor_locator(AutoMinorLocator(5))
axx2.set_xlim(0,nt-1)
axx2.set_ylim(0,201)

axx3.set_ylabel('')
axx3.set_xticks(ttick_pos)
axx3.set_xticklabels([])
axx3.xaxis.set_minor_locator(AutoMinorLocator(4))
axx3.set_xlabel('')
axx3.set_yticks(xtick_pos)
axx3.set_yticklabels([])
axx3.yaxis.set_minor_locator(AutoMinorLocator(5))
axx3.set_xlim(0,nt-1)    
axx3.set_ylim(0,201)
axx3.set_title(r'log($\tau_{500}$) = -4.3', size = 9)

axx4.set_ylabel('')
axx4.set_xticks(ttick_pos)
axx4.set_xticklabels(ttick_lab, fontdict = font)
axx4.xaxis.set_minor_locator(AutoMinorLocator(4))
axx4.set_xlabel('t [min]', fontdict = font)
axx4.set_yticks(xtick_pos)
axx4.set_yticklabels([])
axx4.yaxis.set_minor_locator(AutoMinorLocator(5))
axx4.set_xlim(0,nt-1)
axx4.set_ylim(0,201)


axins1 = inset_axes(axx1,
                    width="3%",  
                    height="90%",
                    loc='center left',
                    bbox_to_anchor=(1.01, 0., 1, 1),
                    bbox_transform=axx1.transAxes,
                    borderpad=0,
)


axins2 = inset_axes(axx2,
                    width="3%",  
                    height="90%",
                    loc='center left',
                    bbox_to_anchor=(1.01, 0., 1, 1),
                    bbox_transform=axx2.transAxes,
                    borderpad=0,
)

axins3 = inset_axes(axx3,
                    width="3%",  
                    height="90%",
                    loc='center left',
                    bbox_to_anchor=(1.01, 0., 1, 1),
                    bbox_transform=axx3.transAxes,
                    borderpad=0,
)


axins4 = inset_axes(axx4,
                    width="3%",  
                    height="90%",
                    loc='center left',
                    bbox_to_anchor=(1.01, 0., 1, 1),
                    bbox_transform=axx4.transAxes,
                    borderpad=0,
)

cbar1 = fig2.colorbar(im_temp0, cax=axins1)
cbar2 = fig2.colorbar(im_vlos0, cax=axins2)
cbar3 = fig2.colorbar(im_temp, cax=axins3)
cbar4 = fig2.colorbar(im_vlos, cax=axins4)

cbar4.set_label(r'$v_{\rm{LOS}}$ [km s$^{-1}$]', size = 8.)
cbar3.set_label(r'$T$ [kK]', size = 8.)

cbar1.ax.tick_params(axis='y', labelsize = 8.)
cbar2.ax.tick_params(axis='y', labelsize = 8.)
cbar3.ax.tick_params(axis='y', labelsize = 8.)
cbar4.ax.tick_params(axis='y', labelsize = 8.)

gs2.update(left=0.045,
          right=0.925,
          wspace=0.15,
          bottom=0.075,
          top=0.965,
          hspace = 0.005,
)

gs.update(left=0.1,
          right=0.98,
          wspace=0.0,
          bottom=0.09,
          top=0.995,
          hspace = 0.05,
)

plt.show()
#stop()
filename = outdir + 'oscillation_curvefit_' +cut_file[-11:-5]+'_non_int.pdf'
fig.savefig(filename, quality = 100)
print 'file saved to: ' + filename

filename2 = outdir + 'oscillation_curvefit_' +cut_file[-11:-5]+'_non_inv.pdf'
fig2.savefig(filename2, quality = 100)
print 'file saved to: ' + filename2


