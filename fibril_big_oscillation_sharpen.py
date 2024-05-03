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
outdir_old = '/home/seki2695/OUTPUT/project2_old/'
invdir = '/home/seki2695/INV/stic/II/slabs_tseries/'

# Difinitions
res = 0.0375 # CHROMIS pixel size in arcse
cad = 8
fr_select = [29,75,135,167] # frames to display
w_pos = 13 # wavelength position to display (K2V, K3, K2R) = (10,13,16)

#functions

# smoothing functin
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


fibdir = datadir+'big_cut/'

cc = 6 # the cut index
cut_file = (file_search(fibdir,'crispex*3950*.csav'))[cc]
cut = restore(fibdir+cut_file)
cube_untransposed = cut.loop_slab
cube = np.transpose(cube_untransposed).squeeze()*1e9
nt = cube.shape[1]
nx = cube.shape[0]
nw = cube.shape[2]
cube_core = cube[:,:,w_pos]
cube_int = np.mean(cube[:,:,w_pos-3:w_pos+3],axis = 2)
# IMAGE ENHANCEMENT
#-----------------------------
wg = 0.15
cube_bg = wg*np.mean(cube[:,:,:w_pos-6]+cube[:,:,w_pos+8:],axis=2)
    
cube_trunc = cube_int +30*(cube_int - cube_bg) # removing the background intensity

cube_blur = spnd.gaussian_filter(cube_trunc, [1,3])

cube_filter_blur = spnd.gaussian_filter(cube_trunc, [3,1])

cube_med = exposure.rescale_intensity(sgnl.medfilt2d(cube_trunc,kernel_size = [1,3]), out_range=(0, 1.))

cube_sharp = exposure.rescale_intensity(sharpen(cube_trunc, sigma =[1,3], unsigma = [3,1]), out_range=(0, 1.))

cube_sharp_gamma = exposure.rescale_intensity(exposure.adjust_gamma(cube_sharp,0.1), out_range=(0, 1.))#exposure.adjust_log(np.abs(np.min(cube_sharp))+cube_sharp)

# Plots
plt.close('all')
fig = plt.figure(figsize = (10,6.5))#4.6))
gs = gridspec.GridSpec(3,3, figure = fig)
ax00 = plt.subplot(gs[0,0])
ax10 = plt.subplot(gs[1,0])
ax20 = plt.subplot(gs[2,0])
ax01 = plt.subplot(gs[0,1])
ax11 = plt.subplot(gs[1,1])
ax21 = plt.subplot(gs[2,1])
ax02 = plt.subplot(gs[0,2])
ax12 = plt.subplot(gs[1,2])
ax22 = plt.subplot(gs[2,2])

# PLOTTING results
ax00.imshow(cube_core,
            cmap = 'gray', origin = 'lower', aspect = 0.6,
            interpolation = 'bilinear',
#          vmin=0.05,
#          vmax = 0.9,
)

ax01.imshow(cube_int,
            cmap = 'gray', origin = 'lower', aspect = 0.6,
            interpolation = 'bilinear',
            #vmin=0.05,
            #vmax = 0.9,                   
)

ax02.imshow(cube_bg,
            cmap = 'gray', origin = 'lower', aspect = 0.6,
            interpolation = 'bilinear',
            #vmin=0.05,
            #vmax = 0.9,                   
)

ax10.imshow(cube_trunc,
            cmap = 'gray', origin = 'lower', aspect = 0.6,
            #interpolation = 'hermite',
            #vmin=0.,
            #vmax = 1.,
)

ax11.imshow(cube_blur,
            cmap = 'gray', origin = 'lower', aspect = 0.6,
            #interpolation = 'hermite',
            #vmin=0.89,
            #vmax = 0.98,
)

ax12.imshow(cube_filter_blur,
                     cmap = 'gray', origin = 'lower', aspect = 0.6,
            #interpolation = 'hermite',
            #vmin=0.89,
            #vmax = 0.98,
)

ax20.imshow(cube_sharp,
            cmap = 'gray', origin = 'lower', aspect = 0.6,
            #interpolation = 'hermite',
            #vmin=0.89,
            #vmax = 0.98,
)

ax21.imshow(cube_sharp_gamma,
            cmap = 'gray', origin = 'lower', aspect = 0.6,
            #interpolation = 'bilinear',
            vmin = 0.89,
            vmax = 0.98,
)

ax22.imshow(cube_sharp_gamma,
            cmap = 'gray', origin = 'lower', aspect = 0.6,
            #interpolation = 'bilinear',
            vmin = 0.89,
            vmax = 0.98,
)


xx = cut.loop_slab.shape[2]

# plotting elements
tt = 30
ttick_pos = np.arange(0,(np.round((nt)/tt)+1)*tt,tt)
ttick_lab = ttick_pos*cad/60

rr = 40/1.5
xtick_pos = np.arange(0,(np.round((xx)/rr)+1)*rr,rr)
xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))

ax00.set_ylabel(r'length [arcsec]', fontdict = font)
ax00.set_xticks(ttick_pos)
ax00.set_xticklabels([])
ax00.xaxis.set_minor_locator(AutoMinorLocator(4))
ax00.set_xlabel(r'')
ax00.set_yticks(xtick_pos)
ax00.set_yticklabels(xtick_lab, fontdict = font)
ax00.yaxis.set_minor_locator(AutoMinorLocator(5))
ax00.set_xlim(0,nt-1)    
ax00.set_ylim(0,nx)
ax00.set_title(r'A: Ca II K linecore with former oscillation paths', fontdict=font)

ax01.set_ylabel('')
ax01.set_xticks(ttick_pos)
ax01.set_xticklabels([])
ax01.xaxis.set_minor_locator(AutoMinorLocator(4))
ax01.set_xlabel('')
ax01.set_yticks(xtick_pos)
ax01.set_yticklabels([])
ax01.yaxis.set_minor_locator(AutoMinorLocator(5))
ax01.set_xlim(0,nt-1)    
ax01.set_ylim(0,nx)
ax01.set_title(r'B: Ca II K wavelength averaged around the K2 peaks', fontdict=font)

ax02.set_ylabel('')
ax02.set_xticks(ttick_pos)
ax02.set_xticklabels([])
ax02.xaxis.set_minor_locator(AutoMinorLocator(4))
ax02.set_xlabel('')
ax02.set_yticks(xtick_pos)
ax02.set_yticklabels([])
ax02.yaxis.set_minor_locator(AutoMinorLocator(5))
ax02.set_xlim(0,nt-1)    
ax02.set_ylim(0,nx)
ax02.set_title('C: background intensity from Ca II K wings', fontdict=font)

ax10.set_ylabel(r'length [arcsec]', fontdict = font)
ax10.set_xticks(ttick_pos)
ax10.set_xticklabels([])
ax10.xaxis.set_minor_locator(AutoMinorLocator(4))
ax10.set_xlabel('', fontdict = font)
ax10.set_yticks(xtick_pos)
ax10.set_yticklabels(xtick_lab, fontdict = font)
ax10.yaxis.set_minor_locator(AutoMinorLocator(5))
ax10.set_xlim(0,nt-1)
ax10.set_ylim(0,nx)
ax10.set_title('D: Cleaned = B + alpha * (B-C)', fontdict=font)

ax11.set_ylabel('')
ax11.set_xticks(ttick_pos)
ax11.set_xticklabels([])
ax11.xaxis.set_minor_locator(AutoMinorLocator(4))
ax11.set_xlabel('')
ax11.set_yticks(xtick_pos)
ax11.set_yticklabels([])
ax11.yaxis.set_minor_locator(AutoMinorLocator(5))
ax11.set_xlim(0,nt-1)    
ax11.set_ylim(0,nx)
ax11.set_title('E: Gaussian filtered with kernel [1,3]', fontdict=font)

ax12.set_ylabel('')
ax12.set_xticks(ttick_pos)
ax12.set_xticklabels([])
ax12.xaxis.set_minor_locator(AutoMinorLocator(4))
ax12.set_xlabel('')
ax12.set_yticks(xtick_pos)
ax12.set_yticklabels([])
ax12.yaxis.set_minor_locator(AutoMinorLocator(5))
ax12.set_xlim(0,nt-1)    
ax12.set_ylim(0,nx)
ax12.set_title('F: Gaussian filtered with kernel [3,1]', fontdict=font)

ax20.set_ylabel(r'length [arcsec]', fontdict = font)
ax20.set_xticks(ttick_pos)
ax20.set_xticklabels(ttick_lab, fontdict = font)
ax20.xaxis.set_minor_locator(AutoMinorLocator(4))
ax20.set_xlabel('t [min]', fontdict = font)
ax20.set_yticks(xtick_pos)
ax20.set_yticklabels(xtick_lab, fontdict = font)
ax20.yaxis.set_minor_locator(AutoMinorLocator(5))
ax20.set_xlim(0,nt-1)
ax20.set_ylim(0,nx)
ax20.set_title('G: Sharpened = E + alpha * (E-F)', fontdict=font)

ax21.set_ylabel('')
ax21.set_xticks(ttick_pos)
ax21.set_xticklabels(ttick_lab, fontdict = font)
ax21.xaxis.set_minor_locator(AutoMinorLocator(4))
ax21.set_xlabel('t [min]', fontdict = font)
ax21.set_yticks(xtick_pos)
ax21.set_yticklabels([])
ax21.yaxis.set_minor_locator(AutoMinorLocator(5))
ax21.set_xlim(0,nt-1)
ax21.set_ylim(0,nx)
ax21.set_title('H: Gamma-corrected & sharpened', fontdict=font)

ax22.set_ylabel('')
ax22.set_xticks(ttick_pos)
ax22.set_xticklabels(ttick_lab, fontdict = font)
ax22.xaxis.set_minor_locator(AutoMinorLocator(4))
ax22.set_xlabel('t [min]', fontdict = font)
ax22.set_yticks(xtick_pos)
ax22.set_yticklabels([])
ax22.yaxis.set_minor_locator(AutoMinorLocator(5))
ax22.set_xlim(0,nt-1)
ax22.set_ylim(0,nx)
ax22.set_title('G: G-corrected and sharpened with new oscillations paths', fontdict=font)

gs.update(left=0.035,
          right=0.99,
          wspace=0.025,
          bottom=0.05,
          top=0.98,
          hspace = 0.1,
)

for uu in range(len(file_search(outdir_old, "cut"+cut_file[-11:-5]+"*txt"))):

    osc_fname_old = outdir_old+file_search(outdir_old, "cut"+cut_file[-11:-5]+"*txt")[uu]
    coord_old = np.loadtxt(osc_fname_old)
    amin, amax = int(np.min(coord_old[:,0])), int(np.max(coord_old[:,0]))+1
    tmin, tmax = int(np.min(coord_old[:,1])), int(np.max(coord_old[:,1]))+1
    trange = np.linspace(tmin, tmax, tmax-tmin+1)
    
    ax00.plot(coord_old[:,1], coord_old[:,0],color = 'wheat', alpha = 0.5, label = 'Manual fit', linestyle = '--', linewidth = 1)
    

for oo in range(len(file_search(outdir, "cut"+cut_file[-11:-5]+"*txt"))):

    osc_fname = outdir+file_search(outdir, "cut"+cut_file[-11:-5]+"*txt")[oo]
    coord = np.loadtxt(osc_fname)
    amin, amax = int(np.min(coord[:,0])), int(np.max(coord[:,0]))+1
    tmin, tmax = int(np.min(coord[:,1])), int(np.max(coord[:,1]))+1
    trange = np.linspace(tmin, tmax, tmax-tmin+1)
    
    ax22.plot(coord[:,1], coord[:,0],color = 'wheat', alpha = 0.5, label = 'Manual fit', linestyle = '--', linewidth = 1)
    
plt.show()

filename = outdir + 'oscillation_curvefit_' +cut_file[-11:-5]+'_enhanced_cleaned.pdf'
fig.savefig(filename, quality = 100)
print 'file saved to: ' + filename

# 3D view of the cuts
if(0):
    #plt.close("all")
    thd = plt.figure()
    ax = thd.add_subplot(111, projection='3d')
    y = np.arange(0, cube.shape[0], 1)
    x = np.arange(0, cube.shape[1], 1)
    X, Y = np.meshgrid(x, y)
    ax.plot_surface(X, Y,
                    cube.squeeze(),
                    cmap = "gray",
                    alpha = 0.75,
                    vmin=1.2,
                    vmax=2.,
                    )
    ax.invert_yaxis()
    plt.tight_layout()
    #plt.show()

