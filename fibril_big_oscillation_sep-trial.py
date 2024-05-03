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
#def test_func(x, a, b, c, d):
 #   return a * np.sin(b * x - c)+d

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
nt = 213
nx = 201

cut_file = (file_search(fibdir,'crispex*3950*.csav'))[cc]
cut = restore(fibdir+cut_file)
cube = np.transpose(cut.loop_slab[w_pos,:,:]).squeeze()*1e9

xx = cut.loop_slab.shape[2]

# plotting elements
tt = 30
ttick_pos = np.arange(0,(np.round((nt)/tt)+1)*tt,tt)
ttick_lab = ttick_pos*cad/60

rr = 40/1.5
xtick_pos = np.arange(0,(np.round((xx)/rr)+1)*rr,rr)
xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))

# intensity oscillation PLOTs
plt.close('all')
fig = plt.figure(figsize = (4,4.6))
gs = gridspec.GridSpec(2,1, figure = fig)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])

ax1.set_ylabel(r'length [arcsec]', fontdict = font)
ax1.set_xticks(ttick_pos)
ax1.set_xticklabels([])
ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
ax1.set_xlabel('')
ax1.set_yticks(xtick_pos)
ax1.set_yticklabels(xtick_lab, fontdict = font)
ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
ax1.set_xlim(0,nt-1)    
ax1.set_ylim(0,nx)

ax2.set_ylabel(r'length [arcsec]', fontdict = font)
ax2.set_xticks(ttick_pos)
ax2.set_xticklabels(ttick_lab, fontdict = font)
ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
ax2.set_xlabel('t [min]', fontdict = font)
ax2.set_yticks(xtick_pos)
ax2.set_yticklabels(xtick_lab, fontdict = font)
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.set_xlim(0,nt-1)
ax2.set_ylim(0,nx)

gs.update(left=0.1,
          right=0.98,
          wspace=0.0,
          bottom=0.09,
          top=0.995,
          hspace = 0.05,
)

# Intensity
#------------
cube_gamma = gamma(cube,0.1)-9. # gamma corrected
cube_med = sci.signal.medfilt2d(cube) # median filtered
cube_gm =  sci.signal.medfilt2d(cube_med) # gamma corrected and median filtered

im_cak = ax1.imshow(cube,
                    cmap = 'gray', origin = 'lower', #aspect=ratio
                    #interpolation = 'bilinear',
                    #vmin=1.2,
                    #vmax = 1.9,
                    aspect = 0.6,
)

im_cak2 = ax2.imshow(cube_gm,#cube_gamma,
                     cmap = 'gray', origin = 'lower', #aspect=ratio
                     #interpolation = 'hermite',
                     vmin=1.2,
                     vmax = 1.9,
                     aspect = 0.6,
)

# 3D view of the cuts
if(1):
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


for oo in range(10, len(file_search(outdir, "cut"+cut_file[-11:-5]+"*txt"))):

    osc_fname = outdir+file_search(outdir, "cut"+cut_file[-11:-5]+"*txt")[oo]
    coord = np.loadtxt(osc_fname)
    tt = coord[:,1]
    xx = coord[:,0]
    yy = cube[np.int_(np.rint(xx)),np.int_(np.rint(tt))].squeeze()*1e9
    amin, amax = int(np.min(xx)), int(np.max(xx))+1
    tmin, tmax = int(np.min(tt)), int(np.max(tt))+1
    trange = np.linspace(tmin, tmax, tmax-tmin+1)
    edge = 5 #pixels
    tt0, ttn  = int(np.max([np.min(tt)-edge,0])), int(np.min([np.max(tt)+edge,nt]))
    xx0, xxn = int(np.max([np.min(xx)-edge,0])), int(np.min([np.max(xx)+edge,nx]))
    display_cube = cube_gm
    osc_cube = display_cube[xx0:xxn, tt0:ttn]

    # interpolating the points in between the clicked fit
    y_interp = np.interp(trange, tt, xx)
    osc_fname_interp = outdir+"interp"+osc_fname[31:]
    np.savetxt(osc_fname_interp, y_interp, fmt='%3.8f', encoding=None)
    
    # track max intensity based on the interp click_fit
    y_imax = np.zeros(tmax-tmin+1)
    y_gmax = np.zeros(tmax-tmin+1)
    dist = 1. # width of the interval in which the max intensity is taken

    for i in range(tmax-tmin+1):

        #Gaussian function 
        def gauss(x, a, x0, sigma): 
            return a*np.exp(-(x-x0)**2/(2*sigma**2)) 


        y_imax[i] = np.argmax(cube[int(y_interp[i]-dist):int(y_interp[i]+dist), i])+int(y_interp[i])-dist
        y_imax_val = np.max(cube[int(y_interp[i]-dist):int(y_interp[i]+dist), i])+int(y_interp[i])
        eps = int(10) # within the range to look for maximum
        curve_int = cube[int(y_imax[i])-eps:int(y_imax[i])+eps, i]
        xdata = np.arange(0,len(curve_int),1.)
        popt, pcov = curve_fit(gauss, xdata, curve_int, maxfev = 1000000)
        curve_gauss = gauss(xdata, *popt)
        y_gmax[i] = np.argmax(curve_gauss)+int(y_imax[i])-eps
        y_gmax_val = np.max(curve_gauss)
        '''plt.figure()
        plt.plot(xdata, curve_int, 'r-')
        plt.plot(xdata, curve_gauss, 'b.')
        #plt.plot(y_imax_val,curve_int[eps],'r.')
        plt.plot(y_gmax_val,curve_gauss[int(y_gmax[i])],'b.')
        plt.show()

        stop()'''
    smth_y = smooth(y_imax, 6)
    smth_y_gauss = smooth(y_gmax, 6)
    osc_fname_smth = outdir+"smth"+osc_fname[31:]
    np.savetxt(osc_fname_smth, smth_y, fmt='%3.8f', encoding=None)

    # Calculation of the velocity of the oscillation in the plane of sky
    xp = trange*cad
    fxp = smth_y*res*725
    der = derivative(xp, fxp, dx = 1e-5)
    x, fx, f1x = der[0]-np.min(der[0]), der[1], der[2]

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
        ax2.plot(trange, smth_y_gauss, color = 'cyan', alpha = 0.5, label = r'gauss fit', linewidth = 2)

    # 3D plot
    #ax.plot(tt, xx, yy, '-', alpha=0.75, color = 'red')

    # PLOTTING each oscillation isolated:
    osc_fig = plt.figure(figsize = (4,1.5))
    gs_osc = gridspec.GridSpec(1,1, figure = osc_fig)
    ax_osc = plt.subplot(gs_osc[0,0])
    p_fact = 1.75
    osc_im = ax_osc.imshow(osc_cube, origin = 'lower',
                       aspect = (0.3)*(ttn-tt0)/(xxn-xx0),
                       cmap = 'gray',
                       vmin = np.mean(osc_cube)-1.75*np.std(osc_cube),
                       vmax = np.mean(osc_cube)+2*np.std(osc_cube),
    )
    ax_osc.plot(tt-tt0, xx-xx0, color = 'wheat', alpha = 0.5, label = 'Manual fit', linestyle = '--', linewidth = 1)
    #ax_osc.plot(trange-tt0, smth_y-xx0, color = 'orange', alpha = 0.5, label = r'smthd I$_{\rm{max}}$', linewidth = 2)
    ax_osc.plot(trange-tt0, y_gmax-xx0, color = 'red', alpha = 0.5, label = r'Guassian fit', linewidth = 2)
    ax_osc.plot(trange-tt0, y_imax-xx0, color = 'orange', alpha = 0.5, label = r'I$_{max}$', linewidth = 2)

    gs_osc.update(left=0.1,
          right=0.98,
          wspace=0.0,
          bottom=0.25,
          top=0.98,
          hspace = 0.0,
            )

    
    # plotting elements
    t_fact = 30
    ntt, nxx = osc_cube.shape[1], osc_cube.shape[0]
    ttick_pos = np.arange(0,(np.round((ntt)/t_fact)+1)*t_fact,t_fact)
    ttick_lab = ttick_pos*cad/60
    rr = 40/1.5
    xtick_pos = np.arange(0,(np.round((nxx)/rr)+1)*rr,rr)
    xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))
    # plot settings
    ax_osc.set_ylabel(r'length [arcsec]', fontdict = font)
    ax_osc.set_xticks(ttick_pos)
    #ax_osc.set_xticklabels(ttick_lab, fontdict = font)
    ax_osc.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax_osc.set_xlabel('t [min]', fontdict = font)
    ax_osc.set_yticks(xtick_pos)
    #ax_osc.set_yticklabels(xtick_lab, fontdict = font)
    ax_osc.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax_osc.set_xlim(0,ntt-1)
    ax_osc.set_ylim(0,nxx-1)

    #plt.tight_layout()
    plt.show()

    stop()
    plt.close(osc_fig)

plt.show()

#filename = outdir + 'oscillation_curvefit_' +cut_file[-11:-5]+'_non_intt.pdf'
#fig.savefig(filename, quality = 100)
#printt 'file saved to: ' + filename
