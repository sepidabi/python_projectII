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
cube = np.transpose(np.mean(cut.loop_slab[w_pos-3:w_pos+3,:,:],axis = 0)).squeeze()*1e9

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
fig = plt.figure(figsize = (4,9.5))#4.6))
gs = gridspec.GridSpec(4,1, figure = fig)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])
ax3 = plt.subplot(gs[2,0])
ax4 = plt.subplot(gs[3,0])

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
ax1.set_title('(a)', fontdict=font)

ax2.set_ylabel(r'length [arcsec]', fontdict = font)
ax2.set_xticks(ttick_pos)
ax2.set_xticklabels([])
ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
ax2.set_xlabel('', fontdict = font)
ax2.set_yticks(xtick_pos)
ax2.set_yticklabels(xtick_lab, fontdict = font)
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.set_xlim(0,nt-1)
ax2.set_ylim(0,nx)
ax2.set_title('(b)', fontdict=font)

ax3.set_ylabel(r'length [arcsec]', fontdict = font)
ax3.set_xticks(ttick_pos)
ax3.set_xticklabels([])
ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
ax3.set_xlabel('', fontdict = font)
ax3.set_yticks(xtick_pos)
ax3.set_yticklabels(xtick_lab, fontdict = font)
ax3.yaxis.set_minor_locator(AutoMinorLocator(5))
ax3.set_xlim(0,nt-1)
ax3.set_ylim(0,nx)
ax3.set_title('(c)', fontdict=font)

ax4.set_ylabel(r'length [arcsec]', fontdict = font)
ax4.set_xticks(ttick_pos)
ax4.set_xticklabels(ttick_lab, fontdict = font)
ax4.xaxis.set_minor_locator(AutoMinorLocator(4))
ax4.set_xlabel('t [min]', fontdict = font)
ax4.set_yticks(xtick_pos)
ax4.set_yticklabels(xtick_lab, fontdict = font)
ax4.yaxis.set_minor_locator(AutoMinorLocator(5))
ax4.set_xlim(0,nt-1)
ax4.set_ylim(0,nx)
ax4.set_title('(d)', fontdict=font)

gs.update(left=0.1,
          right=0.99,
          #wspace=0.05,
          bottom=0.05,
          top=0.98,
          hspace = 0.15,
)

# IMAGE ENHANCEMENT
#-----------------------------
#cube_gamma = exposure.adjust_log(cube,1)#gamma(cube,0.1)-9. # gamma corrected
#cube = np.where(cube>0., cube,1e-4)
# median filtering
# background intensity
uniform = 0
if (uniform==1):
    wg = 4
    cube_bg = wg*np.std(cube)*np.ones(cube.shape)
else:
    wg = 0.15
    cube_bg = np.transpose(wg*np.mean(cut.loop_slab[:w_pos-6,:,:]+cut.loop_slab[w_pos+8:,:,:],axis=0).squeeze()*1e9)

cube_trunc = cube - cube_bg
cube_trunc_med = exposure.rescale_intensity(sgnl.medfilt2d(cube_trunc,kernel_size = [1,3]), out_range=(0, 1.)) # removing the background intensity
cube_sharp = exposure.rescale_intensity(sharpen(cube_trunc, sigma =[1,3], unsigma = [3,1]), out_range=(0, 1.))
cube_sharp_gamma = exposure.rescale_intensity(exposure.adjust_gamma(cube_sharp,0.1), out_range=(0, 1.))

# Panels
panel1 = cube_trunc_med
panel2 = cube_sharp
panel3 = cube_sharp_gamma
panel4 = panel3

# PLOTTING results
im_cak = ax1.imshow(panel1,
                    cmap = 'gray', origin = 'lower', #aspect=ratio
                    interpolation = 'bilinear',
                    vmin=0.05,
                    vmax = 0.9,
                    aspect = 0.6,
)

im_cak2 = ax2.imshow(panel2,#cube_gamma,
                     cmap = 'gray', origin = 'lower', #aspect=ratio
                     #interpolation = 'hermite',
                     #vmin=0.,
                     #vmax = 1.,
                     aspect = 0.6,
)

im_cak3 = ax3.imshow(panel3,#cube_gamma,
                     cmap = 'gray', origin = 'lower', #aspect=ratio
                     #interpolation = 'hermite',
                     vmin=0.89,
                     vmax = 0.98,
                     aspect = 0.6,
)

im_cak4 = ax4.imshow(panel4,#cube_gamma,
                     cmap = 'gray', origin = 'lower', #aspect=ratio
                     #interpolation = 'hermite',
                     vmin=0.89,
                     vmax = 0.98,
                     aspect = 0.6,
)

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

for oo in range(len(file_search(outdir, "cut"+cut_file[-11:-5]+"*txt"))):

    osc_fname = outdir+file_search(outdir, "cut"+cut_file[-11:-5]+"*txt")[oo]
    coord = np.loadtxt(osc_fname)
    tt = coord[:,1]
    xx = coord[:,0]
    yy = cube[np.int_(np.rint(xx)),np.int_(np.rint(tt))].squeeze()*1e9
    amin, amax = int(np.min(xx)), int(np.max(xx))+1
    tmin, tmax = int(np.min(tt)), int(np.max(tt))+1
    trange = np.linspace(tmin, tmax, tmax-tmin+1)
    edge = 5 #pixels
    tt0, ttn  = np.int(tt[0]), np.int(tt[-1])+1#int(np.max([np.min(tt),0])), int(np.min([np.max(tt)+edge,nt]))
    xx0, xxn = int(np.max([np.min(xx)-edge,0])), int(np.min([np.max(xx)+edge,nx]))
    osc_panel1 = panel1[xx0:xxn, tt0:ttn+1]
    osc_panel2 = panel4[xx0:xxn, tt0:ttn+1]

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
        curve_int = panel3[int(y_imax[i])-eps:int(y_imax[i])+eps, i]
        xdata = np.arange(0,len(curve_int),1.)

        # calculate Gaussian fit
        #popt, pcov = curve_fit(gauss, xdata, curve_int, maxfev = 1000000)
        #curve_gauss = gauss(xdata, *popt)
        #y_gmax[i] = np.argmax(curve_gauss)+int(y_imax[i])-eps
        #y_gmax_val = np.max(curve_gauss)
        
        
        '''plt.figure()
        plt.plot(xdata, curve_int, 'r-')
        plt.plot(xdata, curve_gauss, 'b.')
        #plt.plot(y_imax_val,curve_int[eps],'r.')
        plt.plot(y_gmax_val,curve_gauss[int(y_gmax[i])],'b.')
        plt.show()

        stop()'''
    smth_y = smooth(y_imax, 6)
    #smth_y_gauss = smooth(y_gmax, 6)
    osc_fname_smth = outdir+"smth"+osc_fname[31:]
    np.savetxt(osc_fname_smth, smth_y, fmt='%3.8f', encoding=None)

    # Calculation of the velocity of the oscillation in the plane of sky
    xp = trange*cad
    fxp = smth_y*res*725
    der = derivative(xp, fxp, dx = 1e-5)
    x, fx, f1x = der[0]-np.min(der[0]), der[1], der[2]

    
    ax4.plot(coord[:,1], coord[:,0],color = 'wheat', alpha = 0.5, label = 'Manual fit', linestyle = '--', linewidth = 1)
    #if (oo ==7 or oo==5 or oo==12 or oo==6):
        #ax4.plot(trange, smth_y, color = 'red', alpha = 0.5, label = r'smthd I$_{\rm{max}}$', linewidth = 2)
        #if (oo==6):
            #ax4.text(trange[0]+1, smth_y[0]+1.5, "1", color = 'red')
        #if (oo==12):
            #ax4.text(trange[0]+1, smth_y[0]+1.5, "2", color = 'red')
        #if (oo==5):
            #ax4.text(trange[0]+1, smth_y[0]+1.5, "3", color = 'red')
        #if (oo==7):
            #ax4.text(trange[0]+1, smth_y[0]+1.5, "4", color = 'red')

    #elif (oo==6):
        #print("nothing!")
    #else:
        #ax2.plot(coord[:,1], coord[:,0],color = 'wheat', alpha = 0.5, label = 'Manual fit', linestyle = '--', linewidth = 1)
    ax4.plot(trange, smth_y, color = 'orange', alpha = 0.5, label = r'smthd I$_{\rm{max}}$', linewidth = 2)
        #ax4.plot(trange, smth_y_gauss, color = 'cyan', alpha = 0.5, label = r'gauss fit', linewidth = 2)

        
    # 3D plot
    #ax.plot(tt, xx, yy, '-', alpha=0.75, color = 'red')

    # plot single oscillation
    if(0):
        # PLOTTING each oscillation isolated:
        osc_fig = plt.figure(figsize = (6.73,4.))
        gs_osc = gridspec.GridSpec(2,1, figure = osc_fig)
        ax_osc0 = plt.subplot(gs_osc[0,0])
        ax_osc = plt.subplot(gs_osc[1,0])
        p_fact = 1.75
    
        # single oscillation plot
        osc_im0= ax_osc0.imshow(osc_panel1, origin = 'lower',
                                aspect = 0.6,#(0.3)*(ttn-tt0)/(xxn-xx0),
                                cmap = 'gray',
                                #vmin = vmin,
                                #vmax = vmax,
                                )
    
        # normalizing each column in oscillations to the maximum path
        #col=0
        #for col in range (len(y_interp)-1):
        #   osc_cube_norm[:,col] = osc_cube[:,col]/osc_cube[np.int(y_interp[col])-xx0, col]
        #print(col)
        #osc_cube_norm = np.where(osc_cube_norm<=1,osc_cube_norm,1.)
        
        vmin=0.9
        vmax = 0.995
        
        osc_im = ax_osc.imshow(osc_panel2, origin = 'lower',
                               aspect = 0.6,#(0.3)*(ttn-tt0)/(xxn-xx0),
                               cmap = 'gray',
                               vmin = vmin,
                               vmax = vmax,
                               )
        
        #ax_osc.plot(tt-tt0, xx-xx0, color = 'orange', alpha = 0.5, label = 'Manual fit', linestyle = '--', linewidth = 1.5)
        #ax_osc.plot(trange-tt0, smth_y-xx0, color = 'orange', alpha = 0.5, label = r'smthd I$_{\rm{max}}$', linewidth = 2)
        #ax_osc.plot(trange-tt0, smth_y_gauss-xx0, color = 'cyan', alpha = 0.5, label = r'Guassian fit', linewidth = 2)
        #ax_osc.plot(trange-tt0, smth_y-xx0, color = 'orange', alpha = 0.5, label = r'I$_{max}$', linewidth = 2)
        
        
        
        gs_osc.update(left=0.1,
                      right=0.99,
                      wspace=0.0,
                      bottom=0.1,
                      top=0.99,
                      hspace = 0.05,
                      )
        
        
        # plotting elements
        t_fact = 30
        ntt, nxx = osc_panel1.shape[1], osc_panel1.shape[0]
        ttick_pos = np.arange(0,(np.round((ntt)/t_fact)+1)*t_fact,t_fact)
        ttick_lab = ttick_pos*cad/60
        rr = 40/1.5
        xtick_pos = np.arange(0,(np.round((nxx)/rr)+1)*rr,rr)
        xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))
        
        # plot settings
        ax_osc0.set_ylabel(r'length [arcsec]', fontdict = font)
        ax_osc0.set_xticks(ttick_pos)
        ax_osc0.set_xticklabels([])
        ax_osc0.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax_osc0.set_xlabel('')
        ax_osc0.set_yticks(xtick_pos)
        #ax_osc0.set_yticklabels(xtick_lab, fontdict = font)
        ax_osc0.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax_osc0.set_xlim(0,ntt-1)
        ax_osc0.set_ylim(0,nxx-1)
        
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
        
        plt.tight_layout()
        #plt.show()
        
        #stop()
        plt.close(osc_fig)
        
        
plt.show()

filename = outdir + 'oscillation_curvefit_' +cut_file[-11:-5]+'_non_intt.pdf'
fig.savefig(filename, quality = 100)
print 'file saved to: ' + filename
