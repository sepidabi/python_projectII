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

#SCRIPT MODE
# 0 => oscillation sequence
# 1 => image manipulation of the oscillation
# 2 => video of ROI
# 3 => curve fitting to the the oscillation
mode = 3

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
ind = 0
fibdir = datadir+'fr29/'
fib_file = (file_search(fibdir,'crispex*3950*.csav'))[0+ind]
fib = restore(fibdir+fib_file)
cut1_file = (file_search(fibdir,'crispex*3950*.csav'))[1+ind]
cut2_file = (file_search(fibdir,'crispex*3950*.csav'))[2+ind]
cut3_file = (file_search(fibdir,'crispex*3950*.csav'))[3+ind]
#cut4_file = (file_search(fibdir,'crispex*3950*.csav'))[4]

cut1 = restore(fibdir+cut1_file)
cut2 = restore(fibdir+cut2_file)
cut3 = restore(fibdir+cut3_file)

print(np.sqrt((cut1.x_coords[0]-cut1.x_coords[1])**2 + (cut1.y_coords[0]-cut1.y_coords[1])**2),
      cut1.loop_size)
print(np.sqrt((cut2.x_coords[0]-cut2.x_coords[1])**2 + (cut2.y_coords[0]-cut2.y_coords[1])**2),
      cut2.loop_size)
print(np.sqrt((cut3.x_coords[0]-cut3.x_coords[1])**2 + (cut3.y_coords[0]-cut3.y_coords[1])**2),
      cut3.loop_size)

xmin = int(np.min([np.min([cut1.x_coords, cut2.x_coords,
                           #cut3.x_coords,
                           cut3.x_coords]),np.min(fib.x_coords)]) - 5)
ymin = int(np.min([np.min([cut1.y_coords, cut2.y_coords,
                           #cut3.y_coords,
                           cut3.y_coords]),np.min(fib.y_coords)]) - 5)
xmax = int(np.max([np.max([cut1.x_coords, cut2.x_coords,
                           #cut3.x_coords,
                           cut3.x_coords]),np.max(fib.x_coords)]) + 5)
ymax = int(np.max([np.max([cut1.y_coords, cut2.y_coords,
                           #cut3.y_coords,
                           cut3.y_coords]),np.max(fib.y_coords)]) + 5)

print(xmin, xmax)
print(ymin, ymax)

#######################
## plot oscilations
#######################

if (mode==0):
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
    ## PLOTS
    #######################
    vmin = None#2.e-10 #1.9180524e-10
    vmax = None#1.e-8 #1.09262785e-08
    
    ax11.imshow(cube_cak[fr_select[0],0,w_pos,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower',
                vmin = vmin, vmax=vmax)
    ax11.text(1, 1, str(fr_select[0]*8/60)+' m', color = 'red')
    #ax11.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'red')
    #ax11.plot(cut1.x_coords-xmin, cut1.y_coords-ymin, color = 'red')
    #ax11.plot(cut2.x_coords-xmin, cut2.y_coords-ymin, color = 'red')
    #plt.plot(cut3.x_coords-xmin, cut3.y_coords-ymin, color = 'red')
    #ax11.plot(cut3.x_coords-xmin, cut3.y_coords-ymin, color = 'red')
    
    ax12.imshow(cube_cak[fr_select[1],0,w_pos,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower',
                vmin = vmin, vmax=vmax)
    ax12.text(1, 1, str(fr_select[1]*8/60)+' m', color = 'red')
    #ax12.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'red')
    #ax12.plot(cut1.x_coords-xmin, cut1.y_coords-ymin, color = 'red')
    #ax12.plot(cut2.x_coords-xmin, cut2.y_coords-ymin, color = 'red')
    #plt.plot(cut3.x_coords-xmin, cut3.y_coords-ymin, color = 'red')
    #ax12.plot(cut3.x_coords-xmin, cut3.y_coords-ymin, color = 'red')
    
    ax13.imshow(cube_cak[fr_select[2],0,w_pos,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower',
                vmin = vmin, vmax=vmax)
    ax13.text(1, 1, str(fr_select[2]*8/60)+' m', color = 'red')
    #ax13.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'red')
    #ax13.plot(cut1.x_coords-xmin, cut1.y_coords-ymin, color = 'red')
    #ax13.plot(cut2.x_coords-xmin, cut2.y_coords-ymin, color = 'red')
    
    ax14.imshow(cube_cak[fr_select[3],0,w_pos,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower',
                vmin = vmin, vmax=vmax)
    ax14.text(1, 1, str(fr_select[3]*8/60)+' m', color = 'red')
    ax14.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'red')
    ax14.plot(cut1.x_coords-xmin, cut1.y_coords-ymin, color = 'red')
    ax14.plot(cut2.x_coords-xmin, cut2.y_coords-ymin, color = 'red')
    ax14.plot(cut3.x_coords-xmin, cut3.y_coords-ymin, color = 'red')
    ax14.text(cut1.x_coords[0]-xmin-4, cut1.y_coords[0]-ymin-4, '(3)', fontdict=font, color = 'red')
    ax14.text(cut2.x_coords[0]-xmin-4, cut2.y_coords[0]-ymin-4, '(2)', fontdict=font, color = 'red')
    ax14.text(cut3.x_coords[0]-xmin-4, cut3.y_coords[0]-ymin-4, '(1)', fontdict=font, color = 'red')
    
    
    ratio =0.8 # DO NOT CHANGE THIS RATIO otherwise the output image goes bananas!
    ax2.imshow(cut3.loop_slab[0,w_pos,:,:], cmap = 'gray', origin = 'lower', aspect=ratio,
               vmin = vmin, vmax=vmax)
    ax2.text(2,205,'(1)', color = 'red')
    
    ax3.imshow(cut2.loop_slab[0,w_pos,:,:], cmap = 'gray', origin = 'lower', aspect=ratio,
               vmin = vmin, vmax=vmax)
    ax3.text(2,205,'(2)', color = 'red')
    
    ax4.imshow(cut1.loop_slab[0,w_pos,:,:], cmap = 'gray', origin = 'lower', aspect=ratio,
               vmin = vmin, vmax=vmax)
    ax4.text(2,205,'(3)', color = 'red')
    
    plt.show()
    
    # wpos
    if (w_pos == 13): wpos='K3'
    elif (w_pos==10): wpos = 'K2V'
    else: wpos = 'K2R'
    
    filename = outdir+'oscilation'+wpos+'.pdf'
    f.savefig(filename, quality = 100)
    f.tight_layout()
    print 'file saved to: ' + filename



#####################
#  image manipulation
#####################
if(mode==1):

    cut = cut2
    cut_file = cut2_file
    plt.close('all')
    
    amin = 2.7e-9
    amax = 3.e-9
    cube = cut.loop_slab[14,:,:]
    cube_man = spnd.gaussian_filter(cube,1)
    cube_clip = np.clip(cube, a_min = amin,
                        a_max = None
    )
    cube_binary = np.zeros(cube.shape)
    cube_binary[np.where(cube>amin)] = 1
    cube_man_clip = np.clip(cube_man, a_min = amin,
                            a_max = None
    )
    
    cube_man_binary = np.zeros(cube.shape)
    cube_man_binary[np.where(cube_man>amin)] = 1

    #######################
    ## PLOTS
    #######################
    vmin = None#2.e-10 #1.9180524e-10
    vmax = None#1.e-8 #1.09262785e-08

    f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(1, 7, figsize = (15,10))
    ax1.imshow(cube, cmap = 'gray', origin = 'lower', #aspect=ratio
    )
    ax1.plot(np.argmax(cube, axis = 1), np.arange(0,213,1.), color = 'red')
    
    
    ax2.imshow(cube_man, cmap = 'gray', origin = 'lower', #aspect=ratio,
               vmin = vmin,
               vmax = vmax
    )
    ax2.plot(np.argmax(cube_man, axis = 1), np.arange(0,213,1.), color = 'red')
    
    
    ax3.imshow(cube_man, cmap = 'gray', origin = 'lower', #aspect=ratio,
               vmin = vmin,
               vmax = vmax,
               interpolation = 'bilinear'
    )
    #ax3.plot(np.argmax(cube_man, axis = 1), np.arange(0,213,1.), color = 'red')
    
    
    ax4.imshow(cube_clip, cmap = 'gray', origin = 'lower',
    )
    #ax4.plot(np.argmax(cube_clip, axis = 1), np.arange(0,213,1.), color = 'red')
    
    ax5.imshow(cube_binary, cmap = 'gray', origin = 'lower', #aspect=8
    )
    
    ax6.imshow(cube_man_clip, cmap = 'gray', origin = 'lower', #aspect=8
    )
    #ax6.plot(np.argmax(cube_man_clip, axis = 1), np.arange(0,213,1.), color = 'red')
    
    
    ax7.imshow(cube_man_binary, cmap = 'gray', origin = 'lower', #aspect=8
    )
    
    ax1.set_title('obs.')
    ax2.set_title('smoothed')
    ax3.set_title('interpolated')
    ax4.set_title('clipped')
    ax5.set_title('binary')
    ax6.set_title('sm. + cl.')
    ax7.set_title('sm. + bi.')
    
    plt.show()
    filename = outdir + 'oscillation_image_manipulation.pdf'
    #f.savefig(filename, quality = 100)
    print 'file saved to: ' + filename

######################
# make video of the roi
######################

vmin = None#2.e-10 #1.9180524e-10
vmax = None#1.e-8 #1.09262785e-08

if(mode==2):
    for ii in range(nt):
        plt.close('all')
        f, ax13 = plt.subplots()
        rr = 40/1.5
        ytick_pos = np.arange(0,(np.round((ymax-ymin)/rr)+1)*rr,rr)
        ytick_lab = np.intc(np.round(ytick_pos*res, decimals = 0))
        xtick_pos = np.arange(0,(np.round((xmax-xmin)/rr)+1)*rr,rr)
        xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))

        ax13.imshow(cube_cak[ii,0,w_pos,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower',
                    vmin = vmin, vmax=vmax)
        ax13.text(1, 1, str(ii*8)+' s', color = 'red')
        ax13.set_xlabel(r'x [arcsec]', fontdict = font)
        ax13.set_ylabel(r'y [arcsec]', fontdict = font)
        ax13.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax13.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax13.tick_params(which='minor', length=2)
        ax13.set_xticks(xtick_pos)
        ax13.set_xticklabels(xtick_lab, fontdict = font)
        ax13.set_yticks(ytick_pos)
        ax13.set_yticklabels(ytick_lab, fontdict = font)
        
        plt.show()
        
        filename = outdir+'video2/'+'roi_'+str(ii).zfill(4)+'.png'
        #f.savefig(filename, quality=100)
        print 'file saved to: ' + filename
        


########################
# CURVE FITTING
########################
if(mode==3):
    cut = cut2
    cut_file = cut2_file
    cube_raw = np.mean(cut.loop_slab[w_pos-3:w_pos+3,:,:],axis = 0).squeeze()*1e9
    cut_fname = outdir+cut_file[-11:-5]+'.txt'

    # background intensity from the wings
    wg = 0.15
    cube_bg = wg*np.mean(cut.loop_slab[:w_pos-6,:,:]+cut.loop_slab[w_pos+8:,:,:],axis=0).squeeze()*1e9
    cube_trunc = cube_raw - cube_bg # removing the background intensity
    # contrast correction
    cube_med = exposure.rescale_intensity(sgnl.medfilt2d(cube_trunc,kernel_size = [3,1]), out_range=(0, 1.))
    cube_sharp = exposure.rescale_intensity(sharpen(cube_trunc, sigma =[3,1], unsigma = [1,3]), out_range=(0, 1.))
    cube = exposure.rescale_intensity(exposure.adjust_gamma(cube_sharp,0.1), out_range=(0, 1.))
    vmin, vmax = 0.9, 0.97

    
    # to test the curve fitting
    from scipy.optimize import curve_fit
    ti = np.arange(0, nt, 1.)
    yi = np.argmax(cube, axis = 1)
    
    def test_func(x, a, b, c, d):
        return a * np.sin(b * x - c)+d
    
    params, params_covariance = curve_fit(test_func, ti, yi,
                                          p0=[17., 0.075, -10., 17.])
                                          #sigma = None)#np.repeat(2, ti.size),
    print(params)
    
    #boxcar smoothing
    def smooth(y, w):
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
            

    # plotting

    tt = 30
    ttick_pos = np.arange(0,(np.round((213)/tt)+1)*tt,tt)
    ttick_lab = ttick_pos*cad/60
    
    rr = 40/1.5
    xtick_pos = np.arange(0,(np.round((xmax-xmin)/rr)+1)*rr,rr)
    xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))
    xtick_lab2 = np.intc(np.round(xtick_pos*res*(np.sqrt((cut.x_coords[0]-cut.x_coords[1])**2 +
                                                         (cut.y_coords[0]-cut.y_coords[1])**2)/cut.loop_size),
                                  decimals = 0))
    
    plt.close('all')
    f, ax1 = plt.subplots(figsize = (5,10))
    ax1.set_xlabel(r'length [arcsec]', fontdict = font)
    ax1.set_yticks(ttick_pos)
    ax1.set_yticklabels(ttick_lab, fontdict = font)
    ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.set_ylabel(r'$t$ [m]', fontdict = font)
    ax1.set_xticks(xtick_pos)
    ax1.set_xticklabels(xtick_lab2, fontdict = font)
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.set_ylim(0,nt-1)
    
    ax1.imshow(cube,
               cmap = 'gray', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               vmin=vmin,
               vmax = vmax,
               #aspect = 0.45,
    )

    #ax1.plot(yi, ti, color = 'yellow', alpha = 0.75, linewidth = 0.5)
    #ax1.plot(test_func(xi, params[0], params[1], params[2], params[3]), xi, alpha = 0.75, linestyle = '--')
    #ax1.plot(smooth(yi,20), xi, color = 'red', alpha = 0.75)



    # GUI stuff
    click_fit = 1
    
    if(click_fit):
        
        cursor = Cursor(ax1, horizOn=True, vertOn=True, useblit=True,
                        color = 'r', linewidth = 0.75, alpha = 0.75)
        
        coord = []
        
        def onclick(event):
            xx, tt = event.xdata, event.ydata
            global coord
            coord.append((xx, tt))
            print('a(t)=%1.3f, t=%1.2f' %
                  (np.round(xx*res, decimals = 3),
                   np.round(tt*cad, decimals = 2)))
            f.canvas.draw() #redraw the figure
        

        #f.canvas.draw()
        cid = f.canvas.mpl_connect('button_press_event', onclick)


        plt.tight_layout()
        plt.show()

        stop()

        coord = np.array(coord)
        
        np.savetxt(cut_fname, coord, fmt='%3.8f', delimiter=' ', newline='\n', header='a(t) [arcsec], t [s]', footer='', comments='# ', encoding=None)
        
        

    else:
        coord = np.loadtxt(cut_fname, dtype = 'float64')

    amin, amax = int(np.min(coord[:,0])), int(np.max(coord[:,0]))+1
    tmin, tmax = int(np.min(coord[:,1])), int(np.max(coord[:,1]))+1
    #arange = np.linspace(amin, amax, amax-amin+1)
    trange = np.linspace(tmin, tmax, tmax-tmin+1)

    # interpolating in-between points of the click_fit coord
    y_interp = np.interp(ti, coord[:,1], coord[:,0])
    fname = '/home/seki2695/OUTPUT/project2/172118_interp.txt'
    np.savetxt(fname, coord, fmt='%3.8f', encoding=None)
    
    # track max intensity based on the interp click_fit
    y_imax = np.zeros(nt)
    dist = 1. # width of the interval in which the max intensity is taken
    for i in range(nt):
        y_imax[i] = np.argmax(cube[i, int(y_interp[i]-dist):int(y_interp[i]+dist)])+int(y_interp[i])-dist
        #print( int(y_interp[i]-dist), int(y_interp[i]+dist))


    # extracting inversion results
    inv_file = file_search(invdir, '*atmosout*'+cut_file[-11:-5]+'*nc')[0]
    inv_res = sp.model(invdir+inv_file)

    temp_cube = inv_res.temp[0]
    vlos_cube = inv_res.vlos[0]
    vturb = inv_res.vturb[0]
    ltau = inv_res.ltau[0,0,0,:]
    ntau = ltau.size
    
    tau_i = 23 

    # intensity oscillation PLOTs
    plt.close('all')
    fig = plt.figure(figsize = (8,5.5))
    gs = gridspec.GridSpec(3,1)

    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0])
    ax3 = plt.subplot(gs[2,0])
    
    ax1.set_ylabel(r'length [arcsec]', fontdict = font)
    ax1.set_xticks(ttick_pos)
    ax1.set_xticklabels([])
    ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.set_xlabel('')
    ax1.set_yticks(xtick_pos)
    ax1.set_yticklabels(xtick_lab2, fontdict = font)
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.set_xlim(0,nt-1)
    
    ax2.set_ylabel(r'length [arcsec]', fontdict = font)
    ax2.set_xticks(ttick_pos)
    ax2.set_xticklabels([])
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.set_xlabel('')
    ax2.set_yticks(xtick_pos)
    ax2.set_yticklabels(xtick_lab2, fontdict = font)
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.set_xlim(0,nt-1)

    ax3.set_ylabel(r'length [arcsec]', fontdict = font)
    ax3.set_xticks(ttick_pos)
    ax3.set_xticklabels(ttick_lab, fontdict = font)
    ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax3.set_xlabel(r'$t$ [m]', fontdict = font)
    ax3.set_yticks(xtick_pos)
    ax3.set_yticklabels(xtick_lab2, fontdict = font)
    ax3.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax3.set_xlim(0,nt-1)
    
    smth_y = smooth(y_imax, 20)
    fname = '/home/seki2695/OUTPUT/project2/172118_final.txt'
    np.savetxt(fname, smth_y, fmt='%3.8f', encoding=None)


    ax1.imshow(np.transpose(cube),
               cmap = 'gray', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               vmin=2.e-9,
               vmax = 3.5e-9,
               #aspect = 0.45,
    )
    #ax1.plot(coord[:,1], coord[:,0],color = 'red', alpha = 0.75, label = 'Manual fit')
    ax1.plot(ti, y_imax, color = 'yellow', alpha = 0.75, label = r'I$_{\rm{max}}$')
    ax1.plot(ti, smth_y, color = 'red', alpha = 0.75, label = r'smthd I$_{\rm{max}}$')

    ax2.imshow(np.transpose(temp_cube[:,:,tau_i]),
               cmap = 'inferno', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               vmin=4500,
               #aspect = 0.45,
    )
    #ax2.plot(ti, smth_y, color = 'red', alpha = 0.75, label = r'smthd I$_{\rm{max}}$')

    ax3.imshow(np.transpose(vlos_cube[:,:,tau_i]),
               cmap = 'bwr', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               vmin=-8.e5,
               vmax = 8.e5,
               #aspect = 0.45,
    )
    #ax3.plot(ti, smth_y, color = 'red', alpha = 0.75, label = r'smthd I$_{\rm{max}}$')

    # smooth fit to the oscillations
    Imax_smth = np.zeros(nt)
    temp_smth = np.zeros((nt,ntau))
    vlos_smth = np.zeros((nt,ntau))
    for i in range(nt):
        Imax_smth[i] = cube[i,int(np.rint(smth_y[i]))]
        for j in range(ntau):
            temp_smth[i,j] = temp_cube[i,int(np.rint(smth_y[i])),j]
            vlos_smth[i,j] = vlos_cube[i,int(np.rint(smth_y[i])),j]

    plt.subplots_adjust(left=0.05,
                    bottom=0.08, 
                    right=0.99, 
                    top=0.99, 
                    wspace=0.22, 
                    hspace=0.0
    )
    plt.show()
    fname = outdir + 'roi_invres_curvefit_' +cut_file[-11:-5]+'.pdf'
    fig.savefig(fname, quality = 100)

    #stop()
    
    #################
    # physical oscillation PLOTs
    #################
    plt.close('all')
    f = plt.figure(figsize = (8,8))
    gs = gridspec.GridSpec(6,2)
    
    # ltau pos
    upper, lower = 43, 12
    tautick_pos = np.linspace(lower,upper, 5)
    tautick_lab = np.int8(-np.linspace(int(np.abs(ltau[upper])), int(np.abs(ltau[lower])), int(ltau[upper])-int(ltau[lower])+1))
    
    ax0, ax1 ,ax2, ax3, ax4, ax5 = plt.subplot(gs[0,0]), plt.subplot(gs[0,1]), plt.subplot(gs[1:-1,0]), plt.subplot(gs[1:-1,1]), plt.subplot(gs[-1,0]), plt.subplot(gs[-1,1])
    
    ax0.set_ylabel(r'length [arcsec]', fontdict = font)
    ax0.set_xticks(ttick_pos)
    ax0.set_xticklabels([])
    ax0.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax0.set_xlabel('')
    ax0.set_yticks(xtick_pos)
    ax0.set_yticklabels(xtick_lab2, fontdict = font)
    ax0.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax0.set_xlim(0,nt-1)
    ax0.set_title(r'T')
    
    ax1.set_ylabel('')
    ax1.set_xticks(ttick_pos)
    ax1.set_xticklabels([])
    ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.set_xlabel('')
    ax1.set_yticks(xtick_pos)
    ax1.set_yticklabels([])
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.set_xlim(0,nt-1)
    ax1.set_title(r'v$\rm{_{LOS}}$')
    
    ax2.set_ylabel(r'log($\tau_{500}$)', fontdict = font)
    ax2.set_xticks(ttick_pos)
    ax2.set_xticklabels([])
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.set_xlabel('')
    ax2.set_yticks(tautick_pos)
    ax2.set_yticklabels(tautick_lab, fontdict = font)
    ax2.set_ylim(11,44.5)
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.set_xlim(0,nt-1)
    
    ax3.set_ylabel('')
    ax3.set_xticks(ttick_pos)
    ax3.set_xticklabels([])
    ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax3.set_xlabel('')
    ax3.set_yticks(tautick_pos)
    ax3.set_yticklabels([])
    ax3.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax3.set_xlim(0,nt-1)
    ax3.set_ylim(11,44.5)

    ax4.set_ylabel(r'$T$ [$kK$]', fontdict = font)
    ax4.set_xticks(ttick_pos)
    ax4.set_xticklabels(ttick_lab, fontdict = font)
    ax4.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax4.set_xlabel(r'$t$ [m]', fontdict = font)
    #ax4.set_yticks(tautick_pos)
    #ax4.set_yticklabels(tautick_lab, fontdict = font)
    #ax4.set_ylim(11,44.5)
    ax4.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax4.set_xlim(0,nt-1)
    ax4.tick_params(labelsize = 8)
    
    ax5.set_ylabel(r'$v_{\rm{LOS}}$ [$km$ $s^{-1}$]', fontdict = font)
    ax5.set_xticks(ttick_pos)
    ax5.set_xticklabels(ttick_lab, fontdict = font)
    ax5.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax5.set_xlabel(r'$t$ [m]', fontdict = font)
    #ax5.set_yticks(tautick_pos)
    #ax5.set_yticklabels([])
    ax5.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax5.set_xlim(0,nt-1)
    #ax5.set_ylim(11,44.5)
    ax5.tick_params(labelsize = 8)

    ax0.imshow(np.transpose(temp_cube[:,:,tau_i]),
               cmap = 'inferno',
               origin = 'lower', #aspect=ratio
               vmin=4500,
               #vmax = 3.5e-9,
               #aspect = 0.45,
    )
    ax0.plot(coord[:,1], coord[:,0],color = 'red', alpha = 0.75, label = 'Manual fit')
    ax0.plot(ti, y_imax, color = 'yellow', alpha = 0.75, label = r'I$_{\rm{max}}$')
    ax0.plot(ti, smth_y, color = 'cyan', alpha = 0.75, label = r'smthd I$_{\rm{max}}$')
    
    ax1.imshow(np.transpose(vlos_cube[:,:,tau_i]),
               cmap = 'bwr',
               origin = 'lower', #aspect=ratio
               vmin=-8e5,
               vmax = 8e5,
               #aspect = 0.45,
    )
    ax1.plot(coord[:,1], coord[:,0],color = 'red', alpha = 0.75, label = 'fit')
    ax1.plot(ti, y_imax, color = 'yellow', alpha = 0.75, label = r'I$_{\rm{max}}$')
    ax1.plot(ti, smth_y, color = 'cyan', alpha = 0.75, label = r'I$_{\rm{max, sm}$')

    temp_plottingfactor, vlos_plottingfactor= -6, -1
    for j in range(lower,upper):
        ax2.plot(ti, temp_smth[:,j]/1000+upper+lower-j +temp_plottingfactor, color = 'black', alpha = 1-j*0.02, linewidth = 1.)
        ax3.plot(ti, vlos_smth[:,j]/4.5e5+upper+lower-j+vlos_plottingfactor, color = 'black', alpha = 1-j*0.016, linewidth = 1.)
        if(j==23):
            ax2.plot(ti, smooth(temp_smth[:,j]/1000+upper+lower-j +temp_plottingfactor,20), color = 'blue', alpha = 1-j*0.02)
            ax3.plot(ti, smooth(vlos_smth[:,j]/4.5e5+upper+lower-j+vlos_plottingfactor,20), color = 'blue', alpha = 1-j*0.016)
        else:
            ax2.plot(ti, smooth(temp_smth[:,j]/1000+upper+lower-j +temp_plottingfactor,20), color = 'red', alpha = 1-j*0.02)
            ax3.plot(ti, smooth(vlos_smth[:,j]/4.5e5+upper+lower-j+vlos_plottingfactor,20), color = 'red', alpha = 1-j*0.016)
        
    ax4.plot(temp_smth[:,tau_i]/1.e3, color = 'black', alpha = 0.75, linewidth = 1.)
    ax4.plot(smooth(temp_smth[:,tau_i]/1.e3, 20), color = 'blue', alpha = 0.75)
    ax5.plot(vlos_smth[:,tau_i]/1.e5, color = 'black', alpha = 0.75, linewidth = 1.)
    ax5.plot(smooth(vlos_smth[:,tau_i]/1.e5, 20), color = 'blue', alpha = 0.75)

    # Calculation of the velocity of the oscillation in the plane of sky
    xp = trange*cad
    fxp = smth_y[tmin:tmax+1]*res*725
    der = derivative(xp, fxp, dx = 1e-5)
    x, fx, f1x = der[0], der[1], der[2]
    ax5.plot(x/cad, f1x, color = "orangered")
    
    plt.subplots_adjust(left=0.08,
                    bottom=0.05, 
                    right=0.99, 
                    top=0.99, 
                    wspace=0.22, 
                    hspace=0.10
    )


    plt.show()
    filename = outdir + 'oscillation_invres_curvefit_' +cut_file[-11:-5]+'.pdf'
    f.savefig(filename, quality = 100)
    print 'file saved to: ' + filename
    
