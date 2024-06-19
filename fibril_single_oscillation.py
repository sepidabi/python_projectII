# fibril_single_oscillation.py

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

# Font properties
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }

# Functions
def test_func(x, a, b, c, d):
    return a * np.sin(b * x - c)+d

#boxcar smoothing
#def smooth(y, w):
#    N = y.shape[0]
 #   r = np.zeros(N)
  #  for i in range(N):
   #     if(i==0 or i==N-1):
    #        r[i] = y[i]
     #   elif(i>(w-1.)/2. and i<=N-(w+1.)/2.):
      #      r[i] = np.average(y[int(np.rint(i-w/2.)) : int(np.rint(i+w/2.-1))])
       # else:
        #    r[i] = (y[i-1]+y[i]+y[i+1])/3.
    #return r


# boxcar smoothing
"""def smooth(y, w):
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
"""

def smooth(x,window_len=11,window='flat'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

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
    
# Directories
datadir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
outdir = '/home/seki2695/OUTPUT/project2/'
invdir = '/home/seki2695/INV/stic/II/slabs_tseries/'
savedir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/w_calibrated/'

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

cut1 = restore(fibdir+cut1_file)
cut2 = restore(fibdir+cut2_file)
cut3 = restore(fibdir+cut3_file)

#######################
## plot oscilations
#######################

plt.close('all')
f = plt.figure(figsize = (8, 3.1))

# axis info
top = 0.88
bottom = 0.1
left = 0.05
right = 0.99
wspace = 0.01
hspace = 0.001
dx = (right - left-2.*wspace)/3.

ratio = 1.25

gs = gridspec.GridSpec(ncols = 1, nrows = 3, figure = f) # grid scale of the right col.
gs.update(left=left,
          right=left+dx,
          wspace=wspace,
          bottom=bottom,
          top=top,
          hspace = hspace
)

gs1 = gridspec.GridSpec(ncols = 1, nrows = 3, figure = f) # grid scale of the right col.
gs1.update(left=left+dx+wspace,
          right=left+2*dx+wspace,
          wspace=wspace,
          bottom=bottom,
          top=top,
          hspace = hspace
)

gs2 = gridspec.GridSpec(ncols = 1, nrows = 3, figure = f) # grid scale of the right col.
gs2.update(left=left+2*dx+2*wspace,
          right=left+3*dx+2*wspace,
          wspace=wspace,
          bottom=bottom,
          top=top,
          hspace = hspace
)

rr = 40/1.5
xtick_pos = np.arange(0,(np.round((64)/rr)+1)*rr,rr)
tt = 30
ttick_pos = np.arange(0,(np.round((213)/tt)+1)*tt,tt)
ttick_lab = ttick_pos*cad/60
cuts = [cut3, cut2, cut1]
cuts_file = [cut3_file, cut2_file, cut1_file]
cut_no = ['(1)', '(2)', '(3)']

##############
# 2nd plot settings
##############
ff = plt.figure(figsize = (4,6.))
ggs = gridspec.GridSpec(ncols =2, nrows = 3, figure = ff)
ggs.update(left = 0.12,
           bottom = 0.07,
           right = 0.97,
           top = 0.96,
           wspace = 0.07,
           hspace = 0.065
)

#############
#  3rd plot settings
#############
fff = plt.figure(figsize = (4,4))
gggs = gridspec.GridSpec(ncols = 1, nrows = 3, figure = fff)
gggs.update(left = 0.12,
           bottom = 0.1,
           right = 0.88,
           top = 0.92,
           wspace = 0.0,
           hspace = 0.1
)

for gg in range(3):
    cut = cuts[gg]
    cut_file = cuts_file[gg]
    cut_fname = outdir+cut_file[-11:-5]+'.txt'
    cube_cali_file = savedir+file_search(savedir,'slab_obs3950_31_'+cut_file[-11:-5]+'.fits')[0]
    cube_cali = (mf.readfits(cube_cali_file)[1,w_pos,:,:]).squeeze()*1e6
    cube_raw = np.mean(cuts[gg].loop_slab[w_pos-3:w_pos+3,:,:],axis = 0).squeeze()*1e9
    # background intensity from the wings
    wg = 0.15
    cube_bg = wg*np.mean(cut.loop_slab[:w_pos-6,:,:]+cut.loop_slab[w_pos+8:,:,:],axis=0).squeeze()*1e9
    cube_trunc = cube_raw - cube_bg # removing the background intensity
    # contrast correction
    cube_med = exposure.rescale_intensity(sgnl.medfilt2d(cube_trunc,kernel_size = [3,1]), out_range=(0, 1.))
    cube_sharp = exposure.rescale_intensity(sharpen(cube_trunc, sigma =[3,1], unsigma = [1,3]), out_range=(0, 1.))
    cube_final = exposure.rescale_intensity(exposure.adjust_gamma(cube_sharp,0.1), out_range=(0, 1.))
    #vmin, vmax = 0.9, 0.98

    cube = cube_cali
    coord = np.loadtxt(cut_fname, dtype = 'float64')
    amin, amax = int(np.min(coord[:,0])), int(np.max(coord[:,0]))+1
    tmin, tmax = int(np.min(coord[:,1])), int(np.max(coord[:,1]))+1
    trange = np.linspace(tmin, tmax, tmax-tmin+1)
    xlim_max = 176#int(np.rint(np.max(coord[:,1])))
    ylim_min = 4
    ylim_max = -1

    # interpolating in-between points of the click_fit coord
    ti = np.arange(0, nt, 1.)
    yi = np.argmax(cube, axis = 1)
    y_interp = np.interp(ti, coord[:,1], coord[:,0])
    
    # track max intensity based on the interp click_fit
    y_imax = np.zeros(nt)
    dist = 2. # width of the interval in which the max intensity is taken
    sm_fact = 8
    for i in range(nt):
        y_imax[i] = np.argmax(cube[i, int(y_interp[i]-dist):int(y_interp[i]+dist)])+int(y_interp[i])-dist
    smth_y = smooth(y_imax, sm_fact)

    # extracting inversion results
    inv_file = file_search(invdir, '*atmosout*'+cut_file[-11:-5]+'*nc')[0]
    inv_res = sp.model(invdir+inv_file)
    temp_cube = inv_res.temp[0].squeeze()*1e-3
    vlos_cube = inv_res.vlos[0].squeeze()*1e-5
    vturb = inv_res.vturb[0]
    ltau = inv_res.ltau[0,0,0,:]
    ntau = ltau.size
    tau_i = 23
    
    ##########
    # 1st plot axes
    ##########
    ax0 = f.add_subplot(gs[gg,0],adjustable = 'box')#,aspect = 'equal')
    ax1 = f.add_subplot(gs1[gg,0],adjustable = 'box')#,aspect = 'equal')
    ax2 = f.add_subplot(gs2[gg,0],adjustable = 'box')#,aspect = 'equal')
    xticklabs = np.intc(np.round(xtick_pos*res*(np.sqrt((cuts[gg].x_coords[0]-cuts[gg].x_coords[1])**2 +(cuts[gg].y_coords[0]-cuts[gg].y_coords[1])**2)/cuts[gg].loop_size),decimals = 0))
    if  (gg==1):
        ax0.set_ylabel(r'length [arcsec]', fontdict = font)
    ax0.set_xticks(ttick_pos)
    ax1.set_xticks(ttick_pos)
    ax2.set_xticks(ttick_pos)
    if (gg==2):
        ax0.set_xticklabels(ttick_lab, fontdict = font)
        ax0.set_xlabel(r'$t$ [m]', fontdict = font)
        ax1.set_xticklabels(ttick_lab, fontdict = font)
        ax1.set_xlabel(r'$t$ [m]', fontdict = font)
        ax2.set_xticklabels(ttick_lab, fontdict = font)
        ax2.set_xlabel(r'$t$ [m]', fontdict = font)
    else:
        ax0.set_xticklabels([])
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
    ax0.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax0.set_yticks(xtick_pos)
    ax0.set_yticklabels(xticklabs, fontdict = font)
    ax0.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax0.set_xlim(0,xlim_max)
    ax0.set_ylim(0,cube.shape[1]-abs(ylim_min)-abs(ylim_max)-1)
    ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.set_yticks(xtick_pos)
    ax1.set_yticklabels([])
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.set_xlim(0,xlim_max)
    ax1.set_ylim(0,cube.shape[1]-abs(ylim_min)-abs(ylim_max)-1)
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.set_yticks(xtick_pos)
    ax2.set_yticklabels([])
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.set_xlim(0,xlim_max)
    ax2.set_ylim(0,cube.shape[1]-abs(ylim_min)-abs(ylim_max)-1)

    cub = (np.transpose(cube_final))[ylim_min:ylim_max,:xlim_max]
    fact = 3.
    vmin = np.mean(cub) - fact*np.std(cub)
    vmax = np.mean(cub) + fact*np.std(cub)

    # 1st column
    im = ax0.imshow(cub, cmap = 'gray', origin = 'lower', aspect=ratio,
              vmin = vmin,#2.e-9,
              vmax=vmax#3.5e-9
    )

    #if(gg==0):
     #   for v in range(len(fr_select)):
      #      ax0.axvline(fr_select[v], ymin = 0, ymax = 10, linewidth = 0.75, color = 'orangered', linestyle = '--', alpha = 0.5)
    ax0.text(2,2,cut_no[gg], color = 'white')
    ax0.plot(ti, y_imax-ylim_min, color = 'ivory', alpha = 0.75, label = r'I$_{\rm{max}}$', linestyle = '--', linewidth = 1.)
    ax0.scatter(x = ti[fr_select], y = (smth_y-ylim_min)[fr_select], marker = "|", color = 'orangered', zorder = 10)
    ax0.plot(ti, smth_y-ylim_min, color = 'orangered', alpha = 0.75, label = r'smthd I$_{\rm{max}}$')
    np.save('cut'+str(gg)+'.npy', (y_imax)[fr_select]/cube.shape[1])
    
    # 2nd column
    temp_cub = exposure.rescale_intensity(sharpen(np.transpose(temp_cube[:,:,tau_i])[ylim_min:ylim_max,:xlim_max], sigma =[1,3], unsigma = [3,1]), out_range=(0, 1.))

    im1 = ax1.imshow(temp_cub, cmap = 'inferno', origin = 'lower', aspect=ratio,
              #vmin = 4.1,
              #vmax=5.4#3.5e-9
    )
    ax1.plot(ti, smth_y-ylim_min, color = 'white', alpha = 0.5, label = r'smthd I$_{\rm{max}}$')

    # 3rd column
    sigma = 50
    crap = spnd.filters.gaussian_filter(vlos_cube[:,:,tau_i], [sigma, sigma], mode = 'constant')
    vlos_cub = (np.transpose(vlos_cube[:,:,tau_i]-crap))[ylim_min:ylim_max,:xlim_max]
    im2 = ax2.imshow(vlos_cub, cmap = 'bwr', origin = 'lower', aspect=ratio,
              vmin = -8.31,
              vmax=8.31
    )
    ax2.plot(ti, smth_y-ylim_min, color = 'black', alpha = 0.5, label = r'smthd I$_{\rm{max}}$')

    if (gg==0):
        axin = inset_axes(ax0,
                          width="90%",  # width = 10% of parent_bbox width
                          height="7%",  # height : 50%
                          loc='lower center',
                          bbox_to_anchor=(0, 1.03, 1, 1),
                          bbox_transform=ax0.transAxes,
                          borderpad=0,
        )
        axin1 = inset_axes(ax1,
                           width="90%",  # width = 10% of parent_bbox width
                           height="7%",  # height : 50%
                           loc='lower center',
                           bbox_to_anchor=(0, 1.03, 1, 1),
                           bbox_transform=ax1.transAxes,
                           borderpad=0,
        )
        axin2 = inset_axes(ax2,
                           width="90%",  # width = 10% of parent_bbox width
                           height="7%",  # height : 50%
                           loc='lower center',
                           bbox_to_anchor=(0, 1.03, 1, 1),
                           bbox_transform=ax2.transAxes,
                           borderpad=0,
        )

        cb_ticks = [[1.5,2,2.5,3,3.5],
                    [4.3, 4.6, 4.9, 5.2],
                    [-6,-3, 0,3,6]]
        cb_tick_labels = [str(cb_ticks[0][0]), str(cb_ticks[0][1]), str(cb_ticks[0][2]), str(cb_ticks[0][3]),  str(cb_ticks[0][4])]
        cb_tick_labels1 = [str(cb_ticks[1][0]), str(cb_ticks[1][1]), str(cb_ticks[1][2]), str(cb_ticks[1][3]), str(cb_ticks[1][1])]
        cb_tick_labels2 = [str(cb_ticks[2][0]), str(cb_ticks[2][1]), str(cb_ticks[2][2]), str(cb_ticks[2][3]), str(cb_ticks[2][4])]

        # colorbar intensity
        cb = f.colorbar(im, cax = axin, orientation = "horizontal", ticks = cb_ticks[0])
        cb.ax.set_xticklabels(cb_tick_labels, fontdict = font)
        cb.ax.xaxis.tick_top()
        cb.ax.xaxis.set_tick_params(pad = 0)
        cb.ax.set_title(r'I [$\times $10$^{3}$ W m$^{-2}$ Hz$^{-1}$ sr$^{-1}$]', fontdict = font, pad = 17)

        # colorbar Temperature
        cb1 = f.colorbar(im1, cax = axin1, orientation = "horizontal", ticks = cb_ticks[1])
        cb1.ax.set_xticklabels(cb_tick_labels1, fontdict = font)
        cb1.ax.xaxis.tick_top()
        cb1.ax.xaxis.set_tick_params(pad = 0)
        cb1.ax.set_title(r'T [kK]', fontdict = font, pad = 17)

        # colorbar Temperature
        cb2 = f.colorbar(im2, cax = axin2, orientation = "horizontal", ticks = cb_ticks[2])
        cb2.ax.set_xticklabels(cb_tick_labels2, fontdict = font)
        cb2.ax.xaxis.tick_top()
        cb2.ax.xaxis.set_tick_params(pad = 0)
        cb2.ax.set_title(r'$v_{\rm{LOS}}$ [km s$^{-1}$]', fontdict = font, pad = 17)

        
    #plt.show(f)
    #stop()
    ##########
    # 2nd plot axes
    ##########
    axx0 = ff.add_subplot(ggs[gg,0])
    axx1 = ff.add_subplot(ggs[gg,1])

    # ltau pos
    upper, lower = 58, 0
    tautick_pos = 58-np.flip(np.array([12,18,25,32]),0)#np.arange(lower,upper,6)#np.where(np.round(ltau, decimals = 1)%1==0.)[0]#np.linspace(lower,upper, 10)
    tautick_lab = ['-3', '-4', '-5', '-6']#np.flip(np.int8(np.round(ltau[tautick_pos])),0)#np.flip(np.int8(np.round(ltau[tautick_pos])),0)#np.int8(-np.linspace(int(np.abs(ltau[upper])), int(np.abs(ltau[lower])), int(ltau[upper])-int(ltau[lower])+1))

    if (gg==0):
        axx0.set_title(r'$T$')
    axx0.set_xticks(ttick_pos)
    axx0.xaxis.set_minor_locator(AutoMinorLocator(4))
    axx0.set_yticks(tautick_pos)
    axx0.set_ylabel(r'log($\tau_{500}$)', fontdict = font)
    axx0.set_yticklabels(tautick_lab, fontdict = font)
    if (gg==2):
        axx0.set_xlabel(r'$t$ [m]', fontdict = font)
        axx0.set_xticklabels(ttick_lab, fontdict = font)
    else:
        axx0.set_xticklabels([])
        axx0.set_xlabel('')
    axx0.set_ylim(20,47)
    axx0.yaxis.set_minor_locator(AutoMinorLocator(5))
    axx0.set_xlim(0,xlim_max)
    
    if (gg==0):
        axx1.set_title(r'$v\rm{_{LOS}}$')
    axx1.set_xticks(ttick_pos)
    axx1.xaxis.set_minor_locator(AutoMinorLocator(4))
    axx1.set_xlabel(r'$t$ [m]', fontdict = font)
    axx1.set_yticks(tautick_pos)
    if (gg==2):
        axx1.set_xticklabels(ttick_lab, fontdict = font)
        axx1.set_xlabel(r'$t$ [m]', fontdict = font)
    else:
        axx1.set_xticklabels([])
        axx1.set_xlabel('')        
    axx1.set_yticklabels([])
    axx1.set_ylabel('')
    axx1.set_ylim(20,47)
    axx1.yaxis.set_minor_locator(AutoMinorLocator(5))
    axx1.set_xlim(0,xlim_max)


    # CURVE FITTING

    # smooth fit to the oscillations
    #Imax_smth = np.zeros(nt)
    temp_smth = np.zeros((nt,ntau))
    vlos_smth = np.zeros((nt,ntau))
    for i in range(nt):
        #Imax_smth[i] = cube[i,int(np.rint(smth_y[i]))]
        for j in range(ntau):
            temp_smth[i,j] = temp_cube[i,int(np.rint(smth_y[i])),j]
            vlos_smth[i,j] = vlos_cube[i,int(np.rint(smth_y[i])),j]
            

    temp_plottingfactor, vlos_plottingfactor= 0,0#-6, -1
    v_boost = 2.5
    t_boost = 0.5
    for j in range(lower,upper):
        if (j%2==1 and j>11):# and j<30):
            axx0.plot(ti, (temp_smth[:,j]-np.mean(temp_smth[:,j]))/t_boost+upper+lower-j +temp_plottingfactor, color = 'black', alpha = 1-j*0.016, linewidth = 1., linestyle = ':')
            axx1.plot(ti, (vlos_smth[:,j]-np.mean(vlos_smth[:,j]))/v_boost+upper+lower-j+vlos_plottingfactor, color = 'black', alpha = 1-j*0.016, linewidth = 1., linestyle = ':')
            if(j==23):
                axx0.plot(ti, smooth((temp_smth[:,j]-np.mean(temp_smth[:,j]))/t_boost+upper+lower-j +temp_plottingfactor,8), color = 'red', alpha = 1-j*0.016)
                axx1.plot(ti, smooth((vlos_smth[:,j]-np.mean(vlos_smth[:,j]))/v_boost+upper+lower-j+vlos_plottingfactor,8), color = 'red', alpha = 1-j*0.016)
            else:
                axx0.plot(ti, smooth((temp_smth[:,j]-np.mean(temp_smth[:,j]))/t_boost+upper+lower-j +temp_plottingfactor,8), color = 'black', alpha = 1-j*0.016)
                axx1.plot(ti, smooth((vlos_smth[:,j]-np.mean(vlos_smth[:,j]))/v_boost+upper+lower-j+vlos_plottingfactor,8), color = 'black', alpha = 1-j*0.016)
    axx1.text(0.86,0.92, cut_no[gg], color='black', transform=axx1.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha = 0.75, pad = 2))

    
    #############
    # 3rd plot axes
    #############
    axxx0 = fff.add_subplot(gggs[gg,0])
    axxx1 = axxx0.twinx()

    axxx0.set_ylabel(r'[kK]', fontdict = font)
    axxx0.set_xticks(ttick_pos)
    axxx0.xaxis.set_minor_locator(AutoMinorLocator(4))
    #axxx0.set_yticks(tautick_pos)
    #axxx0.set_yticklabels(tautick_lab, fontdict = font)
    #axxx0.set_ylim(11,44.5)
    axxx0.yaxis.set_minor_locator(AutoMinorLocator(2))
    axxx0.set_xlim(0,xlim_max)
    axxx0.tick_params(labelsize = 8)
    axxx0.set_ylim([4.3,5.3])
    axxx0.set_yticks([4.4,4.6,4.8,5.,5.2])
    
    axxx1.set_ylabel(r'[km s$^{-1}$]', fontdict = font)
    #axxx1.set_xticks(ttick_pos)
    #axxx1.xaxis.set_minor_locator(AutoMinorLocator(4))
    #axxx1.set_yticks(tautick_pos)
    #axxx1.set_yticklabels([])
    axxx1.yaxis.set_minor_locator(AutoMinorLocator(2))
    axxx1.set_xlim(0,xlim_max)
    #axxx1.set_ylim(11,44.5)
    axxx1.tick_params(labelsize = 8)
    axxx1.set_ylim([-5,5])
    axxx1.set_ylim([4.3,5.3])
    axxx1.set_yticks([-4,-2,0,2,4])

    if gg==2:
        axxx0.set_xticklabels(ttick_lab, fontdict = font)
        axxx0.set_xlabel(r'$t$ [m]', fontdict = font)
        #axxx1.set_xticklabels(ttick_lab, fontdict = font)
        #axxx1.set_xlabel(r'$t$ [min]', fontdict = font)
    else:
        axxx0.set_xticklabels([])
        axxx0.set_xlabel('')
        #axxx1.set_xticklabels([])
        #axxx1.set_xlabel('')

    # Calculation of the velocity of the oscillation in the plane of sky
    tmin, tmax = 2,177
    trange = np.linspace(tmin, tmax, tmax-tmin+1)
    xp = trange*cad
    fxp = smth_y[tmin:tmax+1]*res*725
    der = derivative(xp, fxp, dx = 1e-5)
    x, fx, f1x = der[0]-np.min(der[0]), der[1], der[2]
    per = vlos_smth.shape[0]
    x_brief = np.zeros(per)
    fx_brief = np.zeros(per)
    f1x_brief = np.zeros(per)
    brief_fact = f1x.shape[0]/per
    sm_fact = 8

    for iii in range(per):
        x_brief[iii] = x[iii*brief_fact]
        fx_brief[iii] = fx[iii*brief_fact]
        f1x_brief[iii] = f1x[iii*brief_fact]
        

    axxx0.plot(temp_smth[:,tau_i]/1., color = 'gold', linestyle = ':', alpha = 0.25)
    axxx1.plot(vlos_smth[:,tau_i]/1., color = 'orangered', linestyle = ':', alpha = 0.25)
    axxx1.plot(x_brief/cad, f1x_brief-np.mean(f1x_brief), color = "gray", linestyle = ':', alpha = 0.25)
    axxx0.plot(smooth(temp_smth[:,tau_i]/1., 8), color = 'gold', alpha = 0.75, label = r'$T$')
    axxx1.plot(smooth(vlos_smth[:,tau_i]/1., 8), color = 'orangered', alpha = 0.75, label = r'$v_{\rm{LOS}}$')
    print(x_brief[-1], tmin, tmax)
    axxx1.axhline(0, linestyle = ":", color = "black", linewidth = 1, alpha = 0.5)
    axxx1.plot(x_brief/cad, smooth(f1x_brief-np.mean(f1x_brief),sm_fact), color = "gray", alpha = 0.75, label = r'$v_{\rm{POS}}$')
    if gg==0:
        fff.legend(loc='upper center', fontsize = 8., ncol = 3)#, handlelength = 0.5, borderpad=0., handletextpad = 0.2, labelspacing = 0.
    
    axxx0.text(0.9,0.85,cut_no[gg], color = 'black', transform=axxx0.transAxes, size = 10)#, rotation = 90)


# file saving 1
if (w_pos == 13): wpos='K3'
elif (w_pos==10): wpos = 'K2V'
else: wpos = 'K2R'
filename = outdir+'oscilation'+wpos+'_invres.pdf'
#stop()
f.savefig(filename, quality = 100)
#f.tight_layout()
print 'file saved to: ' + filename
#plt.close(f)

# file saving 2
filename = outdir+'oscilation_seq_invres.pdf'
ff.savefig(filename, quality = 100)
#ff.tight_layout()
print 'file saved to: ' + filename
#stop()
#plt.close(ff)

# file saving 3
filename = outdir + 'oscilation_single_invres.pdf'
fff.savefig(filename, quality = 100)
#fff.tight_layout()
print 'file saved to: ' + filename
plt.show()
stop()
#plt.close(fff)
#plt.close("all")

        
#################
# physical oscillation PLOTs
#################

        
ax4.plot(temp_smth[:,tau_i]/1., color = 'black', alpha = 0.75, linewidth = 1.)
ax4.plot(smooth(temp_smth[:,tau_i]/1., 8), color = 'blue', alpha = 0.75)
ax5.plot(vlos_smth[:,tau_i]/1., color = 'black', alpha = 0.75, linewidth = 1.)
ax5.plot(smooth(vlos_smth[:,tau_i]/1., 8), color = 'blue', alpha = 0.75)

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

