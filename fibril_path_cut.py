# The following script is meant to read the firbil path
# and the cuts across to record its oscillations

import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import seaborn as sns
import scipy.ndimage as spnd
import matplotlib
from matplotlib.gridspec import GridSpec

# fonts
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }

# Declerations
outdir = '/home/seki2695/OUTPUT/project2/'
datadir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
savedir = datadir+'OUTPUT/'
ncut = 4 # number of cutes per fibril
res = 0.0375 # CHROMIS pixel size in arcse
tres = 8./60 #CHROMIS cadence

# Extracting the data cubes
res = 0.0375 # CHROMIS pixel size in arcsec
pref = ['6173','8542','3950']#,'6563']
pref_name = ['fe','ca8','cak']#,'ha']
file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
file_ca8 = file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak = file_search(datadir,'crispex*'+pref[2]+'*.fcube')
#file_ha = file_search(datadir,'crispex*'+pref[3]+'*.fcube')
cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])#[fr]
cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])#[fr]
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])#[fr]
#cube_ha = lp_read(datadir+file_ha[0],datadir+file_ha[1])[fr]
cak_cont = cube_cak[0,-1,:,:]
fe_cont = cube_fe[0,-1,:,:]

# Extracting the fibrils
fibdir = datadir+'paths/'
n = len(file_search(fibdir, 'crispex*3950*.csav'))/(ncut+1) # total number of the fibrils

for i in range (n):
    fib_file = restore(fibdir+(file_search(fibdir, 'crispex*3950*.csav'))[(ncut+1)*i])
    cut1_file = restore(fibdir+(file_search(fibdir, 'crispex*3950*.csav'))[(ncut+1)*i+1])
    cut2_file = restore(fibdir+(file_search(fibdir, 'crispex*3950*.csav'))[(ncut+1)*i+2])
    cut3_file = restore(fibdir+(file_search(fibdir, 'crispex*3950*.csav'))[(ncut+1)*i+3])
    cut4_file = restore(fibdir+(file_search(fibdir, 'crispex*3950*.csav'))[(ncut+1)*i+4])

    # cuts properties
    w1 = cut1_file.spect_pos # the wavelength that it was chosen upon (K2V)
    cut1 = cut1_file.loop_slab[w1,:,:]
    w2 = cut2_file.spect_pos
    cut2 = cut2_file.loop_slab[w2,:,:]
    w3 = cut3_file.spect_pos
    cut3 = cut3_file.loop_slab[w3,:,:]
    w4 = cut4_file.spect_pos
    cut4 = cut4_file.loop_slab[w2,:,:]
    fb_x, fb_y = fib_file.x_loop_pts, fib_file.y_loop_pts
    cut1_x, cut1_y = cut1_file.x_loop_pts, cut1_file.y_loop_pts
    cut2_x, cut2_y = cut2_file.x_loop_pts, cut2_file.y_loop_pts
    cut3_x, cut3_y = cut3_file.x_loop_pts, cut3_file.y_loop_pts
    cut4_x, cut4_y = cut4_file.x_loop_pts, cut4_file.y_loop_pts

    plt.close('all')


    f = plt.figure(figsize=(8,4))
    plt.subplots_adjust(left=0.07,
                        bottom=0.,
                        right=0.88,
                        top=1.,
                        wspace=0.2,
                        hspace=0.2
    )

    ax00 = plt.subplot2grid((4,9), (0,0), colspan = 2, rowspan = 2)
    ax01 = plt.subplot2grid((4,9), (0,2), colspan = 2, rowspan = 2)
    ax02 = plt.subplot2grid((4,9), (2,0), colspan = 2, rowspan=2)
    ax03= plt.subplot2grid((4,9), (2,2), colspan = 2, rowspan=2)
    ax1 = plt.subplot2grid((4,9), (0,5), colspan = 1, rowspan = 4)
    ax2 = plt.subplot2grid((4,9), (0,6), colspan = 1, rowspan = 4)
    ax3 = plt.subplot2grid((4,9), (0,7), colspan = 1, rowspan = 4)
    ax4 = plt.subplot2grid((4,9), (0,8), colspan = 1, rowspan = 4)


    # plotting the cuts
    edge = 5 # number edgegin pixels (diagonaly)
    
    t = fib_file.t_saved
    xmin = np.int(np.min([np.min(fb_x),np.min(cut1_x), np.min(cut4_x)]))
    xmax = np.int(np.max([np.max(fb_x),np.max(cut1_x), np.max(cut4_x)]))
    ymin = np.int(np.min([np.min(fb_y),np.min(cut1_y), np.min(cut4_y)]))
    ymax = np.int(np.max([np.max(fb_y),np.max(cut1_y), np.max(cut4_y)]))

    # The intensity panels
    # Ca K line center
    cak_cropped_c = cube_cak[:, 0, :, ymin-edge:ymax+edge, xmin-edge:xmax+edge]
    
    # cropped cube dimensions
    nt = cak_cropped_c.shape[0]
    nw = cak_cropped_c.shape[1]
    xlim = cak_cropped_c.shape[3]
    ylim = cak_cropped_c.shape[2]

    ax00.imshow(cube_cak[t, 0, w1,
                        ymin-edge:ymax+edge,
                        xmin-edge:xmax+edge],
               origin = 'lower', cmap = 'gray',
    )

    # overplot the fibril
    ax00.plot(fb_x - xmin+edge, fb_y - ymin+edge,
             color = 'red',
             alpha = 0.5,
    )

    #overplot the cuts
    linestyle = '--'
    linewidth = 1.
    color = 'yellow'
    alpha = 0.5
    ax00.plot(cut1_x - xmin+edge, cut1_y - ymin+edge, linestyle = linestyle, linewidth = linewidth, alpha = alpha, color = color)
    ax00.plot(cut2_x - xmin+edge, cut2_y - ymin+edge, linestyle = linestyle, linewidth = linewidth, alpha = alpha, color = color)
    ax00.plot(cut3_x - xmin+edge, cut3_y - ymin+edge, linestyle = linestyle, linewidth = linewidth, alpha = alpha, color = color)
    ax00.plot(cut4_x - xmin+edge, cut4_y - ymin+edge, linestyle = linestyle, linewidth = linewidth, alpha = alpha, color = color)
    
    ax00.text(cut1_x[0]-xmin+edge-3, cut1_y[0]-ymin+edge, '1', alpha = alpha, color = color)
    ax00.text(cut2_x[0]-xmin+edge-3, cut2_y[0]-ymin+edge, '2', alpha = alpha, color = color)
    ax00.text(cut3_x[0]-xmin+edge-3, cut3_y[0]-ymin+edge, '3', alpha = alpha, color = color)
    ax00.text(cut4_x[0]-xmin+edge-3, cut4_y[0]-ymin+edge, '4', alpha = alpha, color = color)

    # Axis properties
    ytick_pos = np.arange(0,(np.round((ylim)/40)+1)*40,40)
    ytick_lab = ytick_pos*res
    xtick_pos = np.arange(0,(np.round((xlim)/40)+1)*40,40)
    xtick_lab = xtick_pos*res
    ytick_pos_t = np.arange(0,(np.round((nt)/30)+1)*30,30)
    ytick_lab_t = ytick_pos_t*tres
    # cut1 length
    cut1_l = np.sqrt(abs(cut1_x[0]-cut1_x[-1])**2 + abs(cut1_y[0]-cut1_y[-1])**2)
    cut2_l = np.sqrt(abs(cut2_x[0]-cut2_x[-1])**2 + abs(cut2_y[0]-cut2_y[-1])**2)
    cut3_l = np.sqrt(abs(cut3_x[0]-cut3_x[-1])**2 + abs(cut3_y[0]-cut3_y[-1])**2)
    cut4_l = np.sqrt(abs(cut4_x[0]-cut4_x[-1])**2 + abs(cut4_y[0]-cut4_y[-1])**2)
    tt = 20
    xtick_pos_1 = np.array([0,27])#np.arange(0,(np.round((cut1_l)/tt)+1)*tt,tt)
    xtick_lab_1 = np.round(xtick_pos_1*res)
    xtick_pos_2 = np.array([0,27])#np.arange(0,(np.round((cut1_l)/tt)+1)*tt,tt)
    xtick_lab_2 = np.round(xtick_pos_2*res)
    xtick_pos_3 = np.array([0,27])#np.arange(0,(np.round((cut1_l)/tt)+1)*tt,tt)
    xtick_lab_3 = np.round(xtick_pos_3*res)
    xtick_pos_4 = np.array([0,27])#np.arange(0,(np.round((cut1_l)/tt)+1)*tt,tt)
    xtick_lab_4 = np.round(xtick_pos_4*res)


    ax00.set_xlim(0, xlim)
    ax00.set_ylim(0, ylim)
    ax00.set_xlabel('x [arcsec]', fontdict = font)
    ax00.set_ylabel('y [arcsec]', fontdict = font)
    ax00.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax00.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax00.tick_params(which='minor', length=2)
    ax00.set_yticks(ytick_pos)
    ax00.set_yticklabels(ytick_lab, fontdict=font)
    ax00.set_xticks(xtick_pos)
    ax00.set_xticklabels(xtick_lab, fontdict=font)

    ax01.set_xlim(0, xlim)
    ax01.set_ylim(0, ylim)
    ax01.set_xlabel('x [arcsec]', fontdict = font)
    ax01.set_ylabel('y [arcsec]', fontdict = font)
    ax01.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax01.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax01.tick_params(which='minor', length=2)
    ax01.set_yticks(ytick_pos)
    ax01.set_yticklabels(ytick_lab, fontdict=font)
    ax01.set_xticks(xtick_pos)
    ax01.set_xticklabels(xtick_lab, fontdict=font)

    ax02.set_xlim(0, xlim)
    ax02.set_ylim(0, ylim)
    ax02.set_xlabel('x [arcsec]', fontdict = font)
    ax02.set_ylabel('y [arcsec]', fontdict = font)
    ax02.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax02.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax02.tick_params(which='minor', length=2)
    ax02.set_yticks(ytick_pos)
    ax02.set_yticklabels(ytick_lab, fontdict=font)
    ax02.set_xticks(xtick_pos)
    ax02.set_xticklabels(xtick_lab, fontdict=font)

    ax03.set_xlim(0, xlim)
    ax03.set_ylim(0, ylim)
    ax03.set_xlabel('x [arcsec]', fontdict = font)
    ax03.set_ylabel('y [arcsec]', fontdict = font)
    ax03.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax03.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax03.tick_params(which='minor', length=2)
    ax03.set_yticks(ytick_pos)
    ax03.set_yticklabels(ytick_lab, fontdict=font)
    ax03.set_xticks(xtick_pos)
    ax03.set_xticklabels(xtick_lab, fontdict=font)

    # Ca K2V
    w_k2v  = 10
    ax01.imshow(cube_cak[t, 0, w_k2v,
                         ymin-edge:ymax+edge,
                         xmin-edge:xmax+edge],
                origin = 'lower', cmap = 'gray',
    )

    # Ca K2R
    w_k2r = 16
    ax02.imshow(cube_cak[t, 0, w_k2r,
                         ymin-edge:ymax+edge,
                         xmin-edge:xmax+edge],
                origin = 'lower', cmap = 'gray',
    )

    # Ca K2R
    w_k2r = 16
    ax03.imshow(cube_cak[t, 0, w_k2r,
                         ymin-edge:ymax+edge,
                         xmin-edge:xmax+edge],
                origin = 'lower', cmap = 'gray',
    )


    ax1.imshow(cut1, origin='lower',cmap = 'gray')
    ax2.imshow(cut2, origin='lower',cmap = 'gray')
    ax3.imshow(cut3, origin='lower',cmap = 'gray')
    ax4.imshow(cut4, origin='lower',cmap = 'gray')
    
    
    ax1.set_xlim(0, cut1_l)
    ax1.set_ylim(0, nt)
    ax1.set_xlabel('x [arcsec]', fontdict = font)
    ax1.set_ylabel('y [minutes]', fontdict = font)
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.tick_params(which='minor', length=2)
    ax1.set_yticks(ytick_pos_t)
    ax1.set_yticklabels(ytick_lab_t, fontdict=font)
    ax1.set_xticks(xtick_pos_1)
    ax1.set_xticklabels(xtick_lab_1, fontdict=font)

    ax2.set_xlim(0, cut2_l)
    ax2.set_ylim(0, nt)
    #ax2.set_xlabel(xlabel, fontdict = font)
    #ax2.set_ylabel(ylabel, fontdict = font)
    ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.tick_params(which='minor', length=2)
    ax2.set_yticks(ytick_pos_t)
    ax2.set_yticklabels(ytick_lab_t, fontdict=font)
    ax2.set_xticks(xtick_pos_2)
    ax2.set_xticklabels(xtick_lab_2, fontdict=font)

    
    ax3.set_xlim(0, cut3_l)
    ax3.set_ylim(0, nt)
    #ax3.set_xlabel(xlabel, fontdict = font)
    #ax3.set_ylabel(ylabel, fontdict = font)
    ax3.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax3.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax3.tick_params(which='minor', length=2)
    ax3.set_yticks(ytick_pos_t)
    ax3.set_yticklabels(ytick_lab_t, fontdict=font)
    ax3.set_xticks(xtick_pos_3)
    ax3.set_xticklabels(xtick_lab_3, fontdict=font)

    ax4.set_xlim(0, cut4_l)
    ax4.set_ylim(0, nt)
    #ax4.set_xlabel(xlabel, fontdict = font)
    #ax4.set_ylabel(ylabel, fontdict = font)
    ax4.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax4.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax4.tick_params(which='minor', length=2)
    ax4.set_yticks(ytick_pos_t)
    ax4.set_yticklabels(ytick_lab_t, fontdict=font)
    ax4.set_xticks(xtick_pos_4)
    ax4.set_xticklabels(xtick_lab_4, fontdict=font)



    #plt.tight_layout()
    plt.show()
