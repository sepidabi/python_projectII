# A brief script to generate cool images
# for the thesis
# from chromospheric observation at SST
# written by Sepideh Kianfar

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


datadir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
outdir = '/home/seki2695/OUTPUT/project2/video/'
w_cak = 11
w_fe = 7
w_ca8 = 9
fr = 59

# data extract
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
cak_cont = cube_cak[:,0,-1,:,:]
fe_cont = cube_fe[:,0,-1,:,:]

# RoI coords
x1, x2, y1, y2 = 425, 1090, 10, 1110

# photosphere LPtot
#LPtot = np.zeros((cube_fe.shape[2], cube_fe.shape[3]))
#for ii in range(cube_fe.shape[1]-1):
#    LPtot += np.sqrt(cube_fe[1,ii,:,:]**2 + cube_fe[2,ii,:,:]**2)/fe_cont
#LPtot = LPtot/cube_fe.shape[1]

# photosphere CPtot
#CPtot = unsharp(0.5*(np.mean(cube_fe[3,0:7,:,:], axis=0) - np.mean(cube_fe[3,7:-1,:,:], axis=0))/fe_cont, alpha = 0.5, sigma = 2.)

# CaK wavelength-summed and unsharped
wi, wf = 10, 16
cak_sum_un = unsharp(np.mean(cube_cak[:, 0, wi:wf+1, :, :],axis=0), alpha = 0.5, sigma = 10.)

# CaK Core
cak_wc = 13
cak_core = unsharp(cube_cak[:, 0, cak_wc, :, :], alpha = 0.5, sigma = 5.)

# CaK k2v
wv = 11
cak_k2v = unsharp(cube_cak[:, 0, wv, :, :], alpha = 0.5, sigma = 5.)

# CaK k2r
wr = 15
cak_k2r = unsharp(cube_cak[:, 0, wr, :, :], alpha = 0.5, sigma = 5.)


# Ca 8542 Core
ca8_wc = 11
ca8_core = unsharp(cube_ca8[:, 0,ca8_wc,:,:], alpha = 0.5, sigma = 5.)

text = ['', r'a) continuum 4000 $\mathrm{\AA}$', r'b) Ca II K', r'continuum 4000 $\mathrm{\AA}$', r'Ca II 8542 $\mathrm{\AA}$ core', r'Ca II K2V', r'Ca II K2R']
textcolor = ['','white', 'white', 'white', 'white', 'white', 'white']

xx, yy = x2-x1, y2-y1

# axis info
sc_fact = 266.7 # axis scaling factor
ytick_pos = np.arange(0,(np.round(yy/sc_fact)+1)*sc_fact,sc_fact)
ytick_lab = np.round(ytick_pos*res).astype(int)
xtick_pos = np.arange(0,(np.round(xx/sc_fact)+1)*sc_fact,sc_fact)
xtick_lab = np.round(xtick_pos*res).astype(int)

for fr in range(cak_cont.shape[0]):

    # maps
    plt.close('all')

    f = plt.figure(figsize = [7.5,6.1/2.+1])
    gs = gridspec.GridSpec(1,4) # grid scale of the 1st col.
    gs.update(left=0.06,
              right=0.98,
              wspace=0.08,
              bottom=0.07,
              top=1.,
              hspace = 0.01,
    )

    ax3 = plt.subplot(gs[0,0],adjustable = 'box')
    sigma_factor = 3.5
    ax_min = np.mean(cak_cont[fr,y1:y2,x1:x2]) - sigma_factor*np.std(cak_cont[fr,y1:y2,x1:x2])
    ax_max = np.mean(cak_cont[fr,y1:y2,x1:x2]) + sigma_factor*np.std(cak_cont[fr,y1:y2,x1:x2])
    ax3.imshow(cak_cont[fr,y1:y2,x1:x2], origin = 'lower', cmap = 'gray', vmin = ax_min, vmax = ax_max)
    ax3.set_xticks(xtick_pos)
    ax3.set_xlim(0,xx)
    ax3.set_ylim(0,yy)
    ax3.set_yticks(ytick_pos)
    ax3.set_ylabel(r'y [arcsec]', fontdict=font)
    ax3.set_yticklabels(ytick_lab, fontdict = font)
    ax3.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax3.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax3.tick_params(which='minor', length=2)
    ax3.set_xticks(xtick_pos)
    ax3.set_xticklabels(xtick_lab, fontdict = font)
    ax3.set_xlabel(r'x [arcsec]', fontdict=font)
    ax3.text(10, 10 , text[3], fontdict = font, color = textcolor[3])
    
    ax4 = plt.subplot(gs[0,1],adjustable = 'box')
    sigma_factor = 3.5
    ax_min = np.mean(ca8_core[fr,y1:y2,x1:x2]) - sigma_factor*np.std(ca8_core[fr,y1:y2,x1:x2])
    ax_max = np.mean(ca8_core[fr,y1:y2,x1:x2]) + sigma_factor*np.std(ca8_core[fr,y1:y2,x1:x2])
    ax4.imshow(ca8_core[fr,y1:y2,x1:x2], origin = 'lower', cmap = 'gray', vmin = ax_min, vmax = ax_max)
    ax4.set_xticks(xtick_pos)
    ax4.set_xlim(0,xx)
    ax4.set_ylim(0,yy)
    ax4.set_yticks(ytick_pos)
    ax4.set_ylabel('')
    ax4.set_yticklabels([])
    ax4.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax4.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax4.tick_params(which='minor', length=2)
    ax4.set_xticks(xtick_pos)
    ax4.set_xticklabels(xtick_lab, fontdict = font)
    ax4.set_xlabel(r'x [arcsec]', fontdict=font)
    ax4.text(10, 10 , text[4], fontdict = font, color = textcolor[4])
    
    ax5 = plt.subplot(gs[0,2],adjustable = 'box')
    sigma_factor = 3.5
    ax_min = np.mean(cak_k2v[fr,y1:y2,x1:x2]) - sigma_factor*np.std(cak_k2v[fr,y1:y2,x1:x2])
    ax_max = np.mean(cak_k2v[fr,y1:y2,x1:x2]) + sigma_factor*np.std(cak_k2v[fr,y1:y2,x1:x2])
    ax5.imshow(cak_k2v[fr,y1:y2,x1:x2], origin = 'lower', cmap = 'gray', vmin = ax_min, vmax = ax_max)
    ax5.set_xticks(xtick_pos)
    ax5.set_xlim(0,xx)
    ax5.set_ylim(0,yy)
    ax5.set_yticks(ytick_pos)
    ax5.set_ylabel('')
    ax5.set_yticklabels([])
    ax5.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax5.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax5.tick_params(which='minor', length=2)
    ax5.set_xticks(xtick_pos)
    ax5.set_xticklabels(xtick_lab, fontdict = font)
    ax5.set_xlabel(r'x [arcsec]', fontdict=font)
    ax5.text(10, 10 , text[5], fontdict = font, color = textcolor[5])
    
    ax6 = plt.subplot(gs[0,3],adjustable = 'box')
    sigma_factor = 3.5
    ax_min = np.mean(cak_k2r[fr,y1:y2,x1:x2]) - sigma_factor*np.std(cak_k2r[fr,y1:y2,x1:x2])
    ax_max = np.mean(cak_k2r[fr,y1:y2,x1:x2]) + sigma_factor*np.std(cak_k2r[fr,y1:y2,x1:x2])
    ax6.imshow(cak_k2r[fr,y1:y2,x1:x2], origin = 'lower', cmap = 'gray', vmin = ax_min, vmax = ax_max)
    ax6.set_xticks(xtick_pos)
    ax6.set_xlim(0,xx)
    ax6.set_ylim(0,yy)
    ax6.set_yticks(ytick_pos)
    ax6.set_ylabel('')
    ax6.set_yticklabels([])
    ax6.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax6.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax6.tick_params(which='minor', length=2)
    ax6.set_xticks(xtick_pos)
    ax6.set_xticklabels(xtick_lab, fontdict = font)
    ax6.set_xlabel(r'x [arcsec]', fontdict=font)
    ax6.text(10, 10 , text[6], fontdict = font, color = textcolor[6])

    # clipping the image
    #sigma_factor = 4.
    #cak_min = np.mean(cak_scan) - 2.5*np.std(cak_scan)
    #cak_max = np.mean(cak_scan) + 6*np.std(cak_scan)
    
    #plt.tight_layout()
    #plt.show()
    ax3.set_title('frame ' + str(fr).zfill(3), fontdict = font, loc = 'left')
    #break
    filename = 'maps'
    print fr
    plt.savefig(outdir + filename + str(fr).zfill(4) + '.png', dpi = 1000, quality = 100)

# creating the video
os.system("ffmpeg -r 2 -i "+outdir+filename+"%04d.png -vcodec mpeg4 -y "+outdir+filename+".mp4")
