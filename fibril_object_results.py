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
from matplotlib.lines import Line2D
from matplotlib.markers import MarkerStyle
from matplotlib.transforms import Affine2D

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

########
# functions
#########
def gen_arrow_head_marker(rot):
    """generate a marker to plot with matplotlib scatter, plot, ...

    https://matplotlib.org/stable/api/markers_api.html#module-matplotlib.markers

    rot=0: positive x direction
    Parameters
    ----------
    rot : float
        rotation in degree
        0 is positive x direction

    Returns
    -------
    arrow_head_marker : Path
        use this path for marker argument of plt.scatter
    scale : float
        multiply a argument of plt.scatter with this factor got get markers
        with the same size independent of their rotation.
        Paths are autoscaled to a box of size -1 <= x, y <= 1 by plt.scatter
    """
    arr = np.array([[0, 0], [-.9,.2], [-.9, -.2], [0, 0]])  # arrow shape
    angle = rot / 180 * np.pi
    rot_mat = np.array([
        [np.cos(angle), np.sin(angle)],
        [-np.sin(angle), np.cos(angle)]
        ])
    arr = np.matmul(arr, rot_mat)  # rotates the arrow

    # scale
    x0 = np.amin(arr[:, 0])
    x1 = np.amax(arr[:, 0])
    y0 = np.amin(arr[:, 1])
    y1 = np.amax(arr[:, 1])
    scale = np.amax(np.abs([x0, x1, y0, y1]))
    codes = [mpl.path.Path.MOVETO, mpl.path.Path.LINETO,mpl.path.Path.LINETO, mpl.path.Path.CLOSEPOLY]
    arrow_head_marker = mpl.path.Path(arr, codes)
    return arrow_head_marker, scale


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
objresdir = '/home/seki2695/OUTPUT/project2/result_object/'

# Difinitions
res = 0.0375 # CHROMIS pixel size in arcsec
cad = 8
#fr_select = [29,75,135,167] # frames to display
w_pos = 13 # wavelength position to display (K2V, K3, K2R) = (10,13,16)

plt.close("all")

class oscillation_result:
    def __init__(self,
                 per_los,
                 per_los_m,
                 per_los_dom,
                 per_pos,
                 per_pos_m,
                 per_pos_dom,
                 A_los_m,
                 A_pos_m,
                 dphi,
                 ncc_min,
                 ncc_max,
                 mid_coord,
                 cor_rate,
                 theta,
                 per_m,
                 per_dom,
                 ratio,
                 dist,
                 ):
        
        self.per_los = per_los
        self.per_los_m = per_los_m
        self.per_los_dom = per_los_dom
        self.per_pos = per_pos
        self.per_pos_m = per_pos_m
        self.per_pos_dom = per_pos_dom
        self.A_los_m = A_los_m
        self.A_pos_m = A_pos_m
        self.dphi = dphi
        self.ncc_min = ncc_min
        self.ncc_max =  ncc_max
        self.mid_coord = mid_coord
        self.cor_rate = cor_rate
        self.theta = theta
        self.per_m = per_m
        self.per_dom = per_dom
        self.ratio = ratio
        self.dist = dist


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
#cak_int = sharpen(np.mean(cube_cak[fr,0, cak_wc-3:cak_wc+3, :, :],axis = 0), sigma =[3,3], unsigma = [1,1])
cak_int = unsharp(np.mean(cube_cak[fr,0, cak_wc-3:cak_wc+3, :, :],axis = 0), alpha = 0.5, sigma = 2.)
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
A_los_m = np.zeros(0)
A_pos_m = np.zeros(0)
dphi = np.zeros(0)
ncc_min = np.zeros(0)
ncc_max =  np.zeros(0)
mid_coord_x = np.zeros(0)
mid_coord_y = np.zeros(0)
cor_rate = np.zeros(0)
theta = np.zeros(0)
per_m = np.zeros(0)
per_dom = np.zeros(0)
ratio = np.zeros(0)
dist = np.zeros(0)

# Extract the oscillation objects

# User contribution to the specify the oscillation
ind = 0

for i in range(len(file_search(objresdir, "*.obj"))):
    res_obj = file_search(objresdir, "*.obj")[i]
    result = load_obj(objresdir + res_obj)
    
    per_los = np.append(per_los, result.per_los)
    per_los_m = np.append(per_los_m, result.per_los_m)
    per_los_dom = np.append(per_los_dom, result.per_los_dom)
    per_pos = np.append(per_pos, result.per_pos)
    per_pos_m = np.append(per_pos_m, result.per_pos_m)
    per_pos_dom = np.append(per_pos_dom, result.per_pos_dom)
    A_los_m = np.append(A_los_m, result.A_los_m)
    A_pos_m = np.append(A_pos_m, result.A_pos_m)
    dphi = np.append(dphi, result.dphi)
    ncc_min = np.append(ncc_min, result.ncc_min)
    ncc_max =  np.append(ncc_max, result.ncc_max)
    mid_coord_x = np.append(mid_coord_x, result.mid_coord[0])
    mid_coord_y = np.append(mid_coord_y, result.mid_coord[1])
    cor_rate = np.append(cor_rate, result.cor_rate)
    theta = np.append(theta, result.theta)
    per_m = np.append(per_m, result.per_m)
    per_dom = np.append(per_dom, result.per_dom)

# Print statistics for the paper
print('     , min , max , mean , dev')
print('P_pos_m: ', np.round(np.array([np.min(per_pos_m), np.max(per_pos_m), np.mean(per_pos_m), np.std(per_pos_m)]), decimals = 0))
print('P_los_m: ', np.round(np.array([np.min(per_los_m), np.max(per_los_m), np.mean(per_los_m), np.std(per_los_m)]), decimals = 0))
print('A_pos_m: ', np.round(np.array([np.min(A_pos_m), np.max(A_pos_m), np.mean(A_pos_m), np.std(A_pos_m)]), decimals = 2))
print('A_los_m: ', np.round(np.array([np.min(A_los_m), np.max(A_los_m), np.mean(A_los_m), np.std(A_los_m)]), decimals = 2))
print("press 'c' to continue to the plots.")
stop()

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
#plt.show()
#
filename = outdir + 'analys_tot.pdf'
fig.savefig(filename, dpi = 1000)
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
n, bins, patches  = axfi.hist(dphi[(dphi>=-1)], bins = 15, edgecolor='white', linewidth=1., color = 'orangered', alpha = 0.5, range = [-1,1])
bin_centers = 0.5 * (bins[:-1] + bins[1:])
axfi.tick_params(axis = 'both', labelsize = 8.)
axfi.set_xlabel(r'$\Delta \phi$ [rad]', fontdict = font)
axfi.set_ylabel(r'# of oscillation cases', fontdict = font)
axfi.set_xticks([-1,-.75, -.5, -.25, 0,0.25,0.5,0.75,1])
axfi.set_xticklabels([r'-$\pi$', r'-3$\pi$/4', r'-$\pi$/2', r'-$\pi$/4', r'0', r'$\pi$/4', r'$\pi$/2', r'3$\pi$/4', r'$\pi$'], fontdict = font)
#axfi.set_yticks([0,3,6,9,12,15])
axfi.set_ylim(0,100)
axfi.set_xlim(-1.05,1.05)
axfi.yaxis.set_minor_locator(AutoMinorLocator(2))

axins = inset_axes(axfi, width=1.1, height=.91,
                    bbox_to_anchor=(0.64, 0.57),
                    bbox_transform=axfi.transAxes, loc=3, borderpad=0)

nn, bbins, ppatches = axins.hist(np.abs(dphi[(dphi>=-1)]), bins = 7, range = [0,1], edgecolor = 'white', linewidth = 1., color = 'orangered', alpha = 0.5)
axins.tick_params(axis = 'both', labelsize = 8.)
axins.set_xlabel(r'|$\Delta \phi$| [rad]', fontdict = font)
#axins.set_ylabel(r'Number of oscillation cases', fontdict = font)
#axins.set_yticks([0,5,10,15])
axins.set_xticks([0,0.5,1])
axins.set_xlim(-.05,1.05)
axins.set_xticklabels([r'0', r'$\pi$/2', r'$\pi$'], fontdict = font)
#axins.yaxis.set_minor_locator(AutoMinorLocator(2))

plt.tight_layout()
filename = outdir + 'analys_dphi.pdf'
#plt.show()
fifi.savefig(filename, dpi = 1000)
#stop()
#
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
ssize = 900.
# RoI coords
x1, x2, y1, y2 = 425, 900, 400, 1110
xx, yy = x2-x1, y2-y1
plt.close("all")
map_cak = plt.figure(figsize = (4,6))
gs_map = gridspec.GridSpec(ncols = 1,nrows = 1, figure = map_cak)
ax_map =  plt.subplot(gs_map[0,0])
roi = cak_int[y1:y2, x1:x2]
ax_map.imshow(roi, cmap = 'gray', origin = 'lower')
c_levels = [4.5,8] # contour levels
ax_map.contourf(np.abs(CPtot[y1:y2, x1:x2]).squeeze()*1e2, cmap = 'Blues', alpha = 0.25, levels = c_levels)
contour = ax_map.contour(np.abs(CPtot[y1:y2, x1:x2]).squeeze()*1e2, colors = ['white', 'white'], alpha = 0.5, levels = c_levels, linewidths = [0.5,0.5])
#common_style = {k: v for k, v in filled_marker_style.items() if k != 'marker'}
for i in range (len(mid_coord_x)):
    marker, scale = gen_arrow_head_marker(90+theta[i])
    #marker = (3, 0, theta[i])
    markersize = np.power(per_m/np.max(per_m), 2)*ssize
    if (dphi[i]>=-1):
        sc = ax_map.scatter(x = mid_coord_x[i]-x1, y = mid_coord_y[i]-y1, c = np.abs(cor_rate[i]), marker=marker, cmap = cmap_phi, s = markersize[i], alpha = 0.45, vmin =0, vmax = 1)
        dist_i = np.sqrt(contour.find_nearest_contour(mid_coord_x[i]-x1, mid_coord_y[i]-y1)[5])
        dist = np.append(dist, dist_i)
        ratio = np.append(ratio,per_m[i])#np.append(ratio,dist_i/per_m[i])
        #print(np.power(per_m[i]/np.max(per_m), 2)*ssize)

    #elif (dphi[i]==-5):
        #ax_map.scatter(x = mid_coord_x[i]-x1, y = mid_coord_y[i]-y1, marker=marker , s = markersize, alpha = 0.6, color = 'gray')
        #dist_i = np.sqrt(contour.find_nearest_contour(mid_coord_x[i]-x1, mid_coord_y[i]-y1)[5])
        #dist = np.append(dist, dist_i)
        #ratio = np.append(ratio,dist_i/per_m[i])
    #else:
        #ax_map.scatter(x = mid_coord_x[i]-x1, y = mid_coord_y[i]-y1, marker='.' , s = 10, alpha = 0.6, color = 'white')
        #print('-',mid_coord_x[i],mid_coord_y[i],dphi[i])
#plt.show()
#
#ax_map.scatter(x=100,  y = 100, marker = (3,0,0), s = np.power(140/np.max(per_m), 2)*ssize)
hands =np.array([120,230,340,450])
sizes = np.sqrt(np.power(hands/np.max(per_m), 2)*ssize)
marker, scale = gen_arrow_head_marker(0)
legend_elements = [Line2D([0], [0], marker=marker, color='w', label=str(hands[0]),markerfacecolor='gray', markersize=sizes[0], linestyle='None', markeredgecolor = 'darkgray'), #linewidth = 1.),
                   Line2D([], [], marker=marker, color='w', label=str(hands[1]),markerfacecolor='gray', markersize=sizes[1], linestyle='None', markeredgecolor = 'darkgray'),
                   Line2D([], [], marker=marker, color='w', label=str(hands[2]),markerfacecolor='gray', markersize=sizes[2], linestyle='None', markeredgecolor = 'darkgray'),
                   Line2D([], [], marker=marker, color='w', label=str(hands[3]),markerfacecolor='gray', markersize=sizes[3], linestyle='None', markeredgecolor = 'darkgray'),
                   Line2D([], [], marker=marker, color='w', label='',markerfacecolor='gray', markersize=0, linestyle='None'),
                   Line2D([], [], marker=marker, color='w', label='',markerfacecolor='gray', markersize=0, linestyle='None'),
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
#plt.show()
#stop()
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
map_cak.savefig(filename, dpi = 1000)
print('file saved to: '+ filename)
#plt.show()
#stop()
plt.close(map_cak)


# fibril didstance and period correlation
#=======================
fig = plt.figure(figsize = (4,4))
y = ratio
x = dist*res


plt.scatter(x,y,alpha = 0.3, color = 'orangered')

pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
plt.annotate('pearsonr = ' + pearsonr,xy=(0.6,0.92), fontsize = 8., xycoords='axes fraction')

plt.ylabel(r'$\overline{P}$ [s]', fontdict = font)
plt.xlabel(r'$|d|$ [arcsec]', fontdict = font)
plt.xticks(fontsize = 8.)
plt.yticks(fontsize = 8.)

plt.tight_layout()
filename = outdir + 'per_d.pdf'
fig.savefig(filename, dpi = 1000)
print('file saved to: '+ filename)
plt.show()
stop()



##########
# JOINT PLOTS
##########

# Period PLOT
# ========
fontsize = 8.
xy = np.linspace(-0,1000,100)

plt.close('all')
fig = plt.figure(figsize = (4,4))
bins = 30
gs = gridspec.GridSpec(ncols = 5,nrows = 5)

ax_joint = fig.add_subplot(gs[1:5,0:4])
ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_y = fig.add_subplot(gs[1:5,4])

# Period
x = per_pos_m[np.where(dphi>=-1)]
y = per_los_m[np.where(dphi>=-1)]

rmin = 60#np.min([x,y])
rmax = 500#np.max([x,y])

angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
ax_joint.text(xmax - 60,ymax - 55,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

cmapc = mpl.colors.LinearSegmentedColormap.from_list("", ["white","orangered"])
#ax_joint.scatter(x, y, alpha = 0.25, color = 'orangered')
#ax_joint.hist2d(x, y, bins = bins, cmap = cmapc)
sns.kdeplot(x,y, cmap = 'Reds', ax = ax_joint, shade = True)
sns.kdeplot(y, ax = ax_marg_y, shade=True, color = 'orangered', vertical = True)
sns.kdeplot(x, ax = ax_marg_x, shade=True, color = 'orangered')
#ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = bins)
#ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = bins)

# Set ax limits on marginals
ax_joint.set_xlim([rmin,rmax])
ax_joint.set_ylim([rmin,rmax])
ax_marg_y.set_ylim([rmin,rmax])
ax_marg_x.set_xlim([rmin,rmax])

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
#plt.show()
#stop()

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

# Amplitude
x = A_pos_m[np.where(dphi>=-1)]
y = A_los_m[np.where(dphi>=-1)]

rmin = 0#np.min([x,y])
rmax = 3.#np.max([x,y])

angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
ax_joint.text(xmax - 0.6,ymax - 0.55,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

#ax_joint.scatter(x, y, alpha = 0.25, color = 'orangered')
#ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = bins)
#ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = bins/3)
sns.kdeplot(x,y, cmap = 'Reds', ax = ax_joint, shade = True)
sns.kdeplot(y, ax = ax_marg_y, shade=True, color = 'orangered', vertical = True)
sns.kdeplot(x, ax = ax_marg_x, shade=True, color = 'orangered')

# Set ax limits on marginals
ax_joint.set_xlim([rmin,rmax])
ax_joint.set_ylim([rmin,rmax])
ax_marg_y.set_ylim([rmin,rmax])
ax_marg_x.set_xlim([rmin,rmax])

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
#plt.show()
#stop()

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

# POS
x = per_pos_m[np.where(dphi>=-5)]
y = A_pos_m[np.where(dphi>=-5)]

xmin, ymin = 75, 0.13
xmax, ymax = 500,2.5


#angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
#ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
#ax_joint.text(xmax - 60,ymax - 0.55,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

#ax_joint.scatter(x, y, alpha = 0.25, color = 'orangered')
#ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = bins)
#ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = bins)
sns.kdeplot(x,y, cmap = 'Reds', ax = ax_joint, shade = True)
sns.kdeplot(y, ax = ax_marg_y, shade=True, color = 'orangered', vertical = True)
sns.kdeplot(x, ax = ax_marg_x, shade=True, color = 'orangered')


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
#plt.show()
#stop()

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

# LOS
x = per_los_m[np.where(dphi>=-5)]
y = A_los_m[np.where(dphi>=-5)]

xmin, ymin = 70, 0
xmax, ymax = 500,3.

#angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
#ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
#ax_joint.text(xmax - 60,ymax - 0.55,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

#ax_joint.scatter(x, y, alpha = 0.25, color = 'orangered')
#ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = bins)
#ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = bins/2)
sns.kdeplot(x,y, cmap = 'Reds', ax = ax_joint, shade = True)
sns.kdeplot(y, ax = ax_marg_y, shade=True, color = 'orangered', vertical = True)
sns.kdeplot(x, ax = ax_marg_x, shade=True, color = 'orangered')


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
plt.close('all')        
