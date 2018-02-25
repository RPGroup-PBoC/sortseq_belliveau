import os
import glob
import pickle
import re

# Our numerical workhorses
import numpy as np
import pandas as pd

# Import the project utils
import sys
sys.path.insert(0, '../')
import NB_sortseq_utils as utils

# Import matplotlib stuff for plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Logo-generating module
import anylogo

utils.set_plotting_style_emat()
#===============================================================================
# Set output directory based on the graphicspath.tex file to print in dropbox
#===============================================================================
output = 'output_figs/'
#------------------------------------------------------------------------------#
# directory where emat csv files are contained
#------------------------------------------------------------------------------#
datadir = '../sortseq/20160707_purT_xylE_dgoR/'

#------------------------------------------------------------------------------#
# putative xylR
#------------------------------------------------------------------------------#
# Set color scale - I want the colorbar to be symmetric
emat_min=-0.4
emat_max=0.4
mid_val=0.0

# Create background array
gc = .508
background_array =pd.DataFrame( [[(1-gc)/2,gc/2,gc/2,(1-gc)/2]])

# energy_df = pd.read_csv(datadir1 + 'pymc_xylose_mut1_mut2_004_emat_mean_mod.txt')
energy_df = pd.read_csv(datadir + '20160710_xylE_MG1655_M9xylose_na_mut2_4bins_XylR_emat_mean.csv')
energy_df = energy_df[['A','C','G','T']]

energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
energy_df_scaled = energy_df_scaled[['A','C','G','T']]

# create background nucleotide frequencies dataframe
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])

seq = 'ACAGAAAAGACATTACGTAAACGCATTGTAAAAAATGATAA'

plt.figure(figsize=utils.cm2inch((0.18*len(seq) + 0.2,2.5)))

ax = plt.gca()
relative_scale=1.5
relative_spacing=.65
emat_ymin = -2 * (relative_scale + relative_spacing)
emat_ymax = -2 * relative_spacing
yticks = np.linspace(emat_ymin, emat_ymax, 9)[[1, 3, 5, 7]]
yticklabels = list('TGCA')
anylogo.draw(ax, effect_df=energy_df_scaled, logo_type='information',
             background = background_df,
             use_transparency=False)
L = len(seq)
ax.set_xticks([])

im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), seq),
            interpolation='none',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max),
            extent=(-.5, L - .5, emat_ymin, emat_ymax),
            zorder=100,
            aspect='auto')

ax.set_ylim([emat_ymin, 2])
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels, fontsize=5, horizontalalignment='center')
ax.set_ylabel('')
ax.yaxis.set_tick_params(length=0)

# # create an axes on the right side of ax. The width of cax will be 3%
# # of ax and the padding between cax and ax will be fixed at 0.05 inch.
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="3%", pad=0.05)
#
# cbar = plt.colorbar(im, cax=cax, ticks=[-0.5, 0, 0.5])
# cbar.ax.set_yticklabels(['-0.5', '0', '0.5'], fontsize=6, fontname='Arial')
# cbar.outline.set_visible(False)
# cbar.ax.tick_params(axis=u'both', which=u'both',length=0)

y = .5*emat_ymax
for i in range(L):
    ax.text(i, y, seq[i], horizontalalignment='center', verticalalignment='center',
    fontsize=6)
    ax.tick_params(axis='y', pad=7)

plt.tight_layout()
plt.savefig(output + 'fig6_xylE_logo_matrix_right.pdf')




# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(seq)
im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), seq),
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max),
            zorder=100,
            aspect='auto')
ax.axis('off')

plt.savefig(output + 'fig6_xylE_matrix_right_matrixonly.pdf')



# #------------------------------------------------------------------------------#
# # putative CRP
# #------------------------------------------------------------------------------#
#
# # energy_df = pd.read_csv(datadir1 + 'pymc_xylose_mut1_mut2_004_emat_mean_mod.txt')
energy_df = pd.read_csv(datadir + '20160710_xylE_MG1655_M9xylose_na_mut3_4bins_CRP_emat_mean.csv')
energy_df = energy_df[['A','C','G','T']]


energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
energy_df_scaled = energy_df_scaled[['A','C','G','T']]

# create background nucleotide frequencies dataframe
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])

seq = 'ATTATCACAATTAAGATCACAGA'

plt.figure(figsize=utils.cm2inch((0.18*len(seq) + 0.2,2.5)))

ax = plt.gca()
relative_scale=1.5
relative_spacing=.65
emat_ymin = -2 * (relative_scale + relative_spacing)
emat_ymax = -2 * relative_spacing
yticks = np.linspace(emat_ymin, emat_ymax, 9)[[1, 3, 5, 7]]
yticklabels = list('TGCA')
anylogo.draw(ax, effect_df=energy_df_scaled, logo_type='information',
             background = background_df,
             use_transparency=False)
L = len(seq)
ax.set_xticks([])

im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), seq),
            interpolation='none',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max),
            extent=(-.5, L - .5, emat_ymin, emat_ymax),
            zorder=100,
            aspect='auto')

ax.set_ylim([emat_ymin, 2])
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels, fontsize=5, horizontalalignment='center')
ax.set_ylabel('')
ax.yaxis.set_tick_params(length=0)

# # create an axes on the right side of ax. The width of cax will be 3%
# # of ax and the padding between cax and ax will be fixed at 0.05 inch.
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="3%", pad=0.05)
#
# cbar = plt.colorbar(im, cax=cax, ticks=[-0.5, 0, 0.5])
# cbar.ax.set_yticklabels(['-0.5', '0', '0.5'], fontsize=6, fontname='Arial')
# cbar.outline.set_visible(False)
# cbar.ax.tick_params(axis=u'both', which=u'both',length=0)

y = .5*emat_ymax
for i in range(L):
    ax.text(i, y, seq[i], horizontalalignment='center', verticalalignment='center',
    fontsize=6)
    ax.tick_params(axis='y', pad=7)

plt.tight_layout()

plt.savefig(output + 'fig6_xylE_logo_matrix_left.pdf')


# save energy matrix using nearest interpolation

# emat_left = utils.zero_matrix_WT(np.array(energy_df.T))
energy_df_left = pd.read_csv(datadir + '20160710_xylE_MG1655_M9xylose_na_mut3_4bins_CRP_emat_mean.csv')
energy_df_left = energy_df_left[['A','C','G','T']]
emat_left = np.array(energy_df_left.T)
print(len(emat_left[0,2:]))
energy_df_right = pd.read_csv(datadir + '20160710_xylE_MG1655_M9xylose_na_mut2_4bins_XylR_emat_mean.csv')
energy_df_right = energy_df_right[['A','C','G','T']]
emat_right = np.array(energy_df_right.T)
print(len(emat_right[0,5:]))
emat_combined = np.zeros([4,57])
emat_combined[:,:21] = utils.fix_matrix_gauge(emat_left[:,2:])
emat_combined[:,21:] = utils.fix_matrix_gauge(emat_right[:,5:])

seq = 'TATCACAATTAAGATCACAGAAAAGACATTACGTAAACGCATTGTAAAAAATGATAA'
plt.figure()
ax = plt.gca()
L = len(seq)
im = ax.imshow(utils.zero_matrix_WT(emat_combined, seq),
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'fig6_xylE_logo_matrix_combined.pdf')




# energy_df_left = pd.read_csv(datadir + '20160710_xylE_MG1655_M9xylose_na_mut3_4bins_CRP_emat_mean.csv')
# energy_df_right = pd.read_csv(datadir + '20160710_xylE_MG1655_M9xylose_na_mut2_4bins_XylR_emat_mean.csv')
