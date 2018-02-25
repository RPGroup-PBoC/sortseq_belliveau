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

# Create background array
gc = .508
background_array =pd.DataFrame( [[(1-gc)/2,gc/2,gc/2,(1-gc)/2]])

#------------------------------------------------------------------------------#
# putative PurR binding site
#------------------------------------------------------------------------------#

# energy_df = pd.read_csv(datadir1 + '20160710_purT_MG1655_M9glucose_adenine_mut1_4bins_PurR_emat_mean.csv')
energy_df = pd.read_csv('input_data/20160710_purT_MG1655_M9glucose_adenine_mut1_4bins_pymc_purR_MCMC_104_fixed_thermo.csv')

# energy_df = pd.read_csv('purR_RNAP_mult',delim_whitespace=True)
seq = 'ACGCAAACGTTTTCGT'
energy_df = energy_df[['A','C','G','T']]

# for plotting logo
energy_df_logo = utils.estimate_scalefactor(np.array(energy_df.copy()))*energy_df.copy()

# for plotting emat in KbT units
scalefactor = 11.5558
energy_df_scaled = energy_df.copy()
energy_df_scaled = scalefactor * energy_df_scaled[['A','C','G','T']]
energy_df_scaled.to_csv('output_figs/purR_matrix.csv')
emat_scaled = scalefactor * utils.zero_matrix_WT(np.array(energy_df[['A','C','G','T']].T), seq)



# create background nucleotide frequencies dataframe
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_logo), 1)), columns=['A','C','G','T'])


# Set color scale - I want the colorbar to be symmetric
emat_min=-5#emat_scaled.min()
emat_max=5#emat_scaled.min()
# print(energy_df_scaled.min())
# print(energy_df_scaled.max())
mid_val=0.0

plt.figure(figsize=utils.cm2inch((0.18*len(seq) + 0.2,2.5)))

ax = plt.gca()
relative_scale=1.5
relative_spacing=.65
emat_ymin = -2 * (relative_scale + relative_spacing)
emat_ymax = -2 * relative_spacing
yticks = np.linspace(emat_ymin, emat_ymax, 9)[[1, 3, 5, 7]]
yticklabels = list('TGCA')
anylogo.draw(ax, effect_df=energy_df_logo, logo_type='information',
             background = background_df,
             use_transparency=False)
L = len(seq)
ax.set_xticks([])

im = ax.imshow(emat_scaled,
            interpolation='nearest',
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

# create an axes on the right side of ax. The width of cax will be 3%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.05)

cbar = plt.colorbar(im, cax=cax, ticks=[emat_min, 0, emat_max])
cbar.ax.set_yticklabels(['-'+str(emat_min), '0', str(emat_min)], fontsize=6, fontname='Arial')
cbar.outline.set_visible(False)
cbar.ax.tick_params(axis=u'both', which=u'both',length=0)

y = .5*emat_ymax
for i in range(L):
    ax.text(i, y, seq[i], horizontalalignment='center', verticalalignment='center',
    fontsize=5)
    ax.tick_params(axis='y', pad=7)

plt.tight_layout()
plt.savefig(output + 'fig5_thermomodel_PurR_emat_logo.pdf')


# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(seq)
im = ax.imshow(emat_scaled,
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'fig5_thermomodel_PurR_ematonly.pdf')

#------------------------------------------------------------------------------#
# RNAP binding site
#------------------------------------------------------------------------------#

# energy_df = pd.read_csv(datadir2 + '20160824_purT_MG1655deltapurR_M9glucose_adenine_mut1_4bins_RNAP_emat_mean.csv')
energy_df = pd.read_csv('input_data/purR_RNAP_mult.csv')
# energy_df = pd.read_csv(datadir2 + '20160710_purT_MG1655_M9glucose_na_mut1_4bins_pymc_rnap_MCMC_114_thermo.csv')

energy_df = energy_df[['A','C','G','T']]
seq = 'AAAGACACACGCAAACGTTTTCGTTTATACTG'

# for plotting logo
energy_df_logo = utils.estimate_scalefactor(np.array(energy_df.copy()))*energy_df.copy()

# for plotting emat in KbT units (note matrix.csv was provided by Bill in kBT units)
energy_df_scaled = energy_df.copy()
energy_df_scaled = energy_df_scaled[['A','C','G','T']]
# emat_scaled = scalefactor * utils.zero_matrix_WT(np.array(energy_df.T), seq)

# create background nucleotide frequencies dataframe
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])

# Set color scale - I want the colorbar to be symmetric
emat_min = -5
emat_max = 5
mid_val=0

plt.figure(figsize=utils.cm2inch((0.18*len(seq) + 0.2,2.5)))

ax = plt.gca()
relative_scale=1.5
relative_spacing=.65
emat_ymin = -2 * (relative_scale + relative_spacing)
emat_ymax = -2 * relative_spacing
yticks = np.linspace(emat_ymin, emat_ymax, 9)[[1, 3, 5, 7]]
yticklabels = list('TGCA')
anylogo.draw(ax, effect_df=energy_df_logo, logo_type='information',
             background = background_df,
             use_transparency=False)
L = len(seq)
ax.set_xticks([])

im = ax.imshow(energy_df_scaled.T,
            interpolation='nearest',
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

# create an axes on the right side of ax. The width of cax will be 3%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.05)

cbar = plt.colorbar(im, cax=cax, ticks=[emat_min, 0, -emat_min])
cbar.ax.set_yticklabels([str(emat_min), '0', str(emat_max)], fontsize=6, fontname=None)
cbar.outline.set_visible(False)
cbar.ax.tick_params(axis=u'both', which=u'both',length=0)

y = .5*emat_ymax
for i in range(L):
    ax.text(i, y, seq[i], horizontalalignment='center', verticalalignment='center',
    fontsize=5)
    ax.tick_params(axis='y', pad=7)

plt.tight_layout()
plt.savefig(output + 'fig5_thermomodel_RNAP_emat_logo.pdf')


# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(seq)
im = ax.imshow(energy_df_scaled.T,
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'fig5_thermomodel_RNAP_ematonly.pdf')
