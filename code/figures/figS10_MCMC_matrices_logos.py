import os
import glob
import pickle
import re


# Our numerical workhorses
import numpy as np
import pandas as pd
import scipy

# Import the project utils
import sys
sys.path.insert(0, '../')
import NB_sortseq_utils as utils

import anylogo

# Import matplotlib and stuff for plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize

utils.set_plotting_style4()

import seaborn as sns
sns.set_style('white')


import matplotlib.colors as colors

# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
     http://chris35wills.github.io/matplotlib_diverging_colorbar/
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
#===============================================================================
# Set output directory based
#===============================================================================
output = 'output_figs/'
input_data = 'input_data/'

#==============================================================================#
# plots for initialization of MCMC using sql MCMC trace file
#==============================================================================#

energy_df = pd.read_csv(input_data + 'figS10_SI_MCMCinference_intialize_matrix_mean.csv')
energy_df = energy_df[['A','C','G','T']]

energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()

RNAP_seq = 'CTTGCCCTAAGCATGTTGTAGTGCGATACTT'
# emat_min=emat.min()
emat_min=-0.5
# emat_max=emat.max()
emat_max=0.5
mid_val=0.0

plt.figure(figsize=utils.cm2inch((6.5,3)))
ax = plt.gca()
relative_scale=1
relative_spacing=.45
emat_ymin = -2 * (relative_scale + relative_spacing)
emat_ymax = -2 * relative_spacing
yticks = np.linspace(emat_ymin, emat_ymax, 9)[[1, 3, 5, 7]]
yticklabels = list('TGCA')
anylogo.draw(ax, effect_df=energy_df_scaled, logo_type='information',
             use_transparency=False)
L = len(energy_df)
ax.set_xticks([])

im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), RNAP_seq),
			   interpolation='nearest',cmap='RdBu_r',
			   clim=(emat_min, emat_max),
			   norm=MidpointNormalize(midpoint=mid_val,
			   			vmin=emat_min, vmax=emat_max),
						extent=(-.5, L - .5, emat_ymin, emat_ymax),
			   zorder=100, aspect='auto')

ax.set_ylim([emat_ymin, 2])
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels, fontsize=6, horizontalalignment='center')
ax.set_ylabel('')
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.yaxis.set_tick_params(length=0)


# create an axes on the right side of ax. The width of cax will be 3%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.05)

cbar = plt.colorbar(im, cax=cax, ticks=[-0.5, 0, 0.5])
cbar.ax.set_yticklabels(['-0.5', '0', '0.5'], fontsize=6, fontname='Arial')
cbar.outline.set_visible(False)


y = .5*emat_ymax
for i in range(L):
    ax.text(i, y, RNAP_seq[i], horizontalalignment='center', verticalalignment='center',
    fontsize=6)
    ax.tick_params(axis='y', pad=7)

# # make heat map of mean matrix
# plt.figure()
# im = plt.imshow(utils.zero_matrix(emat_mean),interpolation='nearest')
# ax = plt.gca()
# ax.set_axis_off()
# plt.title('Energy Matrix, ' + ', MI: %.5f' % MI)
plt.tight_layout()
plt.savefig(output + 'figS10_SI_MCMCinference_intialize_matrix_mean_logo.pdf')

# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(RNAP_seq)
im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), RNAP_seq),
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'figS10_SI_MCMCinference_intialize_ematonly.pdf')


#==============================================================================#
# plots using energy matrix at burn in point of MCMC
#==============================================================================#
#==============================================================================#
# plots for initialization of MCMC using sql MCMC trace file
#==============================================================================#

energy_df = pd.read_csv(input_data + 'figS10_SI_MCMCinference_burnin_matrix_mean.csv')
energy_df = energy_df[['A','C','G','T']]

energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()

RNAP_seq = 'CTTGCCCTAAGCATGTTGTAGTGCGATACTT'
# emat_min=emat.min()
emat_min=-0.5
# emat_max=emat.max()
emat_max=0.5
mid_val=0.0

plt.figure(figsize=utils.cm2inch((6.5,3)))
ax = plt.gca()
relative_scale=1
relative_spacing=.45
emat_ymin = -2 * (relative_scale + relative_spacing)
emat_ymax = -2 * relative_spacing
yticks = np.linspace(emat_ymin, emat_ymax, 9)[[1, 3, 5, 7]]
yticklabels = list('TGCA')
anylogo.draw(ax, effect_df=energy_df_scaled, logo_type='information',
             use_transparency=False)
# anylogo.draw(ax, effect_df=energy_df, logo_type='information', find_beta=True,
#              use_transparency=False)
L = len(energy_df)
ax.set_xticks([])

im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), RNAP_seq),interpolation='nearest',cmap='RdBu_r',clim=(emat_min, emat_max), norm=MidpointNormalize(midpoint=mid_val,vmin=emat_min, vmax=emat_max), extent=(-.5, L - .5, emat_ymin, emat_ymax), zorder=100, aspect='auto')

ax.set_ylim([emat_ymin, 2])
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels, fontsize=6, horizontalalignment='center')
ax.set_ylabel('')
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.yaxis.set_tick_params(length=0)


# create an axes on the right side of ax. The width of cax will be 3%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.05)

cbar = plt.colorbar(im, cax=cax, ticks=[-0.5, 0, 0.5])
cbar.ax.set_yticklabels(['-0.5', '0', '0.5'], fontsize=6, fontname='Arial')
cbar.outline.set_visible(False)


y = .5*emat_ymax
for i in range(L):
    ax.text(i, y, RNAP_seq[i], horizontalalignment='center', verticalalignment='center',
    fontsize=6)
    ax.tick_params(axis='y', pad=7)

# # make heat map of mean matrix
# plt.figure()
# im = plt.imshow(utils.zero_matrix(emat_mean),interpolation='nearest')
# ax = plt.gca()
# ax.set_axis_off()
# plt.title('Energy Matrix, ' + ', MI: %.5f' % MI)
plt.tight_layout()
plt.savefig(output + 'figS10_SI_MCMCinference_burnin_matrix_mean_logo.pdf')

# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(RNAP_seq)
im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), RNAP_seq),
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'figS10_SI_MCMCinference_burnin_ematonly.pdf')
