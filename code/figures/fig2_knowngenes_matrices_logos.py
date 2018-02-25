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

# Import matplotlib stuff for plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Logo-generating module and utils
import anylogo
import NB_sortseq_utils as utils

# set plotting format options
utils.set_plotting_style_emat()

#===============================================================================
# Set output directory
#===============================================================================
output = 'output_figs/'

#===============================================================================
# directory where emat csv files are contained
#===============================================================================

# lac
datadir_lac = '../sortseq/2011_lacZ/'
# mar
datadir_marA = '../sortseq/20150820_marRmut2only/'
datadir_mar_RNAP = '../sortseq/20150513_marRmut1only_marRdeltaRAB_marRdeltaR/'
# rel
datadir_rel = '../sortseq/20150312_relB/'

#===============================================================================
# plot energy matrices with logos on top.
#===============================================================================

# Set color scale - I want the colorbar to be symmetric and will pick values#
# that seem appropriate for all matrices.
emat_min=-0.4
emat_max=0.4
mid_val=0.0

# Create background dict
gc = .508
background_array =pd.DataFrame( [[(1-gc)/2,gc/2,gc/2,(1-gc)/2]])
#------------------------------------------------------------------------------#
# lacZ: LacI
#------------------------------------------------------------------------------#

energy_df = pd.read_csv(datadir_lac + '2011_lacZ_MG1655_M9glucose_na_mut1_4bins_LacI_O1_emat_mean.csv')
energy_df = energy_df[['A','C','G','T']]

# create background nucleotide frequencies dataframe
energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])

seq = 'AATTGTGAGCGGATAACAATT'

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
plt.savefig(output + 'fig2_lacZ_emat_logo_lacI.pdf')

# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(seq)
im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), seq),
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'fig2_lacZ_emat_logo_lacI_ematonly.pdf')

#------------------------------------------------------------------------------#
# logos for CRP
#------------------------------------------------------------------------------#

energy_df = pd.read_csv(datadir_lac + '2011_lacZ_MG1655_M9glucose_na_mut1_4bins_CRP_emat_mean.csv')

energy_df = energy_df[['A','C','G','T']]

energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
energy_df_scaled = energy_df_scaled[['A','C','G','T']]

# create background nucleotide frequencies dataframe
energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])

seq = 'ATTAATGTGAGTTAGCTCACTCATTA'

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
plt.savefig(output + 'fig2_lacZ_emat_logo_CRP.pdf')

# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(seq)
im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), seq),
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'fig2_lacZ_emat_logo_CRP_ematonly.pdf')

#------------------------------------------------------------------------------#
# lacZ: RNAP
#------------------------------------------------------------------------------#

energy_df = pd.read_csv(datadir_lac + '2011_lacZ_MG1655_M9glucose_na_mut1_4bins_RNAP_emat_mean.csv')

# energy_df['position'] = energy_df['position'] -  63
energy_df = energy_df[energy_df.position != energy_df.position.min()]
energy_df = energy_df[energy_df.position != energy_df.position.max()]

energy_df.reset_index(inplace=True)
energy_df = energy_df[['A','C','G','T']]

energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()

energy_df_scaled.reset_index(inplace=True)
energy_df_scaled = energy_df_scaled[['A','C','G','T']]

# create background nucleotide frequencies dataframe
energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])

seq = 'TTTACACTTTATGCTTCCGGCTCGTATGTT'

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
plt.savefig(output + 'fig2_lacZ_emat_logo_RNAP.pdf')

# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(seq)
im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), seq),
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'fig2_lacZ_emat_logo_RNAP_ematonly.pdf')

#------------------------------------------------------------------------------#
# marRAB: MarA
#------------------------------------------------------------------------------#

energy_df = pd.read_csv(datadir_marA + '20150820_marR_MG1655_LB_na_mut2_4bins_MarA_emat_mean.csv')
# energy_df = energy_df[energy_df.position != energy_df.position.max()]
energy_df = energy_df[['A','C','G','T']]

energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
energy_df_scaled = energy_df_scaled[['A','C','G','T']]

# create background nucleotide frequencies dataframe
energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])


seq = 'ATTTAGCAAAACGTGGCATC'

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
plt.savefig(output + 'fig2_marRAB_emat_logo_marA.pdf')


# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(seq)
im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), seq),
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'fig2_marRAB_emat_logo_marA_ematonly.pdf')

#------------------------------------------------------------------------------#
# marRAB: RNAP
#------------------------------------------------------------------------------#

energy_df = pd.read_csv(datadir_mar_RNAP + '20150513_marR_MG1655_LB_na_mut1_4bins_RNAP_emat_mean.csv')
energy_df = energy_df[energy_df.position != energy_df.position.max()]
energy_df = energy_df[energy_df.position != energy_df.position.max()]
energy_df.reset_index(inplace=True)
energy_df = energy_df[['A','C','G','T']]

energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
energy_df_scaled.reset_index(inplace=True)
energy_df_scaled = energy_df_scaled[['A','C','G','T']]

# create background nucleotide frequencies dataframe
energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])

seq = 'TTGACTTATACTTGCCTGGGCAATATTAT'

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
plt.savefig(output + 'fig2_marRAB_emat_logo_RNAP.pdf')

# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(seq)
im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), seq),
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'fig2_marRAB_emat_logo_RNAP_ematonly.pdf')

#------------------------------------------------------------------------------#
# relB promoter: RelBE
#------------------------------------------------------------------------------#

energy_df = pd.read_csv(datadir_rel + '20150513_relB_MG1655_M9glucose_na_mut1_4bins_RelBE_emat_mean.csv')
energy_df = energy_df[['A','C','G','T']]

energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
energy_df_scaled = energy_df_scaled[['A','C','G','T']]

# create background nucleotide frequencies dataframe
energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])

seq = 'TGTAATGACATTTGTAATTACAA'

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
plt.savefig(output + 'fig2_relB_emat_logo_RelBE.pdf')

# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(seq)
im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), seq),
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'fig2_relB_emat_logo_RelBE_ematonly.pdf')

#------------------------------------------------------------------------------#
# relB promoter: RelBE
#------------------------------------------------------------------------------#

energy_df = pd.read_csv(datadir_rel + '20150513_relB_MG1655_M9glucose_na_mut1_4bins_RNAP_emat_mean.csv')
energy_df = energy_df[energy_df.position != energy_df.position.min()]
energy_df = energy_df[energy_df.position != energy_df.position.max()]
energy_df.reset_index(inplace=True)
energy_df = energy_df[['A','C','G','T']]

energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
energy_df_scaled.reset_index(inplace=True)
energy_df_scaled = energy_df_scaled[['A','C','G','T']]

# create background nucleotide frequencies dataframe
energy_df_scaled = utils.estimate_scalefactor(np.array(energy_df))*energy_df.copy()
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])

seq = 'TTGCCCTAAGCATGTTGTAGTGCGATACT'

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
plt.savefig(output + 'fig2_relB_emat_logo_RNAP.pdf')

# save energy matrix using nearest interpolation
plt.figure()
ax = plt.gca()
L = len(seq)
im = ax.imshow(utils.zero_matrix_WT(np.array(energy_df.T), seq),
            interpolation='nearest',
            cmap='RdBu_r',
            clim=(emat_min, emat_max),
            norm = utils.MidpointNormalize(midpoint = mid_val,
                    vmin = emat_min, vmax = emat_max))
ax.axis('off')

plt.savefig(output + 'fig2_relB_emat_logo_RNAP_ematonly.pdf')
