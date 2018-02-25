import os
import glob
import pickle
import re
import pymc

import ConfigParser
config = ConfigParser.RawConfigParser()

# Our numerical workhorses
import numpy as np
import pandas as pd
import scipy as sp

# Import the project utils
import sys
sys.path.insert(0, '../')
import NB_sortseq_utils as utils

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize

# Seaborn, useful for graphics
import seaborn as sns

utils.set_plotting_style1()

from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size


#===============================================================================
# Set output directory based
#===============================================================================
def calc_MI(emat_0,data_fn):
    """ Function to estimate mutual information between binned sequences and
    their energy matrix predictions.
    """
    # load in the data
    seq_mat_temp, batch_vec_temp = utils.load_seqs_batches(data_fn,mut_region_start,mut_region_length)
    # shuffle the elements of seq_mat and batch_vec. This will prevent
    # spuriously high mutual information values
    index_shuf = range(len(batch_vec_temp))
    sp.random.shuffle(index_shuf)
    seq_mat = sp.zeros_like(seq_mat_temp)
    batch_vec = sp.zeros_like(batch_vec_temp)
    for i, i_s in enumerate(index_shuf):
        seq_mat[:,:,i] = seq_mat_temp[:,:,i_s]
        batch_vec[i] = batch_vec_temp[i_s]

    s=seq_mat
    b=batch_vec
    value=emat_0

    n_seqs = len(b)
    MI, f_reg = utils.compute_MI(s,b,value)
    return MI, f_reg

#===============================================================================
# Set output directory based
#===============================================================================
output = 'output_figs/'

#===============================================================================
#load in the data for rel
#===============================================================================
mut_region_start = 31
mut_region_length = 23
data_fn_delta = 'input_data/20150513_relB_MG1655deltarelBE_M9glucose_na_mut1_4bins_pymc_slim.csv'
data_fn_WT_replicate = 'input_data/20150519_relB_MG1655_M9glucose_15percentile_mut1_4bins_pymc_slim.csv'

data_rel_emat = '../sortseq/20150312_relB/20150513_relB_MG1655_M9glucose_na_mut1_4bins_RelBE_emat_mean.csv'
energy_df = pd.read_csv(data_rel_emat)
energy_df = energy_df[['A','C','G','T']].T
emat_RelBE = np.array(energy_df)

data_relDelta_emat = '../sortseq/20150513_relBWTRNAP_relBdeltaBE/20150513_relB_MG1655deltarelBE_M9glucose_na_mut1_4bins_RelBE_emat_mean.csv'
energy_df_delta = pd.read_csv(data_relDelta_emat)
energy_df_delta = energy_df_delta[['A','C','G','T']].T
emat_deltaRelBE = np.array(energy_df_delta)

#==============================================================================#
# Joint distribution plots
#==============================================================================#

#------------------------------------------------------------------------------#
# wild type data (replicate from experiments testing different sorting)
#------------------------------------------------------------------------------#

# calculate mutual information and joint distribution
MI, f_reg = calc_MI(emat_RelBE,data_fn_WT_replicate)

fig, (ax1) = plt.subplots(1,1,figsize=(5,4))
im = ax1.imshow(f_reg,interpolation='nearest',aspect='auto', clim=(0.0002,0.0003), cmap = 'Blues')
ax1.set_xlabel('rank order of energy prediction')
ax1.set_ylabel('bin number')
ax1.grid(b=False)
ylim = ax1.get_ylim()
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
yticks = np.linspace(0,3,4)
ax1.set_yticks(yticks)
ax1.set_yticklabels(['1','2','3','4'])

# create an axes on the right side of ax. The width of cax will be 3%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="3%", pad=0.05)
cbar = plt.colorbar(im, cax=cax, ticks=[0.0002,0.0003])
cbar.ax.set_yticklabels(['0.0002', '0.00034'], fontname='Arial')
cbar.outline.set_visible(False)

for ax, color in zip([ax1], ['k']):
    plt.setp(ax.spines.values(), color=color, linewidth=0.3)

plt.tight_layout()
plt.savefig(output + 'figS7_SI_wildtype_replicate_jointDist_bin4.pdf')

#------------------------------------------------------------------------------#
# wild type data (replicate from experiments testing different sorting)
#------------------------------------------------------------------------------#

# calculate mutual information and joint distribution
MI_delta_model, f_reg = calc_MI(emat_RelBE,data_fn_delta)

fig, (ax1) = plt.subplots(1,1,figsize=(5,4))
im = ax1.imshow(f_reg,interpolation='nearest',aspect='auto', clim=(0.0002,0.0003), cmap = 'Blues')
ax1.set_xlabel('rank order of energy prediction')
ax1.set_ylabel('bin number')
ax1.grid(b=False)
ylim = ax1.get_ylim()
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
yticks = np.linspace(0,3,4)
ax1.set_yticks(yticks)
ax1.set_yticklabels(['1','2','3','4'])

# create an axes on the right side of ax. The width of cax will be 3%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="3%", pad=0.05)
cbar = plt.colorbar(im, cax=cax, ticks=[0.0002,0.0003])
cbar.ax.set_yticklabels(['0.0002', '0.00034'], fontname='Arial')
cbar.outline.set_visible(False)

for ax, color in zip([ax1], ['k']):
    plt.setp(ax.spines.values(), color=color, linewidth=0.3)

plt.tight_layout()
plt.savefig(output + 'figS7_SI_delta_relBE_jointDist_bin4.pdf')


#==============================================================================#
# Now lets only consider delta data
#==============================================================================#

#------------------------------------------------------------------------------#
# estimate maximum information (using matrix generated with delta relBE data)
#------------------------------------------------------------------------------#
MI_max_delta_relBE, f_reg = calc_MI(emat_deltaRelBE,data_fn_delta)

#------------------------------------------------------------------------------#
# estimate maximum information from 20 random matrices (gauge fixed for
# propter comparison with other matrices)
#------------------------------------------------------------------------------#
# # make random matrices
# MI_rand_delta_data = np.zeros(20)
#
# for i in MI_rand_delta_data:
#     emat_rand = np.random.normal(size=[4, 23])
#     emat_rand = utils.fix_matrix_gauge(emat_rand)
#     MI_rand_delta_data[i], f_reg = calc_MI(emat_rand,data_fn_delta)

MI_rand_delta_data = [0.000203556344573, 0.000286124670955, 9.03629024909e-05,
                      4.65683350724e-05, 0.000249101783341, 0.000475455457806,
                      0.00017929444908, 9.76428055725e-05, 8.49672871268e-05,
                      7.88130373761e-05, 0.000259838937893, 0.000186330367669,
                      0.000215818418913, 0.000120932149623, 0.000132778911397,
                      0.00067400017891, 0.000244297353263, 5.70352349856e-05,
                      0.000165716920104, 9.79446090617e-05]

#------------------------------------------------------------------------------#
# plot MI values for comparisons against delta relBE data
#------------------------------------------------------------------------------#
# comparing relBE matrices to delta data

Info_y = np.array([MI_max_delta_relBE, MI_delta_model, np.mean(MI_rand_delta_data)])*1000
Info_yerr = np.array([0.0, 0.0, np.std(MI_rand_delta_data)])*1000

fig1 = plt.figure(figsize = (5,4))

plt.bar(np.arange(0,3),Info_y, yerr=Info_yerr, ecolor='r')
plt.xticks([])
plt.yticks(np.arange(0, 6,2))

plt.ylabel('mutual information\n(mbits)')
plt.tight_layout()

fig1.savefig(output + 'figS7_SI_relBE_matrix_compare.pdf')


#==============================================================================#
# Now lets switch to mar promoter; only consider delta data
#==============================================================================#

mut_region_start = 11
mut_region_length = 14
data_fn_WT = 'input_data/20150513_marR_MG1655_LB_na_mut1_4bins_pymc.csv'
data_fn_delta = 'input_data/20150513_marR_MG1655deltamarR_LB_na_mut1_4bins_pymc.csv'

data_MarR_emat = '../sortseq/20150513_marRmut1only_marRdeltaRAB_marRdeltaR/20150513_marR_MG1655_LB_na_mut1_4bins_MarR_left_emat_mean.csv'
energy_df_marR = pd.read_csv(data_MarR_emat)
energy_df_marR = np.array(energy_df_marR[['A','C','G','T']].T)[:,:-2]

data_MarRDelta_emat = '../sortseq/20150513_marRmut1only_marRdeltaRAB_marRdeltaR/20150513_marR_MG1655deltamarR_LB_na_mut1_4bins_MarR_left_emat_mean.csv'
energy_df_marR_delta = pd.read_csv(data_MarRDelta_emat)
energy_df_marR_delta = np.array(energy_df_marR_delta[['A','C','G','T']].T)[:,:-2]

#------------------------------------------------------------------------------#
# estimate maximum information (using matrix generated with delta relBE data)
#------------------------------------------------------------------------------#
MI_max_delta_marR, f_reg = calc_MI(energy_df_marR_delta,data_fn_delta)

#------------------------------------------------------------------------------#
# estimate maximum information from 20 random matrices (gauge fixed for
# propter comparison with other matrices)
#------------------------------------------------------------------------------#
# # make random matrices
# MI_rand_delta_data_marR = np.zeros(20)
#
# for i in MI_rand_delta_data_marR:
#     emat_rand = np.random.normal(size=[4, 14])
#     emat_rand = utils.fix_matrix_gauge(emat_rand)
#     MI_rand_delta_data_marR[i], f_reg = calc_MI(emat_rand,data_fn_delta)

MI_rand_delta_data_marR = [0.000103859348938, 6.76679274067e-05, 6.03309979553e-05,
                           6.01277720217e-05, 9.0278879905e-05, 7.29150402355e-05,
                           5.43521103489e-05, 2.08105431957e-05, 5.26956831726e-05,
                           4.69295707964e-05, 0.000127103831841, 6.30123212317e-05,
                           7.12755765161e-05, 5.70322290925e-05, 8.79042238412e-05,
                           7.29015868737e-05, 5.80020209679e-05, 8.53626007545e-05,
                           0.000106762561939, 4.2431094977e-05]

#------------------------------------------------------------------------------#
# estimate mutual information between marR matrix (left binding site
# in between -10 and -35 of RNAP)
#------------------------------------------------------------------------------#
MI_marR_model, f_reg = calc_MI(energy_df_marR,data_fn_delta)

#------------------------------------------------------------------------------#
# plot MI values for comparisons against delta marR data
#------------------------------------------------------------------------------#
# comparing marR matrices to delta data

Info_y = 1000*np.array([MI_max_delta_marR, MI_marR_model, np.mean(MI_rand_delta_data_marR)])
Info_yerr = 1000*np.array([0.0, 0.0, np.std(MI_rand_delta_data_marR)])

fig1 = plt.figure(figsize = (5,4))

plt.bar(np.arange(0,3),Info_y, yerr=Info_yerr, ecolor='r')
plt.xticks([])
plt.yticks(np.arange(0, 2,0.5))

plt.ylabel('mutual information\n(mbits)')
plt.tight_layout()

fig1.savefig(output + 'figS7_SI_marR_matrix_compare.pdf')
