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
import scipy

# Import the project utils
import sys
sys.path.insert(0, '../')
import NB_sortseq_utils as utils

# Import matplotlib and stuff for plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize

utils.set_plotting_style4()

import seaborn as sns
sns.set_style('white')
#===============================================================================
# Set output directory based
#===============================================================================
output = 'output_figs/'
output2 = 'input_data/'
#===============================================================================
# define all details for the rel dataset
#===============================================================================
############## everything that needs to be changed goes here ##################

config.read('../sortseq/20150312_relB/cfg_files/20150312_relB_MG1655_M9glucose_na_mut1_4bins.cfg')
emat_choice = 'RNAP'

date = config.get('Input','date')
promoter = config.get('Input','promoter')
mutregion = config.get('Input', 'mutregion')
strain = config.get('Input','strain')
media = config.get('Input','media')
condition = config.get('Input','condition')
bincount = config.getint('Input','bincount')
Seq = config.get('Input','seq')

emat_dir_seq = '../../data/sortseq_pymc_dump/'
emat_dir_sql = '../../data/sortseq_MPAthic_MCMC/'
emat_dir_out = config.get('Input','emat_dir_out')
emat_fseqname = config.get('Input','emat_fseqname')

mut_region_start = config.getint(emat_choice,'mut_region_start')
mut_region_length = config.getint(emat_choice,'mut_region_length')
TF = config.get(emat_choice,'TF')
TF_type = config.getint(emat_choice,'TF_type')

# target file
fname_sql = date + '_' + promoter + '_' + strain + '_' + media + \
                   '_' + condition +  '_' + mutregion + '_' + str(bincount)\
                    + 'bins' + '_pymc_' + TF +'_MCMC'
fname_seq = date + '_' + promoter + '_' + strain + '_' + media + \
                   '_' + condition +  '_' + mutregion + '_' + str(bincount)\
                    + 'bins' + '_pymc'

burn_in = 1000

#------------------------------------------------------------------------------#
#load in the sequence data
#------------------------------------------------------------------------------#
data_fn = os.path.join(emat_dir_seq, emat_fseqname)
seq_mat, batch_vec = utils.load_unique_seqs_batches(data_fn, mut_region_start, mut_region_length)
energies = np.zeros(len(batch_vec))

#------------------------------------------------------------------------------#
# find sql MCMC trace files
#------------------------------------------------------------------------------#
listing_emat = glob.glob(os.path.join(emat_dir_sql, (fname_sql + '*.sql')))

#==============================================================================#
# plots for initialization of MCMC using sql MCMC trace file
#==============================================================================#

# initialize array to hold trace values and counter for finding average of
# energy matrices across MCMC traces
trace_emats = []
count = 0
for fn in listing_emat:
    run_number = fn[-8:-6]
    # if run_number != '02':
    #     continue
    db = pymc.database.sqlite.load(fn)

    emat_mean_temp = db.trace('emat')[0]
    emat_mean_temp = utils.fix_matrix_gauge(emat_mean_temp)

    if trace_emats == []:
        trace_emats = emat_mean_temp
        count += 1
    else:
        trace_emats += emat_mean_temp
        count += 1

#------------------------------------------------------------------------------#
# Step 5:  Calculate average energy matrix from all MCMC traces analyzed
# Combine all data in a tidy pandas and save csv
# Also calculate MI value and joint probability using everage energy matrix
#------------------------------------------------------------------------------#
emat_mean = trace_emats/count
emat_mean = utils.fix_matrix_gauge(emat_mean)


# # initialize a DataFrame to save the experimental details
df = pd.DataFrame()
pos = np.arange(mut_region_start,mut_region_start + mut_region_length)
# read the files and compute the mean YFP value
for i in enumerate(pos):
    try:
        df = df.append([[date, promoter, strain, media, condition, mutregion, \
                bincount, emat_choice, i[1]]], ignore_index=True)
    except:
        pass
#
# rename the columns of the data_frame
df.columns = ['date', 'promoter', 'strain', 'media', 'condition', 'mutregion', \
              'bincount', 'TF', 'position']

emat_df = pd.DataFrame(emat_mean.T,columns=list('ACGT'))

df_summary = pd.concat([emat_df,df], axis=1)

fname_out = date + '_' + promoter + '_' + strain + '_' + media + '_' + \
                             condition +  '_' + mutregion + '_' + str(bincount)\
                              + 'bins_' + emat_choice
df_summary.to_csv(output2 + 'figS10_SI_MCMCinference_intialize_matrix_mean.csv',
                              index=False)

# compute mutual information and joint pdf
MI,f_reg = utils.compute_MI(seq_mat,batch_vec,emat_mean)
print(MI, 'initialize')

fig, (ax1, ax2) = plt.subplots(2,1,figsize=utils.cm2inch((5,5)),
            gridspec_kw = {'height_ratios':[4, 1]})
ax1.imshow(f_reg,interpolation='nearest',aspect='auto', clim=(0.0, 0.001), cmap = 'Blues')
# ax1.set_xlabel('Rank order')
# ax1.set_ylabel('Batch number')
ax1.grid(b=False)
ylim = ax1.get_ylim()
print(ylim)
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
yticks = np.linspace(0,3,4)
ax1.set_yticks(yticks)
ax1.set_yticklabels(['1','2','3','4'], fontsize=8)

y = f_reg[3,:]
x = np.arange(0,len(y))
ax2.plot(x,y)
ax2.fill_between(x, 0, y, where=y >= 0, interpolate=True)
ax2.set_ylim(0,0.001)
ax2.set_yticks([])
ax2.set_axis_bgcolor('white')

for ax, color in zip([ax1, ax2], ['k', 'k']):
    plt.setp(ax.spines.values(), color=color, linewidth=0.3)

plt.tight_layout()
plt.savefig(output + 'figS10_SI_MCMCinference_intialize_jointDist_bin4.pdf')

#==============================================================================#
# plots using energy matrix at burn in point of MCMC
#==============================================================================#

# initialize array to hold trace values and counter for finding average of
# energy matrices across MCMC traces
trace_emats = []
count = 0
for fn in listing_emat:
    run_number = fn[-8:-6]
    # if run_number != '00':
    #     continue
    db = pymc.database.sqlite.load(fn)
    trace_temp = db.trace('emat')[burn_in]

    # Step 1: determine matrix entried, fix gauge, and check sign of energy values
    emat_mean_temp = trace_temp.copy()
    emat_mean_temp = utils.fix_matrix_gauge(emat_mean_temp)
    for i in range(len(batch_vec)):
        energies[i] = np.sum(seq_mat[:,:,i]*emat_mean_temp)
    r = scipy.stats.pearsonr(energies,batch_vec)[0]
    if r>0:
        trace_temp = -trace_temp * TF_type
        emat_mean_temp = -emat_mean_temp * TF_type
    else:
        trace_temp = trace_temp * TF_type
        emat_mean_temp = emat_mean_temp * TF_type

    if trace_emats == []:
        trace_emats = emat_mean_temp
        count += 1
    else:
        trace_emats += emat_mean_temp
        count += 1

    # If this is first loop, replace trace = [] with
    # trace_temp AND don't append energy matrix if the correlation r is poor.
    # I'm going to re-calculate the energy matrix based on sign corrected trace
    # and re-fix gauge.
    # emat_mean = trace_temp
    # emat_mean = utils.fix_matrix_gauge(emat_mean_temp)
    emat_mean = trace_emats/count
    emat_mean = utils.fix_matrix_gauge(emat_mean)


# # initialize a DataFrame to save the experimental details
df = pd.DataFrame()
pos = np.arange(mut_region_start,mut_region_start + mut_region_length)
# read the files and compute the mean YFP value
for i in enumerate(pos):
    try:
        df = df.append([[date, promoter, strain, media, condition, mutregion, \
                bincount, emat_choice, i[1]]], ignore_index=True)
    except:
        pass
#
# rename the columns of the data_frame
df.columns = ['date', 'promoter', 'strain', 'media', 'condition', 'mutregion', \
              'bincount', 'TF', 'position']

emat_df = pd.DataFrame(emat_mean.T,columns=list('ACGT'))

df_summary = pd.concat([emat_df,df], axis=1)

fname_out = date + '_' + promoter + '_' + strain + '_' + media + '_' + \
                             condition +  '_' + mutregion + '_' + str(bincount)\
                              + 'bins_' + emat_choice
df_summary.to_csv(output2 + 'figS10_SI_MCMCinference_burnin_matrix_mean.csv',
                              index=False)


# compute mutual information and joint pdf
MI,f_reg = utils.compute_MI(seq_mat,batch_vec,emat_mean)
print(MI, 'burn-in')

fig, (ax1, ax2) = plt.subplots(2,1,figsize=utils.cm2inch((5,5)),
            gridspec_kw = {'height_ratios':[4, 1]})
ax1.imshow(f_reg,interpolation='nearest',aspect='auto', clim=(0.0, 0.001), cmap = 'Blues')
# ax1.set_xlabel('Rank order')
# ax1.set_ylabel('Batch number')
ax1.grid(b=False)
ylim = ax1.get_ylim()
print(ylim)
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
yticks = np.linspace(0,3,4)
ax1.set_yticks(yticks)
ax1.set_yticklabels(['1','2','3','4'], fontsize=8)

y = f_reg[3,:]
x = np.arange(0,len(y))
ax2.plot(x,y)
ax2.fill_between(x, 0, y, where=y >= 0, interpolate=True)
ax2.set_ylim(0,0.001)
ax2.set_yticks([])
ax2.set_axis_bgcolor('white')

for ax, color in zip([ax1, ax2], ['k', 'k']):
    plt.setp(ax.spines.values(), color=color, linewidth=0.3)
    # plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=color)

plt.tight_layout()
fig.savefig(output + 'figS10_SI_MCMCinference_burnin_jointDist_bin4.pdf')
