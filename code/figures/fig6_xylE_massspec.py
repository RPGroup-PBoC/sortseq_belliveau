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
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Seaborn, useful for graphics
import seaborn as sns

sns.set_palette("deep", color_codes=True)

utils.set_plotting_style_MS()

#===============================================================================
# Set output directory based on the graphicspath.tex file to print in dropbox
#===============================================================================
output = 'output_figs/'

#===============================================================================
# Read the data
#===============================================================================

datadir = '../mass_spec/*/'
files = glob.glob(datadir+'*_xylE_longDNA*.csv')

df = pd.DataFrame()

for f in enumerate(files):
    print(f[1])
    df_temp = pd.DataFrame()
    df_temp = pd.read_csv(f[1])

    # append data to df
    df = df.append(df_temp)


#===============================================================================
# determine 95% probability density bounds using all data
#===============================================================================

# calculate the average and standard deviation of each protein detected
by_row_index = df.groupby(df.gene)
df_means = by_row_index['maxquant_ratio_medianshift'].mean()
df_std = by_row_index['maxquant_ratio_medianshift'].std()

# convert the groupby arrays back into dataframes
df_means_new =  pd.DataFrame(df_means).reset_index()
df_means_new.columns = ['gene', 'maxquant_ratio_normalized_avg']
df_std_new = pd.DataFrame(df_std).reset_index()
df_means_new['maxquant_ratio_normalized_std'] = df_std_new['maxquant_ratio_medianshift']

# determine 95% bounds on data
Y = df_means_new['maxquant_ratio_normalized_avg'].dropna().values
x_2_5 = np.log(np.percentile(Y,2.5))
x_97_5 = np.log(np.percentile(Y,97.5))


#===============================================================================
# remove proteins which were not identified in each replicate.
# for identification of TFs we will consider those proteins which were
# found in each of the three replicates.
#===============================================================================

# group by gene names and remove those which do not have a ratio in
# each replicate
by_row_index = df.groupby(df.gene)

# discard enrichment ratios not based on three measurements
for i in by_row_index['maxquant_ratio_medianshift']:
    if len(i[1].dropna()) < 3:
        df = df[df.gene != i[0]]

#===============================================================================
# plot the enrichment ratios as scatter points
#===============================================================================

#-------------------------------------------------------------------------------
# consider only the proteins with predicted DNA binding motif for plotting.
# follow same procedure and find mean and standard deviation
#-------------------------------------------------------------------------------
# make DataFrame with only proteins with known or predicted DNA binding domain
df_TF = df[df['TF_check'] == 1]
# group data by protein and calculate average log ratio and std. deviation
by_row_index = df_TF.groupby(df_TF.gene)
df_means_TF_log = by_row_index['maxquant_logratio_medianshift'].mean()
df_std_TF_log = by_row_index['maxquant_logratio_medianshift'].std()

# Make DataFrame of mean values
df_means_new_TF =  pd.DataFrame(df_means_TF_log).reset_index()
df_means_new_TF.columns = ['gene', 'maxquant_ratio_normalized_avg_log']

# Make DataFrame of std. deviation values
df_std_new_TF =  pd.DataFrame(df_std_TF_log).reset_index()
df_std_new_TF.columns = ['gene', 'maxquant_ratio_normalized_std_log']

# Merge average and std. deviation values into summary DataFrame
df_summary_TF = pd.merge(df_means_new_TF,df_std_new_TF,on='gene')

#-------------------------------------------------------------------------------
# calculate average enrichment ratio (not log) as well as lower and upper bounds
# based on log std. deviation.
#-------------------------------------------------------------------------------

df_summary_TF['maxquant_ratio_normalized_avg'] = \
            np.exp(df_summary_TF['maxquant_ratio_normalized_avg_log'])

# determine error bar - lower bound
df_summary_TF['bottom'] = np.exp(df_summary_TF['maxquant_ratio_normalized_avg_log'] -
                          df_summary_TF['maxquant_ratio_normalized_std_log']/np.sqrt(3))
# determine error bar - upper bound
df_summary_TF['top'] = np.exp(df_summary_TF['maxquant_ratio_normalized_avg_log'] +
                       df_summary_TF['maxquant_ratio_normalized_std_log']/np.sqrt(3))

#-------------------------------------------------------------------------------
# plot mean enrichment ratio and standard deviation for each protein
#-------------------------------------------------------------------------------

fig = plt.figure(figsize=utils.cm2inch(3, 5))
ax1 = fig.add_subplot(111)

ax1.fill_between([0.95,1.05],np.exp(x_2_5),np.exp(x_97_5), color = 'k',alpha = 0.25)

yvals_TF = df_summary_TF
vals = np.random.uniform(0.99,1.01,size=len(yvals_TF))

yvals_err = [df_summary_TF['maxquant_ratio_normalized_avg'] - df_summary_TF['bottom'],
            df_summary_TF['top'] - df_summary_TF['maxquant_ratio_normalized_avg']]

ax1.errorbar(vals, yvals_TF['maxquant_ratio_normalized_avg'],
             yerr = yvals_err,
            linestyle='none',fmt='o', alpha = 0.6, markersize=4)
ax1.set_xlim(0.95,1.05)
ax1.set_xticklabels(['xylE'])
ax1.set_ylabel('protein enrichment', fontsize=8)
ax1.xaxis.grid(False)

ax1.set_yscale('log')
ax1.tick_params(which='minor', length=2, color='#ffffff', direction = 'in')
ax1.set_ylim(1E-3,1E3)
plt.tight_layout()

plt.savefig(output+'fig6_xylE_massspec_scatter.pdf', format='pdf')
print('# with predected DNA-binding motif',len(yvals_TF))
