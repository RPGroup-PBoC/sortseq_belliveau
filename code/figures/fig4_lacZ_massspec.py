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
files = glob.glob(datadir+'*_lacZ_*.csv')

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

# grab data for strain HG104 (wild-type copy number of LacI)
df = df[df.strain == 'HG104deltalysA']

# consider O3 data
df_O3 = df[df.promoter == 'lacZ_O3']

# determine 95% bounds on data
Y = df_O3['maxquant_ratio_medianshift'].dropna().values
x_2_5_O3 = np.log(np.percentile(Y,2.5))
x_97_5_O3 = np.log(np.percentile(Y,97.5))

# consider Oid data
df_Oid = df[df.promoter == 'lacZ_Oid']

# determine 95% bounds on data
Y = df_Oid['maxquant_ratio_medianshift'].dropna().values
x_2_5_Oid = np.log(np.percentile(Y,2.5))
x_97_5_Oid = np.log(np.percentile(Y,97.5))

#===============================================================================
# plot the enrichment ratios as scatter points
#===============================================================================

# consider only the proteins with predicted DNA binding motif for plotting.
# follow same procedure and find mean and standard deviation
df_TF_O3 = df_O3[df_O3['TF_check'] == 1]
df_TF_Oid = df_Oid[df_Oid['TF_check'] == 1]

# plot mean enrichment ratio and SEM for each protein
f, (ax1,ax2) = plt.subplots(1, 2, sharey=True, figsize=utils.cm2inch(4.5, 5))
f.subplots_adjust(wspace=0.1)

ax1.fill_between([0.95,1.05],np.exp(x_2_5_O3),np.exp(x_97_5_O3), color = 'k',alpha = 0.25)

yvals_TF_O3 = df_TF_O3.replace([-np.inf,np.inf], np.nan).dropna()

vals = np.random.uniform(0.99,1.01,size=len(yvals_TF_O3))
ax1.errorbar(vals, yvals_TF_O3['maxquant_ratio_medianshift'].values,
            linestyle='none',fmt='o', alpha = 0.6, markersize=4)
ax1.set_xlim(0.95,1.05)
ax1.set_xticklabels(['O3'])
ax1.set_ylabel('protein enrichment', fontsize=8)
ax1.xaxis.grid(False)

ax1.set_yscale('log')
ax1.tick_params(which='minor', length=2, color='#ffffff', direction = 'in')
ax1.set_ylim(1E-3,1E3)

ax2.fill_between([0.95,1.05],np.exp(x_2_5_Oid),np.exp(x_97_5_Oid), color = 'k',alpha = 0.25)

yvals_TF_Oid = df_TF_Oid.replace([-np.inf,np.inf], np.nan).dropna()
vals = np.random.uniform(0.99,1.01,size=len(yvals_TF_Oid))

ax2.errorbar(vals, yvals_TF_Oid['maxquant_ratio_medianshift'].values,
            linestyle='none',fmt='o', alpha = 0.6, markersize=4)
ax2.set_xlim(0.95,1.05)
ax2.set_xticklabels(['Oid'])
ax2.xaxis.grid(False)

ax2.set_yscale('log')
ax2.tick_params(which='minor', length=2, color='#ffffff', direction = 'in')
ax2.set_ylim(1E-2,1E3)

plt.tight_layout()


plt.savefig(output+'fig4_lacZ_massspec_scatter.pdf', format='pdf')
print('# with predected DNA-binding motif; O3',len(yvals_TF_O3))
print('# with predected DNA-binding motif; Oid',len(yvals_TF_Oid))
