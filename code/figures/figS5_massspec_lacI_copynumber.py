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
# Set output directory
#===============================================================================

output = 'output_figs/'

#===============================================================================
# Read the data
#===============================================================================


datadir = '../mass_spec/*/'
files = glob.glob(datadir + '*.csv')

df = pd.DataFrame()

for f in enumerate(files):
    # if "titration" in f[1]:
    #     continue
    # grab only the lac data.
    if "lac" not in f[1]:
        continue
    df_temp = pd.DataFrame()
    df_temp = pd.read_csv(f[1])
    # for those with reverse labeling, take inver ratio as ratio
    if df_temp['labeling'].unique() == 'reverse':
        df_temp['maxquant_ratio_medianshift'] = 1.0 / df_temp['maxquant_ratio_medianshift']
        df_temp['maxquant_logratio_medianshift'] = - df_temp['maxquant_logratio_medianshift']

    # append data to df
    df = df.append(df_temp)

#===============================================================================
# plot the data
#===============================================================================

fig = plt.figure(figsize=utils.cm2inch(3, 5))
ax1 = fig.add_subplot(111)

# copy numbers associated with each strain (number of tetramers per cell)
R_copy = {'HG104deltalysA': 11,
          'HG105deltalysA_ybcN<>3*1_RBS1027_lacI': 130,
          'HG105deltalysA_ybcN<>3*1lacI': 870}

df_lacI = df[(df.gene == 'lacI') & (df.strain != 'MG1655deltalysA')]
op_color  = {'lacZ_O3': 'b', 'lacZ_Oid': 'r'}
op_name  = {'lacZ_O3': 'O3', 'lacZ_Oid': 'Oid'}
marker_style = {'HG104deltalysA': 'o',
                'HG105deltalysA_ybcN<>3*1_RBS1027_lacI': 's',
                'HG105deltalysA_ybcN<>3*1lacI': 'D'}

strains = ['HG104deltalysA', 'HG105deltalysA_ybcN<>3*1_RBS1027_lacI', 'HG105deltalysA_ybcN<>3*1lacI']
check = []

for op in df_lacI['promoter'].unique():
    for strain in strains:
        for ratio in df_lacI['maxquant_ratio_medianshift'][(df_lacI.strain == strain) & (df_lacI.promoter == op)]:
            # use if statement so that label is only created once for each
            # strain.
            if check == [strain, op]:
                plt.errorbar(1, ratio, color = op_color[op], markersize = 4,
                       marker = marker_style[strain], alpha = 0.6)
            else:
                plt.errorbar(1, ratio, color = op_color[op], markersize = 4,
                            label = str(op_name[op]) + ', R = ' +
                            str(R_copy[strain]) + ' LacI / cell',
                            marker = marker_style[strain], alpha = 0.6,
                            linewidth=0)
            check = [strain, op]

# ax1.set_ylabel('protein enrichment', fontsize=8)
ax1.xaxis.grid(False)

ax1.set_yscale('log')

ax1.tick_params(which='minor', length=2, color='#ffffff', direction = 'in')

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 8)
ax1.set_ylim(1E0,1E3)
ax1.xaxis.grid(False)
ax1.set_xticklabels([])
ax1.set_ylabel('protein enrichment', fontsize=8)

plt.tight_layout()
plt.savefig(output + 'figS5_massspec_lacI_copynumber.pdf', format='pdf',
        bbox_extra_artists=(lgd,), bbox_inches='tight')
