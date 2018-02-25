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
from IPython.core.pylabtools import figsize

# Seaborn, useful for graphics
import seaborn as sns

sns.set_palette("deep", color_codes=True)
utils.set_plotting_style4()

#===============================================================================
# Set output directory
#===============================================================================
output = 'output_figs/'

#------------------------------------------------------------------------------#
# Load in the summary csv file
#------------------------------------------------------------------------------#

data_dir = '../sortseq/20160519_relB_sortingtests/'

fname = '20150519_relB_MG1655_M9glucose_15percentile_mut1_4bins_summary.csv'
df_1 = pd.read_csv(data_dir + fname)

fname = '20150519_relB_MG1655_M9glucose_22percentile_mut1_4bins_summary.csv'
df_2 = pd.read_csv(data_dir + fname)

fname = '20150519_relB_MG1655_M9glucose_8bins_mut1_8bins_summary.csv'
df_3 = pd.read_csv(data_dir + fname)

seqLength = len(df_1.position.unique())
#------------------------------------------------------------------------------#
# make shared axis plots of information footprint and delta bin shift
#------------------------------------------------------------------------------#
colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]

fig, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=utils.cm2inch(8,7))#, gridspec_kw = {'height_ratios':[1.75, 1]})
# make array for x positions - note that this will be relative to TSS
ind = np.arange(seqLength) - 35

# footprint 1
ax1.bar(ind, df_1.sort_values(by='position').MI, alpha=1, width=0.8, \
        linewidth = 0, color=colours[2])
ax1.set_xlim(ind[0],ind[-1])
ax1.set_ylabel('mutual\ninformation')
ax1.set_xlabel('position')
ax1.grid(b=False)

# # Hide the right and top spines
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
#
# # Only show ticks on the left and bottom spines
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')

#plt.setp(ax1.get_xticklabels(), visible=False)
# plt.setp(ax1.get_yticklabels(), visible=False)
ax1.set_yticklabels(['-','+'], fontname='Arial')

# footprint 2
ax2.bar(ind, df_2.sort_values(by='position').MI, alpha=1, width=0.8, \
        linewidth = 0, color=colours[2])
ax2.set_xlim(ind[0],ind[-1])
ax2.set_ylabel('mutual\ninformation')
ax2.set_xlabel('position')
ax2.grid(b=False)

# # Hide the right and top spines
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
#
# # Only show ticks on the left and bottom spines
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')

#plt.setp(ax2.get_xticklabels(), visible=False)
# plt.setp(ax2.get_yticklabels(), visible=False)
ax2.set_yticklabels(['-','+'], fontname='Arial')

# footprint 3
ax3.bar(ind, df_3.sort_values(by='position').MI, alpha=1, width=0.8, \
        linewidth = 0, color=colours[2])
ax3.set_xlim(ind[0],ind[-1])
ax3.set_ylabel('mutual\ninformation')
ax3.grid(b=False)

# # Hide the right and top spines
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
#
# # Only show ticks on the left and bottom spines
ax3.yaxis.set_ticks_position('left')
ax3.xaxis.set_ticks_position('bottom')

#plt.setp(ax3.get_xticklabels(), visible=False)
# plt.setp(ax3.get_yticklabels(), visible=False)


# make these tick labels invisible
#plt.setp(ax1.get_xticklabels(), fontsize=6)
# plt.setp(ax3.get_yticklabels(), visible=False)
ax3.set_yticklabels(['-','+'], fontname='Arial')
ax3.set_xlabel('position')

fig.subplots_adjust(hspace=-1)
plt.tight_layout()

figname_out = 'figS3_relB_sorting_MI.pdf'
fig.savefig(output + figname_out, format='pdf')
