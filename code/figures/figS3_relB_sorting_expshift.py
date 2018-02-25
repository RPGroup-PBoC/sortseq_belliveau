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
ax1.bar(ind, df_1.sort_values(by='position').expshift, alpha=1, width=0.8, \
        linewidth = 0, color='k')
ax1.set_xlim(ind[0],ind[-1])
ax1.set_ylabel('expression\nshift')
ax1.set_xlabel('position')
ax1.grid(b=False)
ax1.set_facecolor('white')

# Hide the right and top spines
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)

# Annotate y axis
ylim = ax1.get_ylim()
#ax.set_yticks([])
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
ax1.set_yticks(yticks)
# fontProperties = {'family':'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 14}
ax1.set_yticklabels(['-','+'], fontname='Arial')
ax1.set_ylabel('expression\nshift')

# Make grid lines
xticks = [x for x in ind if x%10 == 0]
for n, x in enumerate(xticks):
        ax1.axvline(x, color='lightgray', linewidth=0.36, zorder=-1)


# footprint 2
ax2.bar(ind, df_2.sort_values(by='position').expshift, alpha=1, width=0.8, \
        linewidth = 0, color='k')
ax2.set_xlim(ind[0],ind[-1])
ax2.set_ylabel('expression\nshift')
ax2.set_xlabel('position')
ax2.grid(b=False)
ax2.set_facecolor('white')

# Hide the right and top spines
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)

# Only show ticks on the left and bottom spines
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')

# Annotate y axis
ylim = ax2.get_ylim()
#ax.set_yticks([])
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
ax2.set_yticks(yticks)
# fontProperties = {'family':'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 14}
ax2.set_yticklabels(['-','+'], fontname='Arial')
ax2.set_ylabel('expression\nshift')

# Make grid lines
xticks = [x for x in ind if x%10 == 0]
for n, x in enumerate(xticks):
        ax2.axvline(x, color='lightgray', linewidth=0.36, zorder=-1)
# footprint 3
ax3.bar(ind, df_3.sort_values(by='position').expshift, alpha=1, width=0.8, \
        linewidth = 0, color='k')
ax3.set_xlim(ind[0],ind[-1])
ax3.set_ylabel('expression\nshift')
ax3.grid(b=False)
ax3.set_facecolor('white')

# Hide the right and top spines
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)

# Only show ticks on the left and bottom spines
ax3.yaxis.set_ticks_position('left')
ax3.xaxis.set_ticks_position('bottom')

# Annotate y axis
ylim = ax3.get_ylim()
#ax.set_yticks([])
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
ax3.set_yticks(yticks)
# fontProperties = {'family':'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 14}
ax3.set_yticklabels(['-','+'], fontname='Arial')
ax3.set_ylabel('expression\nshift')

# Make grid lines
xticks = [x for x in ind if x%10 == 0]
for n, x in enumerate(xticks):
        ax3.axvline(x, color='lightgray', linewidth=0.36, zorder=-1)

ax3.set_xlabel('position')

fig.subplots_adjust(hspace=-1)
plt.tight_layout()

figname_out = 'figS3_relB_sorting_expshift.pdf'
fig.savefig(output + figname_out, format='pdf')
