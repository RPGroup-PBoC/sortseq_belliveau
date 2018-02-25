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
# Set output directory based on the graphicspath.tex file to print in dropbox
#===============================================================================
output = 'output_figs/'

#------------------------------------------------------------------------------#
# Load in the summary csv file
#------------------------------------------------------------------------------#

data_dir = '../sortseq/20150312_relB/'
fname = '20150312_relB_MG1655_M9glucose_na_mut1_4bins_summary.csv'
df = pd.read_csv(data_dir + fname)

seqLength = len(df.position.unique())
#------------------------------------------------------------------------------#
# make shared axis plots of information footprint and delta bin shift
#------------------------------------------------------------------------------#
colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]

fig = plt.figure(figsize=utils.cm2inch((9,3)))
ax1 = fig.add_subplot(111)
# make array for x positions - note that this will be relative to TSS
ind = np.arange(seqLength) - 35

# delta bin
ax1.bar(ind, df.sort_values(by='position').expshift, width=0.8, \
        linewidth = 0, color = 'k')
ax1.set_xlim(ind[0]-0.5,ind[0]+132+0.5)
ax1.set_xlabel('position', fontsize=8)
ax1.grid(b=False)
ax1.set_facecolor('white')

# Hide the right and top spines
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
# ax1.set_xlim([ind[0]-.5, ind[-1]+.5])

# Annotate y axis
ylim = ax1.get_ylim()
#ax.set_yticks([])
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
ax1.set_yticks(yticks)
# fontProperties = {'family':'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 14}
ax1.set_yticklabels(['-','+'], fontname='Arial')
ax1.set_ylabel('expression\nshift')

# Make grid lines
xticks = [x for x in ind if x%20 == 0]
for n, x in enumerate(xticks):
        ax1.axvline(x, color='lightgray', linewidth=0.36, zorder=-1)


plt.tight_layout()

figname_out = 'fig2_relB_main_expshift.pdf'
fig.savefig(output + figname_out, format='pdf')
