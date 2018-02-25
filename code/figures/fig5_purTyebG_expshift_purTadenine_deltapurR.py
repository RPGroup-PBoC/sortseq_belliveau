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
utils.set_plotting_style4()

#===============================================================================
# Set output directory
#===============================================================================
output = 'output_figs/'

#------------------------------------------------------------------------------#
# Load in the summary csv file
#------------------------------------------------------------------------------#

# load in purT files M9 glucose + adenine in delta purR strain.
data_dir = '../sortseq/20160902_purTdelta_dgoRdelta/'
fname = '20160824_purT_MG1655deltapurR_M9glucose_adenine_mut1_4bins_summary.csv'
df_purT = pd.read_csv(data_dir + fname)

#------------------------------------------------------------------------------#
# make expression shift plots
#------------------------------------------------------------------------------#

seqLength = 60
colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]

fig1 = plt.figure(figsize=utils.cm2inch((5,3.5)))
ax1 = fig1.add_subplot(111)
# fig= plt.figure(figsize=cm2inch((8.4,4.6)))
# gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1],
#          wspace=0.0, hspace=0, top=0.95, bottom=0.5, left=0.17, right=0.845)
# ax1= plt.subplot(gs[0,0])

# make array for x positions - note that this will be relative to start codon
ind_purT = np.arange(seqLength) - 75 #

y_purT = df_purT[df_purT['MI']!=0].groupby(['position'])['expshift'].mean().values
#purT expression shift
ax1.bar(ind_purT, y_purT, width=0.8, \
        linewidth = 0, color = 'k')
#ax1.set_xlim(ind_purT[0],ind_purT[-1])
ax1.set_ylabel('expression\nshift', fontsize=8)
ax1.set_xlabel('position\n(relative to $purT$ gene)', fontsize=8)
ax1.grid(b=False)
ax1.set_facecolor('white')

# Hide the right and top spines
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(True)
ax1.spines['left'].set_visible(False)
# ax1.spines['bottom'].set_linewidth(0.36)

# Only show ticks on the left and bottom spines
# ax1.xaxis.set_ticks_position('bottom')
ax1.tick_params(axis='y', left = 'off')
# ax1.xaxis.set_tick_params(width=0.36)

#plt.setp(ax1.spines.values(),color ='k')
# for spine in ax1.spines.values():
#     spine.set_edgecolor('k')
    #plt.setp(ax1.get_xticklabels(), visible=False)
#plt.setp(ax1.get_xticklabels(), fontsize=6)

#plt.setp(ax1.get_yticklabels(), visible=False)

# Adjust xlim to accomodate entire bars on ends

ax1.set_xlim([ind_purT[5]-.5, ind_purT[-1]+.5])

# Annotate y axis
ylim = ax1.get_ylim()
#ax.set_yticks([])
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
ax1.set_yticks(yticks)
# fontProperties = {'family':'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 14}
ax1.set_yticklabels(['-','+'], fontname='Arial')
ax1.set_ylabel('expression\nshift')

# Make grid lines
xticks = [x for x in ind_purT if x%20 == 0]
for n, x in enumerate(xticks):
        ax1.axvline(x, color='lightgray', linewidth=0.36, zorder=-1)

plt.tight_layout()
figname_out = 'fig5_purTyebG_expshift_purTadenine_deltapurR.pdf'
fig1.savefig(output + figname_out, format='pdf')
