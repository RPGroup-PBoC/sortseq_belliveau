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

# load in purT files M9 glucose
data_dir_1 = '../sortseq/20160707_purT_xylE_dgoR/'
fname_1 = '20160710_purT_MG1655_M9glucose_na_mut1_4bins_summary.csv'
df_purT_1 = pd.read_csv(data_dir_1 + fname_1)

data_dir_2 = '../sortseq/20160707_purT_xylE_dgoR/'
fname_2 = '20160710_purT_MG1655_M9glucose_na_mut2_4bins_summary.csv'
df_purT_2 = pd.read_csv(data_dir_2 + fname_2)

df_purT = df_purT_1.append(df_purT_2)

# load in yebG files M9 glucose
data_dir_1 = '../sortseq/20170717_yebG/'
fname_1 = '20170717_yebG_MG1655_M9glucose_na_mut1_4bins_summary.csv'
df_yebG_1 = pd.read_csv(data_dir_1 + fname_1)

data_dir_2 = '../sortseq/20170717_yebG/'
fname_2 = '20170717_yebG_MG1655_M9glucose_na_mut2_4bins_summary.csv'
df_yebG_2 = pd.read_csv(data_dir_2 + fname_2)

df_yebG = df_yebG_1.append(df_yebG_2)


#------------------------------------------------------------------------------#
# make shared axis plots of information footprint and delta bin shift
#------------------------------------------------------------------------------#

seqLength = 105
colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]

fig1 = plt.figure(figsize=utils.cm2inch((5,3.5)))
ax1 = fig1.add_subplot(111)
# fig= plt.figure(figsize=cm2inch((8.4,4.6)))
# gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1],
#          wspace=0.0, hspace=0, top=0.95, bottom=0.5, left=0.17, right=0.845)
# ax1= plt.subplot(gs[0,0])

# make array for x positions - note that this will be relative to start codon
ind_purT = np.arange(seqLength) - 120 #

y_purT = df_purT[df_purT['MI']!=0].groupby(['position'])['expshift'].mean().values
#purT expression shift
ax1.bar(ind_purT[55:], y_purT[55:], width=0.8, \
        linewidth = 0, color = 'k')
#ax1.set_xlim(ind_purT[0],ind_purT[-1])
ax1.set_ylabel('expression\nshift', fontsize=8)
ax1.set_xlabel('position\n(relative to $purT$ gene)', fontsize=8)
ax1.grid(b=False)
# ax1.set_axis_bgcolor('white')
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

ax1.set_xlim([ind_purT[55]-.5, ind_purT[-1]+.5])

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
figname_out = 'fig5_purTyebG_expshift_purT.pdf'
fig1.savefig(output + figname_out, format='pdf')

fig2 = plt.figure(figsize=utils.cm2inch((5,3.5)))
ax2 = fig2.add_subplot(111)
# information footprint yebG
# try different indexing approach
#ind_yebG = df_yebG[df_yebG['MI']!=0].groupby(['position'])['position'].mean().values

# ax2= plt.subplot(gs[1,0])
ind_yebG = np.arange(seqLength) - 120 # assumed to be at position 6 relative to -10, so approx
#59
y_yebG = df_yebG[df_yebG['MI']!=0].groupby(['position'])['expshift'].mean().values
ax2.bar(ind_yebG[50:], y_yebG[50:], width=0.8, \
        linewidth = 0, color = 'k')
#ax2.bar(df_yebG.groupby(['position'])['position'].mean(),df_yebG.groupby(['position'])['MI'].mean(), width=0.8,linewidth = 0, color = 'k')


#ax2.set_xlim(ind_yebG[0],ind_yebG[-1])
#ax2.set_xlim(ind_yebG[0],ind_yebG[-1])
ax2.set_ylabel('expression\nshift', fontsize=8)
ax2.set_xlabel('position\n(relative to $yebG$ gene)', fontsize=8)
ax2.grid(b=False)

# ax2.set_axis_bgcolor('white')
ax2.set_facecolor('white')

# Hide the right and top spines
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(True)
ax2.spines['left'].set_visible(False)
# ax2.spines['bottom'].set_linewidth(0.36)

# Only show ticks on the left and bottom spines
ax2.xaxis.set_ticks_position('bottom')
ax2.tick_params(axis='y', left = 'off')
# ax2.xaxis.set_tick_params(width=0.36)

# Adjust xlim to accomodate entire bars on ends
ax2.set_xlim([ind_yebG[55]-.5, ind_yebG[-1]+.5])
ax2.invert_xaxis()
# Annotate y axis
ylim = ax2.get_ylim()
#ax.set_yticks([])
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
ax2.set_yticks(yticks)
# fontProperties = {'family':'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 14}
ax2.set_yticklabels(['-','+'], fontname='Arial')
ax2.set_ylabel('expression\nshift')

# Make grid lines
xticks = [x for x in ind_yebG if x%20 == 0]
for n, x in enumerate(xticks):
        ax2.axvline(x, color='lightgray', linewidth=0.36, zorder=-1)

plt.tight_layout()
figname_out = 'fig5_purTyebG_expshift_yebG.pdf'
fig2.savefig(output + figname_out, format='pdf')
