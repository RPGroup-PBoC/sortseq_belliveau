import os
import glob
# Our numerical workhorses
import numpy as np
import pandas as pd
import scipy
# Import matplotlib stuff for plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize

from Bio import SeqIO

# Seaborn, useful for graphics
import seaborn as sns
#
# Import the project utils
import sys
sys.path.insert(0, '../../../')
import NB_sortseq_utils as utils
plt.style.use('ggplot')


# pboc_rc = { ##########SIZES##################
#                 #'lines.linewidth'       : 2,
#                 'axes.titlesize'        : 18,
#                 'axes.labelsize'        : 16,
#                 #'font.family'           : 'Lucida Sans Unicode',
#
#                 ##########COLORS#################
#                 'axes.facecolor'        :'#E3DCD0',
#
#                 #########GRIDS/TICKS############
#                 'xtick.labelsize'       : 12,
#                 'ytick.labelsize'       : 12,
#
#                 #########LEGEND#################
#                 'legend.numpoints'      : 1,
#                 'legend.fontsize'       : 13,
#                 'legend.loc'            : 'best',
#                 }

#Define the colorscheme.
# r, b, m, g, orange
# pboc  = sns.color_palette(['#d46c55', '#7aa874','#728ec1',
#                            '#aa85ab','#e08d14'])
#
#
# sns.set_context('notebook', rc=pboc_rc)
# sns.set_style('dark', rc=pboc_rc)
# sns.set_palette("deep", color_codes=True)



##########
# Fill these items in.

date = '20160710'
promoter = 'dgoR'
mutregion = 'mut3'
strain = 'MG1655'
media = 'M9glucose'
condition = 'na'
bincount = 4
TSS = 140

data_dir = '../../data/sortseq_raw/20160902_NB_004_deltas_lacRNA_pool/seq/'

# entire wild type sequence in file
Seq = "TTTTATGCATTGTTCTTTTTGTGATCTAAATTGTAGTACAACAATATAAGTTTGTACTACATTACACGCACGGCAAACGCGAACGTCATCACGCTGGTACTACAAAGTTGCCGCGTTATGCATCGATCGGGGTAAAGTAGAGAAGAACATACAGAG"

seqLength = len(Seq)

colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]
ind = np.arange(seqLength)              # the x locations for the groups

############
# generate pandas dataframe
# columns: date, strain, media, condition, bincount, position, MI, MI_error, deltabin, avg_fluor_bin1,..., WT_bp, mutation_rate, number_reads,
# initialize the DataFrame to save the mean expression levels
df = pd.read_csv('20160710_dgoR_MG1655_M9glucose_na_mut3_4bins_summary.csv')

# Mutation rate
fig1 = plt.figure(1, figsize(12.5, 3))
ax1 = plt.subplot(111, axisbg=colours[3])             # the first subplot in the first figure
plt.bar(df.sort_values(by='position').position, \
        df.sort_values(by='position').mutation_rate, alpha=1, width=0.8, \
        label="Mutation Rate")
ax1.set_xlim(df.sort_values(by='position').position[0],df.sort_values(by='position').position[len(ind)-1]+1)
ax1.set_ylabel('mutation rate')
ax1.set_xlabel('position')
plt.legend(loc="upper left");
ax1.grid(b=False)

# Hide the right and top spines
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')

plt.tight_layout()

figname_out = date + '_' + promoter + '_' + strain + '_' + media + '_' + \
                             condition + '_' + mutregion + '_' + str(bincount) + 'bins' + \
                             '_mutrate.pdf'
fig1.savefig(figname_out, format='pdf')


# footprint

fig2 = plt.figure(2,figsize(12.5, 3))
ax2 = plt.subplot(111, axisbg=colours[3])             # the first subplot in the first figure
ax2.bar(ind, \
        df.sort_values(by='position').MI, alpha=1, \
        width=0.8, color=colours[2])#, yerr=MI_prel_var) label="Raw calculation",
#ax2.set_xlim(df.sort_values(by='position').position[0],df.sort_values(by='position').position[len(ind)-1]+1)
ax2.set_xlim(ind[0],ind[-1])
ax2.set_ylabel('mutual information')
ax2.set_xlabel('position')
#plt.legend(loc="upper left");
#ax2.grid(b=False)

# # Hide the right and top spines
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
#
# # Only show ticks on the left and bottom spines
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')

plt.tight_layout()

figname_out = date + '_' + promoter + '_' + strain + '_' + media + '_' + \
                             condition + '_' + mutregion + '_' + str(bincount) + 'bins' + \
                             '_Footprint.pdf'
fig2.savefig(figname_out, format='pdf')

###########
# Delta bin

fig3 = plt.figure(3,figsize(12.5, 3))
ax3 = plt.subplot(111, axisbg=colours[3])

plt.bar(ind, df.sort_values(by='position').delta_bin, alpha=1, width=0.8, color=colours[2])
ax3.set_xlim(ind[0],ind[-1])
ax3.set_ylabel('average bin shift \n (relative to wild-type; a.u.)')
ax3.set_xlabel('position', fontsize=20)
#plt.legend(loc="upper left");
#ax3.grid(b=False)

# Hide the right and top spines
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax3.yaxis.set_ticks_position('left')
ax3.xaxis.set_ticks_position('bottom')


plt.tight_layout()

figname_out = date + '_' + promoter + '_' + strain + '_' + media + '_' + \
              condition + '_' + mutregion + '_' + str(bincount) + 'bins' + \
                             '_deltabin.pdf'
fig3.savefig(figname_out, format='pdf')


# ###########
# # Delta bin
#
# ind2 = np.arange(len(avgBin_3bpavg))
#
#
# fig4 = plt.figure(4,figsize(12.5, 3))
# ax4 = plt.subplot(111, axisbg=colours[3])
#
# plt.bar(ind2, avgBin_3bpavg, alpha=1, width=0.8, color=colours[2])
# ax4.set_xlim(ind2[0],ind2[len(ind2)-1]+1)
# ax4.set_ylabel('average bin shift \n (relative to wild-type; a.u.)')
# ax4.set_xlabel('position', fontsize=20)
# #plt.legend(loc="upper left");
# ax4.grid(b=False)
#
# # Hide the right and top spines
# ax4.spines['right'].set_visible(False)
# ax4.spines['top'].set_visible(False)
#
# # Only show ticks on the left and bottom spines
# ax4.yaxis.set_ticks_position('left')
# ax4.xaxis.set_ticks_position('bottom')
#
#
# plt.tight_layout()
#
# figname_out_3bpavg = date + '_' + promoter + '_' + strain + '_' + media + '_' + \
#               condition + '_' + mutregion + '_' + str(bincount) + 'bins' + \
#                              '_deltabin_3bpavg.pdf'
# fig4.savefig(figname_out_3bpavg, format='pdf')

# ###########
# # conditional MI
#
# def forceAspect(ax,aspect=1):
#     im = ax.get_images()
#     extent =  im[0].get_extent()
#     ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
#
#
# MI_mu_corr_cond = np.loadtxt(\
#     '20160824_MG1655_M9glucose_adenine_4bins_condMI_temp.csv', delimiter=',')
#
#
# fig4 = plt.figure(4,figsize(14, 8))
# ax4 = plt.subplot(111)
# ax4.imshow(MI_mu_corr_cond[75:,75:],cmap=cm.Spectral_r, interpolation='none',vmin=0,vmax = 0.05)# vmax=np.max(MI_mu_corr_cond)) interpolation='none', vmin=np.min(emat), vmax=np.max(emat)*0.1)
# #plt.colorbar()
# forceAspect(ax4,aspect=1)
#
# figname_out = date + '_' + strain + '_' + media + '_' + \
#                              condition + '_' + str(bincount) + 'bins' + \
#                              '_conditionalMI.pdf'
# fig4.savefig(figname_out, format='pdf')
