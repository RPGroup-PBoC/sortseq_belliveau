import os
import glob

# Our numerical tools
import numpy as np
import pandas as pd
import scipy

# Import the project utils
import sys
sys.path.insert(0, '../')
sys.path.insert(0, '../../')
import NB_sortseq_utils as utils

# import useful plotting tools
import seaborn as sns
colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize
utils.set_plotting_style1()

import ConfigParser
config = ConfigParser.RawConfigParser()

# Useful function to make directories (and pass if dir already exists)
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError:
        if os.path.isdir(path):
            pass
        else: raise

# make output folder for plots
output = 'output_figs/'
mkdir_p(output)
#------------------------------------------------------------------------------#
# load in Sort-Seq details from config file
#------------------------------------------------------------------------------#
cfg_dir = sys.argv[1]
config.read(cfg_dir)

date = config.get('Input','date')
promoter = config.get('Input','promoter')
mutregion = config.get('Input','mutregion')
strain = config.get('Input','strain')
media = config.get('Input','media')
condition = config.get('Input','condition')
bincount = config.getint('Input','bincount')

Seq = config.get('Input','seq')
data_dir = config.get('Input','data_dir')
seqLength = len(Seq)


#------------------------------------------------------------------------------#
#  Load pandas dataframe of processed summary file
#------------------------------------------------------------------------------#

fname = date + '_' + promoter + '_' + strain + '_' + media + '_' + \
        condition + '_' + mutregion + '_' + str(bincount) +\
        'bins' + '_summary.csv'
df = pd.read_csv(fname)

# ind = np.arange(seqLength) # the x locations
ind = df.sort_values(by='position').position.values
#------------------------------------------------------------------------------#
# plot mutation rate
#------------------------------------------------------------------------------#

fig1 = plt.figure(1, figsize(12.5, 3))
ax1 = plt.subplot(111, axisbg=colours[3])
plt.bar(ind, df.sort_values(by='position').mutation_rate, alpha=1, width=0.8, \
        label="Mutation Rate")
ax1.set_xlim(ind[0],ind[-1])
ax1.set_ylabel('mutation rate')
ax1.set_xlabel(r'position (relative to start of $' + promoter + '$ start codon)')
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
fig1.savefig(output + figname_out, format='pdf')


#------------------------------------------------------------------------------#
# plot information footprint
#------------------------------------------------------------------------------#

fig2 = plt.figure(2,figsize(12.5, 3))
ax2 = plt.subplot(111, axisbg=colours[3])
ax2.bar(ind, \
        df.sort_values(by='position').MI, alpha=1, \
        width=0.8, color=colours[2])
ax2.set_xlim(ind[0],ind[-1])
ax2.set_ylabel('mutual information\n(bits)')
ax2.set_xlabel(r'position (relative to start of $' + promoter + '$ start codon)')
ax2.grid(b=False)

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
fig2.savefig(output + figname_out, format='pdf')

#------------------------------------------------------------------------------#
# plot expression shifts
#------------------------------------------------------------------------------#

fig3 = plt.figure(4,figsize(12.5, 3))
ax3 = plt.subplot(111, axisbg='#FFFFFF')

plt.bar(ind, df.sort_values(by='position').expshift, alpha=1,
        width=0.8, color='k')
ax3.set_xlim(ind[0],ind[-1])
ax3.set_ylabel('expresion shift\n(relative to WT; a.u.)')
ax3.set_xlabel(r'position (relative to start of $' + promoter + '$ start codon)')
ax3.grid(b=False)

# Hide the right and top spines
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax3.yaxis.set_ticks_position('left')
ax3.xaxis.set_ticks_position('bottom')

# Make grid lines
xticks = [x for x in ind if x%20 == 0]
for n, x in enumerate(xticks):
        ax3.axvline(x, color='gray', linewidth=0.36, zorder=-1)

plt.tight_layout()

figname_out = date + '_' + promoter + '_' + strain + '_' + media + '_' + \
              condition + '_' + mutregion + '_' + str(bincount) + 'bins' + \
              '_expshift.pdf'
fig3.savefig(output + figname_out, format='pdf')

#------------------------------------------------------------------------------#
# plot expression shifts (3 bp running average)
#------------------------------------------------------------------------------#

fig4 = plt.figure(4,figsize(12.5, 3))
ax4 = plt.subplot(111, axisbg='#FFFFFF')

plt.bar(ind[:-3], df.sort_values(by='position').expshift_3bpavg[:-3],
        alpha=1, width=0.8, color='k')
ax4.set_xlim(ind[0],ind[-3])
ax4.set_ylabel('expresion shift\n(relative to WT; a.u.)')
ax4.set_xlabel(r'position (relative to start of $' + promoter + '$ start codon)')
ax4.grid(b=False)

# Hide the right and top spines
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax4.yaxis.set_ticks_position('left')
ax4.xaxis.set_ticks_position('bottom')

# Make grid lines
xticks = [x for x in ind if x%20 == 0]
for n, x in enumerate(xticks):
        ax4.axvline(x, color='gray', linewidth=0.36, zorder=-1)

plt.tight_layout()

figname_out = date + '_' + promoter + '_' + strain + '_' + media + '_' + \
              condition + '_' + mutregion + '_' + str(bincount) + 'bins' + \
              '_expshift_3bpavg.pdf'
fig4.savefig(output + figname_out, format='pdf')
