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
from matplotlib.ticker import FormatStrFormatter

# Seaborn, useful for graphics
import seaborn as sns

sns.set_palette("deep", color_codes=True)
utils.set_plotting_style1()

#===============================================================================
# Set output directory
#===============================================================================
output = 'output_figs/'

#------------------------------------------------------------------------------#
# Load in the summary csv file
#------------------------------------------------------------------------------#

df = pd.read_csv('RegulonDB_20170630_BindingSiteSet.csv',comment='#')




df['len'] = df['TF_sequence'].str.len()
df['len_half'] = df['TF_sequence'].str.len()/2
df['right_side_bs'] = df['center_position_bs'] + df['len_half']

hist = np.zeros(1)
for index, row in df.iterrows():
    # print(row['right_side_bs']-row['len'])
    if pd.isnull(row['right_side_bs']):
        continue
    hist = np.append(hist,np.arange(row['right_side_bs']-row['len'],row['right_side_bs']), axis=0)

data = hist
# data = df.right_side_bs.dropna()

fig, ax = plt.subplots()
counts, bins, patches = plt.hist(data,
                                 bins=np.arange(min(data), max(data) + 1, 1),
                                 linewidth = 0)#, histtype = 'stepfilled')

# calculate the fraction of binding sites positioned within a specific window
# We want to capture the region downstream of the start site, however, lets
# assume we will go as far as the approximate rbs region since doing
# Sort-Seq probabily will probably not be useful there.
window_150 = 0.0
for i in data:
    if -134 < i < 16:
        window_150 += 1.0
window_150_frac = window_150/len(data)
# Change the colors of bars at the edges...
twentyfifth, seventyfifth = np.percentile(data, [25, 75])
for patch, rightside, leftside in zip(patches, bins[1:], bins[:-1]):
    if leftside > 16:
        patch.set_facecolor('grey')
    if rightside < -134:
        patch.set_facecolor('grey')
from matplotlib.patches import Rectangle

#create legend
handles = [Rectangle((0,0),1,1,ec="none")]
labels= ["150 bp window ("+ str(window_150_frac*100)[:4] + "% \n known TF binding sites)"]
plt.legend(handles, labels)

ax.set_xlim(-400,300)
ax.set_xlabel('position relative to TSS')
ax.set_ylabel('number of binding sites \n overlapping at a base pair')

plt.tight_layout()
figname_out = 'figS6_RegulonDB_bindingsites_summary.pdf'
fig.savefig(output + figname_out, format='pdf')
