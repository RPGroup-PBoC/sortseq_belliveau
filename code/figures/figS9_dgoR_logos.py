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
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize

# Logo-generating module
import anylogo

utils.set_plotting_style1()
#===============================================================================
# Set output directory based on the graphicspath.tex file to print in dropbox
#===============================================================================
output = 'output_figs/'
#------------------------------------------------------------------------------#
# directory where emat csv files are contained
#------------------------------------------------------------------------------#

datadir1 = '../sortseq/20160920_dgoRJK10_xylEJK10/'

# Create background array
gc = .508
background_array =pd.DataFrame( [[(1-gc)/2,gc/2,gc/2,(1-gc)/2]])

#------------------------------------------------------------------------------#
# Ezrich media JK10 RNAP_CRP logo; 500 uM cAMP galactonate
#------------------------------------------------------------------------------#

# Load CRP energy model and tidy up
energy_df = pd.read_csv(datadir1 + '20160921_dgoR_JK10_EZrichgalactonate_500cAMP_mut3_4bins_RNAP_emat_mean.csv')
energy_df = energy_df[['A','C','G','T']]
# energy_df.columns = [col.split('_')[-1] for col in energy_df.columns]

energy_df = utils.estimate_scalefactor(np.array(energy_df))*energy_df

# create background nucleotide frequencies dataframe
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df), 1)), columns=['A','C','G','T'])

# Draw logos using energy matrix
logo_types = ['information']
num_types = len(logo_types)

fig = plt.figure(figsize=[8,2*num_types])
for n, logo_type in enumerate(logo_types):

    ax = fig.add_subplot(num_types,1,n+1)
    anylogo.draw(ax,effect_df=energy_df,logo_type=logo_type,
                use_transparency=False,
                background = background_df)
    ax.grid(False)
    ax.xlabel = 'position'
plt.tight_layout()
plt.savefig(output + 'figS9_dgoR_JK10_logos_500cAMP.pdf')

#------------------------------------------------------------------------------#
# Ezrich media JK10 RNAP_CRP logo; 0 uM cAMP glucose
#------------------------------------------------------------------------------#

# Load CRP energy model and tidy up
energy_df = pd.read_csv(datadir1 + '20160921_dgoR_JK10_EZrichglucose_0cAMP_mut3_4bins_RNAP_emat_mean.csv')
energy_df = energy_df[['A','C','G','T']]
# energy_df.columns = [col.split('_')[-1] for col in energy_df.columns]

energy_df = utils.estimate_scalefactor(np.array(energy_df))*energy_df

# Draw logos using energy matrix
logo_types = ['information']
num_types = len(logo_types)

fig = plt.figure(figsize=[8,2*num_types])
for n, logo_type in enumerate(logo_types):

    ax = fig.add_subplot(num_types,1,n+1)
    anylogo.draw(ax,effect_df=energy_df,logo_type=logo_type,
                use_transparency=False,
                background = background_df)
    ax.grid(False)
    ax.xlabel = 'position'
plt.tight_layout()
plt.savefig(output + 'figS9_dgoR_JK10_logos_0cAMP.pdf')
