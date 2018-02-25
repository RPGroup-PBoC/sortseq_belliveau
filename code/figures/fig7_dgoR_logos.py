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

datadir1 = '../sortseq/20160707_purT_xylE_dgoR/'
datadir2 = '../sortseq/20160902_purTdelta_dgoRdelta/'

# Create background array
gc = .508
background_array =pd.DataFrame( [[(1-gc)/2,gc/2,gc/2,(1-gc)/2]])

#------------------------------------------------------------------------------#
# galactonate dgoR RNAP_CRP logo
#------------------------------------------------------------------------------#

# Load CRP energy model and tidy up
energy_df = pd.read_csv(datadir1 + '20160707_dgoR_MG1655_M9galactonate_na_mut3_4bins_RNAP_CRP_emat_mean.csv')
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
plt.savefig(output + 'fig7_dgoR_logo_M9galactonate_RNAP_CRP.pdf')

#------------------------------------------------------------------------------#
# glucose dgoR RNAP_CRP logo ; delta dgoR strain
#------------------------------------------------------------------------------#

# Load CRP energy model and tidy up
energy_df = pd.read_csv(datadir2 + '20160824_dgoR_MG1655deltadgoR_M9glucose_na_mut3_4bins_RNAP_CRP_emat_mean.csv')
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
    anylogo.draw(ax,effect_df=energy_df,logo_type=logo_type,use_transparency=False)
    ax.grid(False)
    ax.xlabel = 'position'
plt.tight_layout()
plt.savefig(output + 'fig7_dgoR_deltadgoR_logo_M9glucose_RNAP_CRP.pdf')
