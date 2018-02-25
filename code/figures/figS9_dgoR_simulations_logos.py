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

# datadir = '../sortseq/'
# # lac
# datadir_lac = '../sortseq/2011_Djones_lacZ/'
# # rel
# datadir_rel = '../sortseq/20150312_rel_dataanalysis/'
# # mar
# datadir_marA = '../sortseq/20150820_marMut2_marrelsublib/'
# # pur
# datadir_pur = '../sortseq/20160707_NB_003_pD_X_P/'


#------------------------------------------------------------------------------#
# function for PWM
#------------------------------------------------------------------------------#


# Function that can get the counts matrix for any TF
def get_counts_matrix(tf_name, db_file, padding=2):

    # Read in database
    column_indices = [1,11]
    column_names = ['name','site']
    dtype_dict = {'name':str, 'site':str}
    tf_df = pd.read_csv(db_file,sep='\t',comment='#',usecols=column_indices,names=column_names,dtype=dtype_dict)

    # Make sure TF is in database
    tf_names = set(tf_df['name'])
    assert tf_name in tf_names, 'Error! No TF named %s in database.'%tf_name

    # Extract rows for TF of interest
    df = tf_df[tf_df['name']==tf_name].dropna(axis=0)

    # Extract uppercase site with 2 lowercase surrounding sites
    pattern = re.compile('[acgt]*([acgt]{%d}[ACGT]*[acgt]{%d})[acgt]*'%(padding,padding))
    df['site'] = [pattern.match(site).group(1).upper() for site in df['site']]

    # Compute lengths for each site
    df['length'] = [len(site) for site in df['site']]

    # Only take rows with sites of the right size
    length_mode = df['length'].mode()[0]
    df = df[df['length']==length_mode]

    # Fill counts matrix
    counts_matrix = np.zeros([length_mode,4])
    bases = 'ACGT'
    for s in df['site']:
        for i in range(length_mode):
            for b, base in enumerate(bases):
                counts_matrix[i,b] += (s[i] == base)

    # Fill counts_df
    counts_df = pd.DataFrame(data=counts_matrix,columns=list(bases))
    counts_df.index.name='pos'

    return counts_df



def get_counts_seq(db_file):

    # Read in database
    column_indices = [1]
    column_names = ['site']
    dtype_dict = {'site':str}
    df = pd.read_csv(db_file,comment='#',names=column_names,dtype=dtype_dict)

    # Compute lengths for each site
    df['length'] = [len(site) for site in df['site']]

    # Only take rows with sites of the right size
    length_mode = df['length'].mode()[0]
    df = df[df['length']==length_mode]

    # Fill counts matrix
    counts_matrix = np.zeros([length_mode,4])
    bases = 'ACGT'
    for s in df['site']:
        for i in range(length_mode):
            for b, base in enumerate(bases):
                counts_matrix[i,b] += (s[i] == base)

    # Fill counts_df
    counts_df = pd.DataFrame(data=counts_matrix,columns=list(bases))
    counts_df.index.name='pos'

    return counts_df

# Create background array
gc = .508
background_array =pd.DataFrame( [[(1-gc)/2,gc/2,gc/2,(1-gc)/2]])

#------------------------------------------------------------------------------#
# find RNAP emats
#------------------------------------------------------------------------------#

# # if I want to plot all emats in simulations folder:
# root = 'input_data/simulations/'
# filename = '*emat_mean.csv'
# fileslist = []
# for directory,subdirs,files in os.walk(root):
#     for i,j in enumerate(files):
#         if '_emat_mean.csv' in j:
#             fileslist.append(directory+'/'+j)

fileslist = ['../sortseq/20160902_purTdelta_dgoRdelta/20160824_dgoR_MG1655deltadgoR_M9glucose_na_mut3_4bins_RNAP_CRP_emat_mean.csv',
'input_data/simulations/20170822_dgoR_sim_1RNAP_relB/pymc_pdgoR_sim_RNAP_relB/pymc_pdgoR_sim_RNAP_relB_emat_mean.csv',
'input_data/simulations/20170822_dgoR_sim_2RNAP_relB/20170822_pdgor_2RNAP_relB_shifted_pymc_scalefactor1_all_MCMC_100_0/20170822_pdgor_2RNAP_relB_shifted_pymc_scalefactor1_all_MCMC_100_0_emat_mean.csv',
'input_data/simulations/20170822_dgoR_sim_3RNAP_relB/pymc_pdgoR_sim_RNAP_relB/pymc_pdgoR_sim_RNAP_relB_emat_mean.csv']

print(fileslist)

# determine number of plots in figure
rnap_count=0
for i, fname in enumerate(fileslist):
    if 'delta' in fname:
        continue
    if 'binning' in fname:
        continue
    rnap_count += 1

rnap_count = len(fileslist)
count = 0
fig, ax = plt.subplots(rnap_count,1,figsize=[8,rnap_count*2])
for i, fname in enumerate(fileslist):
    energy_df = pd.read_csv(fname)
    if '2RNAP' in fname:
        energy_df = energy_df[energy_df.position >= 5].reset_index()
    # if '1RNAP' in fname:
    #     energy_df_temp = energy_df[ energy_df.index<=17]*0
    #
    #     energy_df_temp[['A','C','G','T']] =10**-10
    #
    #     energy_df.index = energy_df.index + 18
    #     energy_df = energy_df_temp.append(energy_df)
    #     energy_df = pd.DataFrame(energy_df[['A','C','G','T']].values, columns=['A','C','G','T'],dtype='float64')


    # energy_df = energy_df[energy_df.position <= 31]


    energy_df = energy_df[['A','C','G','T']]
    energy_df = utils.estimate_scalefactor(np.array(energy_df))*energy_df

    # create background nucleotide frequencies dataframe
    background_df = pd.DataFrame(pd.np.tile(background_array,
                        (len(energy_df), 1)), columns=['A','C','G','T'])


    # Draw logos using energy matrix
    logo_types = ['information']
    num_types = len(logo_types)

    anylogo.draw(ax[i],effect_df=energy_df,logo_type='information',
                use_transparency=False,
                background = background_df,
                ylabel='information\n(bits)')
    ax[i].grid(False)
    ax[i].set_facecolor('white')
    plt.setp(ax[i].get_xticklabels(), visible=False)

    #
    # anylogo.draw(ax2,effect_df=-pwm_df,logo_type='information',use_transparency=False,
    #             ylabel='information\n(bits)')
    # ax2.grid(False)
    # ax2.set_facecolor('white')
    # plt.setp(ax2.get_xticklabels(), visible=False)
    count +=1
plt.tight_layout()
fig.savefig(output + 'figS9_dgoR_simulations_compare.pdf')
