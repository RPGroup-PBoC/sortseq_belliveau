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

# lac
datadir_lac = '../sortseq/2011_lacZ/'
# rel
datadir_rel = '../sortseq/20150312_relB/'
# mar
datadir_marA = '../sortseq/20150820_marRmut2only/'
# pur
datadir_pur = '../sortseq/20160707_purT_xylE_dgoR/'
# lexA
datadir_lexA = '../sortseq/20170717_yebG/'

# fis; mar mut2 MG1655
datadir_MG1655 = '../sortseq/20150820_marRmut2only/'

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

# def mean2(x):
#     y = np.sum(x) / np.size(x);
#     return y
#
# def corr2(a,b):
#     a = a - mean2(a)
#     b = b - mean2(b)
#
#     r = (a*b).sum() / np.sqrt((a*a).sum() * (b*b).sum());
#     return r

def corr(a,b):

    # a = a - a.mean()
    # b = b - b.mean()
    r = (a*b).sum()# / np.sqrt((a*a).sum() * (b*b).sum());
    return r



#------------------------------------------------------------------------------#
# logos for LacI
#------------------------------------------------------------------------------#

# Get counts matrix and turn into pwm
counts_df = get_counts_seq('input_data/regulondb_PSSMSet_MEME_LacI.txt')

# Get number of CRP sites
num_sites = counts_df.iloc[0,:].sum()
print('num_sites = %d'%num_sites)

# Compute prob matrix
prob_df = (counts_df+0.025)/(num_sites+0.1)

# Create background dict
gc = .508
background = {'A':(1-gc)/2,'C':gc/2,'G':gc/2,'T':(1-gc)/2}

# Compute pwm
background_df = prob_df.copy()
for i in range(len(prob_df)):
    background_df.loc[i,:] = background

pwm_df = np.log(prob_df/background_df)
pwm_df = pwm_df[['A','C','G','T']]

# print(pwm_df)
pwm_df_fix = -utils.fix_matrix_gauge(np.array(pwm_df.copy()).T)
# print(pwm_df_fix)
# print(pwm_df)
# pwm_df_fix2 = -utils.fix_matrix(np.array(pwm_df.copy()))
# print(pwm_df_fix2)


#-----

energy_df = pd.read_csv(datadir_lac + '2011_lacZ_MG1655_M9glucose_na_mut1_4bins_LacI_O1_emat_mean.csv')
energy_df = energy_df[['A','C','G','T']]

# energy_df_fix = utils.fix_matrix(np.array(energy_df.copy()))
# energy_df_fix = utils.fix_matrix_gauge(np.array(energy_df.copy()).T)
# print(energy_df)
energy_df_fix = utils.fix_matrix_gauge(np.array(energy_df.copy()).T,)
# print(energy_df_fix)

# print(energy_df_fix2)

energy_df = utils.estimate_scalefactor(np.array(energy_df))*energy_df


# Draw logos using energy matrix
logo_types = ['information']
num_types = len(logo_types)

fig1, (ax1, ax2) = plt.subplots(2,1,figsize=[8,4])

anylogo.draw(ax1,effect_df=energy_df,logo_type='information',
            use_transparency=False,
            background = background_df,
            ylabel='information\n(bits)')
ax1.grid(False)
ax1.set_facecolor('white')
plt.setp(ax1.get_xticklabels(), visible=False)

anylogo.draw(ax2,effect_df=-pwm_df,logo_type='information',
            use_transparency=False,
            ylabel='information\n(bits)')
ax2.grid(False)
ax2.set_facecolor('white')
plt.setp(ax2.get_xticklabels(), visible=False)

plt.tight_layout()
fig1.savefig(output + 'figS4_lacI_logo_compare.pdf')

correlation_LacI = corr(pwm_df_fix,energy_df_fix)

#------------------------------------------------------------------------------#
# logos for CRP
#------------------------------------------------------------------------------#

# Get counts matrix and turn into pwm
counts_df = get_counts_seq('input_data/regulondb_PSSMSet_MEME_CRP.txt')

# Get number of CRP sites
num_sites = counts_df.iloc[0,:].sum()
print('num_sites = %d'%num_sites)

# Compute prob matrix
prob_df = (counts_df+0.1)/(num_sites+0.4)

# Compute pwm
background_df = prob_df.copy()
for i in range(len(prob_df)):
    background_df.loc[i,:] = background
pwm_df = np.log(prob_df/background_df)

pwm_df_fix = -utils.fix_matrix_gauge(np.array(pwm_df.copy()).T)
#-----

energy_df = pd.read_csv(datadir_lac + '2011_lacZ_MG1655_M9glucose_na_mut1_4bins_CRP_emat_mean.csv')
energy_df = energy_df[energy_df.position != energy_df.position.min()]
energy_df = energy_df[energy_df.position != energy_df.position.min()]
energy_df = energy_df[energy_df.position != energy_df.position.min()]
energy_df = energy_df[energy_df.position != energy_df.position.min()].reset_index()
energy_df = energy_df[['A','C','G','T']]

energy_df_fix = utils.fix_matrix_gauge(np.array(energy_df.copy()).T)

energy_df = utils.estimate_scalefactor(np.array(energy_df))*energy_df

fig2, (ax1, ax2) = plt.subplots(2,1,figsize=[8,4])

anylogo.draw(ax1,effect_df=energy_df,logo_type='information',use_transparency=False, background = background_df,
            ylabel='information\n(bits)')
ax1.grid(False)
ax1.set_facecolor('white')
plt.setp(ax1.get_xticklabels(), visible=False)

anylogo.draw(ax2,effect_df=-pwm_df,logo_type='information',use_transparency=False,
            ylabel='information\n(bits)')
ax2.grid(False)
ax2.set_facecolor('white')
plt.setp(ax2.get_xticklabels(), visible=False)


plt.tight_layout()
fig2.savefig(output + 'figS4_CRP_logo_compare.pdf')

correlation_CRP = corr(pwm_df_fix,energy_df_fix)

#------------------------------------------------------------------------------#
# logos for MarA
#------------------------------------------------------------------------------#

# Get counts matrix and turn into pwm
counts_df = get_counts_seq('input_data/regulondb_PSSMSet_MEME_MarA.txt')

# Get number of CRP sites
num_sites = counts_df.iloc[0,:].sum()
print('num_sites = %d'%num_sites)

# Compute prob matrix
prob_df = (counts_df+0.1)/(num_sites+0.4)

# Compute pwm
background_df = prob_df.copy()
for i in range(len(prob_df)):
    background_df.loc[i,:] = background
pwm_df = np.log(prob_df/background_df)

pwm_df_fix = -utils.fix_matrix_gauge(np.array(pwm_df.copy()).T)
#-----

energy_df = pd.read_csv(datadir_marA + '20150820_marR_MG1655_LB_na_mut2_4bins_MarA_emat_mean.csv')
energy_df = energy_df[energy_df.position != energy_df.position.max()]

energy_df = energy_df[['A','C','G','T']]

energy_df_fix = utils.fix_matrix_gauge(np.array(energy_df.copy()).T)

energy_df = utils.estimate_scalefactor(np.array(energy_df))*energy_df


fig3, (ax1, ax2) = plt.subplots(2,1,figsize=[8,4])

anylogo.draw(ax1,effect_df=energy_df,logo_type='information',use_transparency=False, background = background_df,
            ylabel='information\n(bits)')
ax1.grid(False)
ax1.set_facecolor('white')
plt.setp(ax1.get_xticklabels(), visible=False)

anylogo.draw(ax2,effect_df=-pwm_df,logo_type='information',use_transparency=False,
            ylabel='information\n(bits)')
ax2.grid(False)
ax2.set_facecolor('white')
plt.setp(ax2.get_xticklabels(), visible=False)

plt.tight_layout()
fig3.savefig(output + 'figS4_marA_logo_compare.pdf')

correlation_MarA = corr(pwm_df_fix,energy_df_fix)

#------------------------------------------------------------------------------#
# logos for PurR
#------------------------------------------------------------------------------#

# Get counts matrix and turn into pwm
counts_df = get_counts_seq('input_data/regulondb_PSSMSet_MEME_PurR.txt')

# Get number of CRP sites
num_sites = counts_df.iloc[0,:].sum()
print('num_sites = %d'%num_sites)

# Compute prob matrix
prob_df = (counts_df+0.1)/(num_sites+0.4)

# Compute pwm
background_df = prob_df.copy()
for i in range(len(prob_df)):
    background_df.loc[i,:] = background
pwm_df = np.log(prob_df/background_df)

pwm_df_fix = -utils.fix_matrix_gauge(np.array(pwm_df.copy()).T)
#-----

energy_df = pd.read_csv(datadir_pur + '20160710_purT_MG1655_M9glucose_adenine_mut1_4bins_PurR_emat_mean.csv')
energy_df = energy_df[energy_df.position != energy_df.position.max()]

energy_df = energy_df[['A','C','G','T']]

energy_df_fix = utils.fix_matrix_gauge(np.array(energy_df.copy()).T)

energy_df = utils.estimate_scalefactor(np.array(energy_df))*energy_df


fig4, (ax1, ax2) = plt.subplots(2,1,figsize=[8,4])

anylogo.draw(ax1,effect_df=energy_df,logo_type='information',use_transparency=False,  background = background_df,
            ylabel='information\n(bits)')
ax1.grid(False)
ax1.set_facecolor('white')
plt.setp(ax1.get_xticklabels(), visible=False)

anylogo.draw(ax2,effect_df=-pwm_df,logo_type='information',use_transparency=False,
            ylabel='information\n(bits)')
ax2.grid(False)
ax2.set_facecolor('white')
plt.setp(ax2.get_xticklabels(), visible=False)

plt.tight_layout()
fig4.savefig(output + 'figS4_purR_logo_compare.pdf')

correlation_PurR = corr(pwm_df_fix,energy_df_fix)

#------------------------------------------------------------------------------#
# logos for xylR
#------------------------------------------------------------------------------#

# Get counts matrix and turn into pwm
counts_df = get_counts_seq('input_data/regulondb_PSSMSet_MEME_XylR_trimmed.txt')

# Get number of CRP sites
num_sites = counts_df.iloc[0,:].sum()
print('num_sites = %d'%num_sites)

# Compute prob matrix
prob_df = (counts_df+0.1)/(num_sites+0.4)

# Compute pwm
background_df = prob_df.copy()
for i in range(len(prob_df)):
    background_df.loc[i,:] = background
pwm_df = np.log(prob_df/background_df)

pwm_df_fix = -utils.fix_matrix_gauge(np.array(pwm_df.copy()).T)
#-----

# energy_df_left = pd.read_csv(datadir_pur + 'xylrxylose_0_emat_mean_temp_left.csv')
energy_df_left = pd.read_csv(datadir_pur + '20160710_xylE_MG1655_M9xylose_na_mut2_4bins_XylR_emat_mean_left.csv')
energy_df_left = energy_df_left[['A','C','G','T']]

energy_df_fix_left = utils.fix_matrix_gauge(np.array(energy_df_left.copy()).T)

energy_df_left = utils.estimate_scalefactor(np.array(energy_df_left))*energy_df_left

# energy_df_right = pd.read_csv(datadir_pur + 'xylrxylose_0_emat_mean_temp_right.csv')
energy_df_right = pd.read_csv(datadir_pur + '20160710_xylE_MG1655_M9xylose_na_mut2_4bins_XylR_emat_mean_right.csv')
energy_df_right = energy_df_right[['A','C','G','T']]

energy_df_fix_right = utils.fix_matrix_gauge(np.array(energy_df_right.copy()).T)

energy_df_right = utils.estimate_scalefactor(np.array(energy_df_right))*energy_df_right

fig5, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=[8,6])

anylogo.draw(ax1,effect_df=energy_df_left,logo_type='information',use_transparency=False, background = background_df,
            ylabel='information\n(bits)')
ax1.grid(False)
ax1.set_facecolor('white')
plt.setp(ax1.get_xticklabels(), visible=False)

anylogo.draw(ax2,effect_df=energy_df_right,logo_type='information',use_transparency=False, background = background_df,
            ylabel='information\n(bits)')
ax2.grid(False)
ax2.set_facecolor('white')
plt.setp(ax2.get_xticklabels(), visible=False)

anylogo.draw(ax3,effect_df=-pwm_df,logo_type='information',use_transparency=False,
            ylabel='information\n(bits)')
ax3.grid(False)
ax3.set_facecolor('white')
plt.setp(ax3.get_xticklabels(), visible=False)

plt.tight_layout()
fig5.savefig(output + 'figS4_xylR_logo_compare.pdf')

correlation_XylR_left = corr(pwm_df_fix,energy_df_fix_left)

correlation_XylR_right = corr(pwm_df_fix,energy_df_fix_right)

correlation_XylR_compare = corr(energy_df_fix_left,energy_df_fix_right)



#------------------------------------------------------------------------------#
# logos for LexA
#------------------------------------------------------------------------------#

# Get counts matrix and turn into pwm
counts_df = get_counts_seq('input_data/regulondb_PSSMSet_MEME_LexA_trimmed.txt')

# Get number of CRP sites
num_sites = counts_df.iloc[0,:].sum()
print('num_sites = %d'%num_sites)

# Compute prob matrix
prob_df = (counts_df+0.1)/(num_sites+0.4)

# Compute pwm
background_df = prob_df.copy()
for i in range(len(prob_df)):
    background_df.loc[i,:] = background
pwm_df = np.log(prob_df/background_df)

# for correlation calc. ignore bp that overlap -10 RNAP site
pwm_df_fix = pwm_df.copy()
pwm_df_fix = pwm_df_fix[pwm_df_fix.index >= 4]
pwm_df_fix = -utils.fix_matrix_gauge(np.array(pwm_df_fix.copy()).T)

#-----

energy_df = pd.read_csv(datadir_lexA + '20170717_yebG_MG1655_M9glucose_na_mut2_4bins_LexA_emat_mean.csv')
energy_df = energy_df[energy_df.position != energy_df.position.max()]
energy_df_fix = energy_df.copy()

energy_df = energy_df[['A','C','G','T']]

# re-fix gauge - for correlation calc. ignore bp that overlap -10 RNAP site

energy_df_fix = energy_df_fix[energy_df_fix.index >= 4]
energy_df_fix = energy_df_fix[['A','C','G','T']]
energy_df_fix = utils.fix_matrix_gauge(np.array(energy_df_fix.copy()).T)

energy_df = utils.estimate_scalefactor(np.array(energy_df))*energy_df



fig6, (ax1, ax2) = plt.subplots(2,1,figsize=[8,4])

anylogo.draw(ax1,effect_df=energy_df,logo_type='information',use_transparency=False, background = background_df,
            ylabel='information\n(bits)')
ax1.grid(False)
ax1.set_facecolor('white')
plt.setp(ax1.get_xticklabels(), visible=False)

anylogo.draw(ax2,effect_df=-pwm_df,logo_type='information',use_transparency=False,
            ylabel='information\n(bits)')
ax2.grid(False)
ax2.set_facecolor('white')
plt.setp(ax2.get_xticklabels(), visible=False)

plt.tight_layout()
fig6.savefig(output + 'figS4_lexA_logo_compare.pdf')

# calc. correlation
correlation_LexA = corr(pwm_df_fix,energy_df_fix)

#------------------------------------------------------------------------------#
# logos for fis in MG1655 strain
#------------------------------------------------------------------------------#

# Get counts matrix and turn into pwm
counts_df = get_counts_seq('input_data/regulondb_PSSMSet_MEME_fis.txt')

# Get number of CRP sites
num_sites = counts_df.iloc[0,:].sum()
print('num_sites = %d'%num_sites)

# Compute prob matrix
prob_df = (counts_df+0.1)/(num_sites+0.4)

# Create background dict
gc = .508
background = {'A':(1-gc)/2,'C':gc/2,'G':gc/2,'T':(1-gc)/2}

# Compute pwm
background_df = prob_df.copy()
for i in range(len(prob_df)):
    background_df.loc[i,:] = background
pwm_df = np.log(prob_df/background_df)

pwm_df_fix = -utils.fix_matrix_gauge(np.array(pwm_df.copy()).T)
#-----

energy_df = pd.read_csv(datadir_MG1655 + '20150820_marR_MG1655_LB_na_mut2_4bins_Fis_emat_mean.csv')
# energy_df = energy_df[energy_df.position != energy_df.position.max()]
# energy_df = energy_df[energy_df.position != energy_df.position.max()]
energy_df = energy_df[energy_df.index >= 2]
energy_df = energy_df[energy_df.index <= (energy_df.position.max()-2)].reset_index()

energy_df = energy_df[['A','C','G','T']]

energy_df_fix = utils.fix_matrix_gauge(np.array(energy_df.copy()).T)

energy_df = utils.estimate_scalefactor(np.array(energy_df))*energy_df


fig6, (ax1, ax2) = plt.subplots(2,1,figsize=[8,4])

anylogo.draw(ax1,effect_df=energy_df,logo_type='information',use_transparency=False, background = background_df,
            ylabel='information\n(bits)')
ax1.grid(False)
ax1.set_facecolor('white')
plt.setp(ax1.get_xticklabels(), visible=False)

anylogo.draw(ax2,effect_df=-pwm_df,logo_type='information',use_transparency=False,
            ylabel='information\n(bits)')
ax2.grid(False)
ax2.set_facecolor('white')
plt.setp(ax2.get_xticklabels(), visible=False)

plt.tight_layout()
fig6.savefig(output + 'figS4_fis_logo_compare.pdf')

correlation_fis = corr(pwm_df_fix,energy_df_fix)

#------------------------------------------------------------------------------#
# logos for -10 RNAP - compare to relB promoter RNAP binding site
#------------------------------------------------------------------------------#

# Get counts matrix and turn into pwm
counts_df = get_counts_seq('input_data/PromoterPredictionSigma70Set_SeqOnly_10.txt')

# Get number of CRP sites
num_sites = counts_df.iloc[0,:].sum()
print('num_sites = %d'%num_sites)

# Compute prob matrix
prob_df = (counts_df+0.1)/(num_sites+0.4)

# Compute pwm
background_df = prob_df.copy()
for i in range(len(prob_df)):
    background_df.loc[i,:] = background
pwm_df = np.log(prob_df/background_df)

pwm_df_fix = pwm_df.copy()

pwm_df_fix = pwm_df_fix[pwm_df_fix.index <= 5]
pwm_df_fix = -utils.fix_matrix_gauge(np.array(pwm_df_fix.copy()).T)


#-----

energy_df = pd.read_csv(datadir_rel + '20150513_relB_MG1655_M9glucose_na_mut1_4bins_RNAP_emat_mean.csv')
energy_df = energy_df[energy_df.position != energy_df.position.max()].reset_index()
# make background matrix for 6bp -10 and -35 regions
energy_df = energy_df[['A','C','G','T']]

background_df_emat = pwm_df.copy()
for i in range(30):
    background_df_emat.loc[i,:] = background

# for -10:
energy_df_fix = energy_df.copy()
energy_df_fix = energy_df_fix[energy_df_fix.index >= 24]
energy_df_fix = energy_df_fix[['A','C','G','T']]
energy_df_fix = utils.fix_matrix_gauge(np.array(energy_df_fix.copy()).T)


energy_df = utils.estimate_scalefactor(np.array(energy_df))*energy_df
energy_df.index.name = 'pos'


fig6, (ax1, ax2) = plt.subplots(2,1,figsize=[8,4])

anylogo.draw(ax1,effect_df=energy_df,logo_type='information',use_transparency=False, background = background_df_emat,
            ylabel='information\n(bits)')
ax1.grid(False)
ax1.set_facecolor('white')
plt.setp(ax1.get_xticklabels(), visible=False)

anylogo.draw(ax2,effect_df=-pwm_df,logo_type='information',use_transparency=False,
            ylabel='information\n(bits)')
ax2.grid(False)
ax2.set_facecolor('white')
plt.setp(ax2.get_xticklabels(), visible=False)

plt.tight_layout()
fig6.savefig(output + 'figS4_RNAP10_rel_logo_compare.pdf')

correlation_RNAP_10 = corr(pwm_df_fix,energy_df_fix)

#------------------------------------------------------------------------------#
# logos for RNAP - compare to relB promoter RNAP binding site
#------------------------------------------------------------------------------#

# Get counts matrix and turn into pwm
counts_df = get_counts_seq('input_data/PromoterPredictionSigma70Set_SeqOnly_35.txt')

# Get number of CRP sites
num_sites = counts_df.iloc[0,:].sum()
print('num_sites = %d'%num_sites)

# Compute prob matrix
prob_df = (counts_df+0.1)/(num_sites+0.4)

# Compute pwm
background_df = prob_df.copy()
for i in range(len(prob_df)):
    background_df.loc[i,:] = background
pwm_df = np.log(prob_df/background_df)

pwm_df_fix = pwm_df.copy()

pwm_df_fix = pwm_df_fix[pwm_df_fix.index <= 5]
pwm_df_fix = -utils.fix_matrix_gauge(np.array(pwm_df_fix.copy()).T)


#-----

energy_df = pd.read_csv(datadir_rel + '20150513_relB_MG1655_M9glucose_na_mut1_4bins_RNAP_emat_mean.csv')
energy_df = energy_df[energy_df.position != energy_df.position.max()].reset_index()
# make background matrix for 6bp -10 and -35 regions
energy_df = energy_df[['A','C','G','T']]

background_df_emat = pwm_df.copy()
for i in range(30):
    background_df_emat.loc[i,:] = background

# for -35:
energy_df_fix = energy_df.copy()
energy_df_fix = energy_df_fix[energy_df_fix.index <= 6]
energy_df_fix = energy_df_fix[energy_df_fix.index >= 1]
energy_df_fix = energy_df_fix[['A','C','G','T']]
energy_df_fix = utils.fix_matrix_gauge(np.array(energy_df_fix.copy()).T)



energy_df = utils.estimate_scalefactor(np.array(energy_df))*energy_df
energy_df.index.name = 'pos'


fig6, (ax1, ax2) = plt.subplots(2,1,figsize=[8,4])

anylogo.draw(ax1,effect_df=energy_df,logo_type='information',use_transparency=False, background = background_df_emat,
            ylabel='information\n(bits)')
ax1.grid(False)
ax1.set_facecolor('white')
plt.setp(ax1.get_xticklabels(), visible=False)

anylogo.draw(ax2,effect_df=-pwm_df,logo_type='information',use_transparency=False,
            ylabel='information\n(bits)')
ax2.grid(False)
ax2.set_facecolor('white')
plt.setp(ax2.get_xticklabels(), visible=False)

plt.tight_layout()
fig6.savefig(output + 'figS4_RNAP35_rel_logo_compare.pdf')

correlation_RNAP_35 = corr(pwm_df_fix,energy_df_fix)

#------------------------------------------------------------------------------#
# print correlation coefficients
#------------------------------------------------------------------------------#

print(correlation_LacI, 'lacI')
print(correlation_CRP, 'CRP')
print(correlation_MarA, 'MarA')
print(correlation_PurR, 'purR')
print(correlation_XylR_left, 'XylR_left')
print(correlation_XylR_right, 'XylR_right')
print(correlation_XylR_compare, 'XylR_left_vs_XylR_right')
print(correlation_LexA, 'LexA')
print(correlation_fis, 'Fis')
print(correlation_RNAP_10, 'RNAP -10')
print(correlation_RNAP_35, 'RNAP -35')
