import os
import glob
# Our numerical workhorses
import numpy as np
import pandas as pd
import scipy
# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pymc

# Seaborn, useful for graphics
#import seaborn as sns

#
# Import the project utils
import sys
sys.path.insert(0, '../../')
sys.path.insert(0, '../../sortseq_pymc/')

import NB_sortseq_utils as utils
utils.set_plotting_style1()

import ConfigParser
config = ConfigParser.RawConfigParser()

#=============================================================================#
# Load in file information. This will find all .sql files.
# Since each MCMC was run many times (i.e. 20 times),
# this will use the file names to identify each set of runs to process
# each set together.
#=============================================================================#

############## everything that needs to be changed goes here ##################

config.read(sys.argv[1])
emat_choice = sys.argv[2]

date = config.get('Input','date')
promoter = config.get('Input','promoter')
mutregion = config.get('Input', 'mutregion')
strain = config.get('Input','strain')
media = config.get('Input','media')
condition = config.get('Input','condition')
bincount = config.getint('Input','bincount')
Seq = config.get('Input','seq')

emat_dir_seq = config.get('Input','emat_dir_seq')
emat_dir_sql = config.get('Input','emat_dir_sql')
emat_dir_out = config.get('Input','emat_dir_out')
emat_fseqname = config.get('Input','emat_fseqname')

mut_region_start = config.getint(emat_choice,'mut_region_start')
mut_region_length = config.getint(emat_choice,'mut_region_length')
TF = config.get(emat_choice,'TF')
TF_type = config.getint(emat_choice,'TF_type')

# target file
fname_sql = date + '_' + promoter + '_' + strain + '_' + media + \
                   '_' + condition +  '_' + mutregion + '_' + str(bincount)\
                    + 'bins' + '_pymc_' + TF +'_MCMC'
fname_seq = date + '_' + promoter + '_' + strain + '_' + media + \
                   '_' + condition +  '_' + mutregion + '_' + str(bincount)\
                    + 'bins' + '_pymc'

burn_in = 1000

#------------------------------------------------------------------------------#
#load in the sequence data
#------------------------------------------------------------------------------#
data_fn = os.path.join(emat_dir_seq, emat_fseqname)
seq_mat, batch_vec = utils.load_unique_seqs_batches(data_fn, mut_region_start, mut_region_length)
energies = np.zeros(len(batch_vec))

# make dir for storing MCMC run data_fn
if not os.path.exists('emat_processing/'):
    os.makedirs('emat_processing')
#------------------------------------------------------------------------------#
# find sql MCMC trace files
#------------------------------------------------------------------------------#
listing_emat = glob.glob(os.path.join(emat_dir_sql, (fname_sql + '*.sql')))
print(listing_emat)
#------------------------------------------------------------------------------#
# Start analysis of sql MCMC trace files
#------------------------------------------------------------------------------#
# if these files have been processed, I don't want it to append the MI.txt
# this will delete the entries in the old file
MI_f = open(os.path.join(emat_dir_out,fname_sql + '_MI.txt'),'w')
MI_f.write('')
MI_f.close()


# initialize array to hold trace values and counter for finding average of
# energy matrices across MCMC traces
trace_emats = []
count = 0
for fn in listing_emat:
    run_number = fn[-8:-6]
    db = pymc.database.sqlite.load(fn)
    trace_temp = db.trace('emat')[burn_in:]

    # Step 1: determine matrix entried, fix gauge, and check sign of energy values
    emat_mean_temp = np.mean(trace_temp,axis=0)
    emat_mean_temp = utils.fix_matrix_gauge(emat_mean_temp)
    for i in range(len(batch_vec)):
        energies[i] = np.sum(seq_mat[:,:,i]*emat_mean_temp)
    r = scipy.stats.pearsonr(energies,batch_vec)[0]
    if r>0:
        trace_temp = -trace_temp * TF_type
        emat_mean_temp = -emat_mean_temp * TF_type
    else:
        trace_temp = trace_temp * TF_type
        emat_mean_temp = emat_mean_temp * TF_type

    # # Step 2: fix gauge (i.e. energy scales) such that entried have mean of
    # # zero and standard deviation equal to one.
    # emat_mean_temp = utils.fix_matrix_gauge(emat_mean_temp)

    # Step 3: Save matrix data from this run (emat, joint dist, and PDF images)
    np.savetxt(os.path.join(emat_dir_out,fname_sql + '_' + run_number \
                +'_emat_mean.txt'),emat_mean_temp)


    # compute mutual information and joint pdf for this MCMC
    MI,f_reg = utils.compute_MI(seq_mat,batch_vec,emat_mean_temp.copy())
    MI_f = open(os.path.join(emat_dir_out,fname_sql + '_MI.txt'),'a')
    MI_f.write(str(MI))
    MI_f.write('\n')

    # make heat map of mean matrix for this MCMC
    plt.clf()
    plt.imshow(utils.zero_matrix_WT(emat_mean_temp.copy(),Seq),
               interpolation='nearest', cmap = 'RdBu_r')
    plt.title('Energy Matrix, ' + run_number + ', MI: %.5f' % MI)
    plt.grid(b=False)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(emat_dir_out,fname_sql + '_' + run_number+'.pdf'))

    # make joint distribution plot for this MCMC
    plt.clf()
    plt.imshow(f_reg,interpolation='nearest',aspect='auto', cmap = 'Blues')
    plt.title('Joint regularized pdf, ' + run_number + ', MI: %.5f' % MI)
    plt.xlabel('Rank order')
    plt.ylabel('Batch number')
    plt.tight_layout()
    plt.savefig(os.path.join(emat_dir_out,fname_sql + '_' + run_number+'_regpdf.pdf'))

    # If this is first loop, replace trace = [] with
    # trace_temp AND don't append energy matrix if the correlation r is poor.
    # I'm going to re-calculate the energy matrix based on sign corrected trace
    # and re-fix gauge.
    emat_mean_temp = np.mean(trace_temp,axis=0)
    emat_mean_temp = utils.fix_matrix_gauge(emat_mean_temp)

    # Step 4: Append energy matrix; here I'm just adding and will divide by
    # 'count' at the end to determine the average energy matrix values across
    # MCMC runs.
    if trace_emats == []:
        trace_emats = emat_mean_temp
        count += 1
    else:
        # #check that the matrix entries have the same sign
        r = scipy.stats.pearsonr(trace_emats.flatten(),emat_mean_temp.flatten())[0]
        print(r)
        if np.absolute(r) <= 0.85:
            print('MCMC trace not used')
            continue
        trace_emats += emat_mean_temp
        count += 1

#------------------------------------------------------------------------------#
# Step 5:  Calculate average energy matrix from all MCMC traces analyzed
# Combine all data in a tidy pandas and save csv
# Also calculate MI value and joint probability using everage energy matrix
#------------------------------------------------------------------------------#
emat_mean = trace_emats/count
emat_mean = utils.fix_matrix_gauge(emat_mean)

#=============================================================================#

# # initialize a DataFrame to save the experimental details
df = pd.DataFrame()
pos = np.arange(mut_region_start,mut_region_start + mut_region_length)
# read the files and compute the mean YFP value

for i in enumerate(pos):
    try:
        df = df.append([[date, promoter, strain, media, condition, mutregion, \
                bincount, emat_choice, i[1], Seq[i[1]]]], ignore_index=True)
    except:
        pass
#
# rename the columns of the data_frame
df.columns = ['date', 'promoter', 'strain', 'media', 'condition', 'mutregion', \
              'bincount', 'TF', 'position', 'WT_sequence']

emat_df = pd.DataFrame(emat_mean.T,columns=list('ACGT'))

df_summary = pd.concat([emat_df,df], axis=1)

fname_out = date + '_' + promoter + '_' + strain + '_' + media + '_' + \
                             condition +  '_' + mutregion + '_' + str(bincount)\
                              + 'bins_' + emat_choice
df_summary.to_csv(fname_out + '_emat_mean.csv',
                              index=False)


#--------------------------------------#
#
# np.savetxt(fname_sql + '_emat_mean.txt',emat_mean)

# compute mutual information and joint pdf
MI,f_reg = utils.compute_MI(seq_mat,batch_vec,emat_mean)

# make heat map of mean matrix
fig1 = plt.figure()
ax1 = plt.subplot(111)
ax1.imshow(utils.zero_matrix_WT(emat_mean,Seq),interpolation='nearest',cmap = 'RdBu_r')
ax1.set_title('Energy Matrix, ' + ', MI: %.5f' % MI)
ax1.axis('off')
ax1.grid(b=False)
fig1.savefig(os.path.join(fname_out + '_emat_mean.pdf'))

fig2 = plt.figure()
ax2 = plt.subplot(111)
ax2.imshow(f_reg,interpolation='nearest',aspect='auto', cmap = 'Blues')
ax2.set_title('Joint regularized pdf, ' + ', MI: %.5f' % MI)
ax2.set_xlabel('Rank order')
ax2.set_ylabel('Batch number')
plt.tight_layout()
fig2.savefig(os.path.join(fname_out + '_regpdf.pdf'))
