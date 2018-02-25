import os
import glob
# Our numerical workhorses
import numpy as np
import pandas as pd
import scipy
# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm


from Bio import SeqIO
# Seaborn, useful for graphics
import seaborn as sns

#
# Import the project utils
import sys
sys.path.insert(0, '../')
sys.path.insert(0, '../../')
import NB_sortseq_utils as utils

import ConfigParser
config = ConfigParser.RawConfigParser()

#=============================================================================#
# Load in sequence file information.
# This extracts the information from a cfg file.
#=============================================================================#

##################
#config.read('cfg.py')
config.read(sys.argv[1])

date = config.get('Input','date')
promoter = config.get('Input','promoter')
mutregion = config.get('Input', 'mutregion')
strain = config.get('Input','strain')
media = config.get('Input','media')
condition = config.get('Input','condition')
bincount = config.getint('Input','bincount')
seq_pos = config.getint('Input','seq_pos')

Seq = config.get('Input','seq')
data_dir = config.get('Input','data_dir')
data_dir_out = config.get('Input','data_dir_out')
data_dir_out_summary = config.get('Input','data_dir_out_summary')
fseqname = config.get('Input','fseqname')

files_seq = glob.glob(data_dir + fseqname)
print(data_dir + fseqname)
len_seq = len(Seq)
print(files_seq)
#=============================================================================#
# Process sequence data

# naming convention:
# freq counts: data + strain + media + condition + (binning descriptor) '_freqCounts.csv'
# MI values: data + strain + media + condition + (binning descriptor) '_MI.csv'

#=============================================================================#
# Determine bin letter frequencies
#=============================================================================#

bin_freqs = np.zeros([len(files_seq), 4, len_seq])

if len(bin_freqs) <= 1:
    raise RuntimeError('Only one or no files found. There must be multiple binned\
    sequence files.')

for seqname in files_seq:
    binnumber = int(seqname[-7]) - 1
    seqtype = seqname[-5:]
    if not ('fasta'==seqtype or 'fastq'==seqtype):
        raise RuntimeError('Sequence is not of type .fastq or .fasta')
    original_reads = SeqIO.parse(seqname, seqtype)
    bin_freqs[binnumber,:,:] = utils.freqcount(original_reads)
    #check_reads = utils.freqcount(original_reads)

    # save frequency counts
    fname_out_template_MIfreqcounts = date + '_' + promoter + '_' + strain + '_' + media + \
                       '_' + condition +  '_' + mutregion + '_' + str(bincount)\
                        + 'bins' + '_bin' + str(binnumber) + '_MI_freqCounts.csv'
    np.savetxt(data_dir_out + fname_out_template_MIfreqcounts, bin_freqs[binnumber,:,:], delimiter=",")

#=============================================================================#
# Calculate information footprint MI in bits
#=============================================================================#
totreads, seqLength , MI_mu = utils.footprint_MI_calc(bin_freqs, finite_corr = True)

# save the total sequence count from each bin
fname_out_template_totreads = date + '_' + promoter + '_' + strain + '_' + \
                              media + '_' + condition +  '_' + mutregion +  \
                              '_' + str(bincount)  +  'bins' + '_bincounts.csv'

np.savetxt(data_dir_out + fname_out_template_totreads, totreads, delimiter=",")

# save the MI values for each position
fname_out_template_MI = date + '_' + promoter + '_' + strain + '_' + media + \
                       '_' + condition +  '_' + mutregion + '_' + str(bincount)\
                        + 'bins' + '_MI.csv'

np.savetxt(data_dir_out + fname_out_template_MI, MI_mu, delimiter=",")

print('Calculated information footprint.')

#=============================================================================#
# Calculate mutation rate of library
#=============================================================================#
mut_rate = utils.calc_mutrate(bin_freqs)

fname_out_template_mutrate = date + '_' + promoter + '_' + strain + '_' + \
                            media + '_' + condition +  '_' +  mutregion + '_' \
                            + str(bincount) + 'bins' + '_mutrate.csv'
np.savetxt(data_dir_out + fname_out_template_mutrate, mut_rate, delimiter=",")
print('Calculated mutation rate.')

#=============================================================================#
# Calculate expression shifts
#=============================================================================#

avgBin_counts, avgBin, avgbin_WT = utils.calc_deltabin(Seq, files_seq, bin_freqs,seqtype)

fname_out_deltabin_counts = date + '_' + promoter + '_' + strain + '_' + media \
                            + '_' + condition +  '_' + mutregion + '_' + \
                            str(bincount) + 'bins' + '_expshiftcounts.csv'
fname_out_deltabin = date + '_' + promoter + '_' + strain + '_' + media + '_' \
                    + condition +  '_' + mutregion + '_' + str(bincount) + \
                    'bins' + '_expshift.csv'
fname_out_deltabin_WTavg = date + '_' + promoter + '_' + strain + '_' + media + '_' \
                    + condition +  '_' + mutregion + '_' + str(bincount) + \
                    'bins' + '_expshiftWTavg.csv'

np.savetxt(data_dir_out + fname_out_deltabin_counts, avgBin_counts, \
            delimiter=",", fmt="%s")
np.savetxt(data_dir_out + fname_out_deltabin, avgBin, delimiter=",", \
            fmt="%s")
np.savetxt(data_dir_out + fname_out_deltabin_WTavg, np.array(avgbin_WT).reshape(1,))

print('Calcualted expression shifts')

#=============================================================================#
# Calculate expression shifts, 3 bp running avg
#=============================================================================#
avgBin_counts_3bpavg, avgBin_3bpavg = utils.calc_deltabin_3bpavg(Seq, \
                                            files_seq, bin_freqs, seqtype)
fname_out_deltabin_counts_3bpavg = date + '_' + promoter + '_' + strain + '_' \
                                   + media + '_' + condition +  '_' + mutregion \
                                   + '_' + str(bincount) + 'bins' + \
                                   '_expshift_3bpavg.csv'
fname_out_deltabin_3bpavg = date + '_' + promoter + '_' + strain + '_' + media \
                            + '_' + condition +  '_' + mutregion + '_' + \
                            str(bincount) + 'bins' + '_expshift_3bpavg.csv'
np.savetxt(data_dir_out + fname_out_deltabin_counts_3bpavg, avgBin_counts_3bpavg, \
            delimiter=",", fmt="%s")
np.savetxt(data_dir_out + fname_out_deltabin_3bpavg, avgBin_3bpavg, delimiter=",", \
            fmt="%s")
print('Calcualted expression shifts, 3bpavg.')


#=============================================================================#
# Combine all data in a tidy pandas and save csv
# generate pandas dataframe
# columns: date, strain, media, condition, bincount, position, MI, MI_error,
# expshift, avg_fluor_bin1,..., WT_bp, mutation_rate, number_reads,
#=============================================================================#

# initialize the DataFrame to save the mean expression levels
df = pd.DataFrame()

# read the files and compute the mean YFP value
for i, let in enumerate(Seq):
    try:
        df = df.append([[date, promoter, strain, media, condition, mutregion, \
                bincount, (seq_pos+ i), Seq[i], MI_mu[i], mut_rate[i], avgBin[i], \
                avgBin_3bpavg[i], files_seq ]], ignore_index=True)
    except:
        pass

# rename the columns of the data_frame
df.columns = ['date', 'promoter', 'strain', 'media', 'condition', 'mutregion', \
              'bincount', 'position', 'WT_bp', 'MI', 'mutation_rate', \
              'expshift', 'expshift_3bpavg', 'seq_files']

df.to_csv(data_dir_out_summary + date + '_' + promoter + '_' + strain + '_' + media + '_' + \
                             condition +  '_' + mutregion + '_' + str(bincount)\
                              + 'bins' + '_summary.csv', index=False)
