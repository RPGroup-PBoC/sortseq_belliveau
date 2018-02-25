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
mutregion = 'mut1'
strain = 'MG1655'
media = 'M9glucose'
condition = 'na'
bincount = 4
TSS = 140

data_dir = '../../data/sortseq_raw/20160707_NB_003_pD_X_P/seq/'

# entire wild type sequence in file
Seq = "TTTTATGCATTGTTCTTTTTGTGATCTAAATTGTAGTACAACAATATAAGTTTGTACTACATTACACGCACGGCAAACGCGAACGTCATCACGCTGGTACTACAAAGTTGCCGCGTTATGCATCGATCGGGGTAAAGTAGAGAAGAACATACAGAG"

seqLength = len(Seq)

colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]
ind = np.arange(seqLength)              # the x locations for the groups

fseqname = '20160710_MG1655_pdgoR_glucose_merged.mut1.bin*.filt.fastq'

files_seq = [['20160710_MG1655_pdgoR_glucose_merged.mut1.bin1.filt.fastq'],\
            ['20160710_MG1655_pdgoR_glucose_merged.mut1.bin2.filt.fastq'],\
            ['20160710_MG1655_pdgoR_glucose_merged.mut1.bin3.filt.fastq'],\
            ['20160710_MG1655_pdgoR_glucose_merged.mut1.bin4.filt.fastq']]#glob.glob(data_dir + fseqname)

############
# generate pandas dataframe
# columns: date, strain, media, condition, bincount, position, MI, MI_error, deltabin, avg_fluor_bin1,..., WT_bp, mutation_rate, number_reads,
# initialize the DataFrame to save the mean expression levels
MI_mu =         np.loadtxt('20160710_dgoR_MG1655_M9glucose_na_mut1_4bins_MI.csv')
mut_rate =      np.loadtxt('20160710_dgoR_MG1655_M9glucose_na_mut1_4bins_mutrate.csv')
bincounts =     np.loadtxt('20160710_dgoR_MG1655_M9glucose_na_mut1_4bins_bincounts.csv')
avgBin =        np.loadtxt('20160710_dgoR_MG1655_M9glucose_na_mut1_4bins_deltabin.csv')
avgBin_3bpavg = np.loadtxt('20160710_dgoR_MG1655_M9glucose_na_mut1_4bins_deltabin_3bpavg.csv')

df = pd.DataFrame()

# read the files and compute the mean YFP value
for i, let in enumerate(Seq):
    try:
        df = df.append([[date, promoter, strain, media, condition, mutregion, \
                bincount, i, Seq[i], MI_mu[i], mut_rate[i], avgBin[i], \
                avgBin_3bpavg[i], files_seq ]], ignore_index=True)
    except:
        pass

# rename the columns of the data_frame
df.columns = ['date', 'promoter', 'strain', 'media', 'condition', 'mutregion', \
              'bincount', 'position', 'WT_bp', 'MI', 'mutation_rate', \
              'delta_bin', 'delta_bin_3bpavg', 'seq_files']

df.to_csv(date + '_' + promoter + '_' + strain + '_' + media + '_' + \
                             condition +  '_' + mutregion + '_' + str(bincount)\
                              + 'bins' + '_summary.csv', index=False)
