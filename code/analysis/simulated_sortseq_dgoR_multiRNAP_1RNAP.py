# Here we will similate a Sort-Seq experiment for pdgoR with the hypothesis that
# there three overlapping RNAP binding sites

#load tools
import scipy as sp
import numpy as np
import pandas as pd
from Bio import SeqIO
import math
import matplotlib.pyplot as plt
from IPython.core.pylabtools import figsize
import sys
sys.path.insert(0, '../')
import NB_sortseq_utils as utils


datadir = 'input_data/'


def seq_mut(seq):
    """
    Input: sequence to be mutated
    Output: mutated sequence where mutation rate is defined by
    p, where the probability of each transition is 25% (with the
    transition to itself not allowed).
    """
    seq_dict = {0:'A',1:'C',2:'G',3:'T'}
    seq_mut = ''
    for i,bp in enumerate(seq):
        p = np.random.uniform(low=0.0, high=1.0)
        if p >= 0.9:
            seq_temp = seq[i]
            while seq[i] == seq_temp:
                t = np.random.randint(low =0,high=4)
                seq_temp = seq_dict[t]
            seq_mut += seq_temp
        else:
            seq_mut += seq[i]
    return seq_mut


def freqcount(records):
    """counts the frequency of each nucleotide at each position along the sequences considered
    """
    count = 0
    len_seq = 64#len(records[0].seq)
    global countletter
    countletter = np.zeros((4, len_seq))

    for record in records:
        count += 1
        for x in range(0, len_seq):
            if record.seq[x] == "A":
               countletter[0,x] += 1
            elif record.seq[x] == "C":
               countletter[1,x] += 1
            elif record.seq[x] == "G":
               countletter[2,x] += 1
            elif record.seq[x] == "T":
               countletter[3,x] += 1
            else:
                print("error: unexpected letter")
#------------------------------------------------------------------------------#
# make mutated library
#------------------------------------------------------------------------------#

RNAP_full = 'TTTTTATGCATTGTTCTTTTTGTGATCTAAATTGTAGTACAACAATATAAGTTTGTACTACATTA'

RNAP_WT = 'TTTTTATGCA'
RNAP_mut_region = 'TTGTTCTTTTTGTGATCTAAATTGTAGTACAACAATATAAGTTTGTACTACATTA'

wt_pdgor_mut_array = np.loadtxt(datadir +
       '20170814_dgoR_sequences_10percentMut.txt' ,dtype='str')


#------------------------------------------------------------------------------#
# calculate energies based on RNAP matrix
# #------------------------------------------------------------------------------#
ematdir = 'input_data/'
emat = '20150513_relB_MG1655_M9glucose_na_mut1_4bins_RNAP_emat_mean_plusspacer1bp.csv'
emat_promoters = 'relB'

# Load in energy matrix
emat_RNAP = pd.read_csv(ematdir + emat)

emat_RNAP = scalefactor*utils.fix_matrix_gauge(np.array(emat_RNAP[['A','C','G','T']].T)[:,1:31])
# we're going to ignore the bp on each end (padding)
len_emat = 30

# collect 5*10**6 per bin
N_seq = 5*10**6

# compute energies for each operator sequence
#seqs1 = sp.zeros((4,30,N_seq),dtype=int)
# seqs2 = sp.zeros((4,30,N_seq),dtype=int)
seqs3 = sp.zeros((4,30,N_seq),dtype=int)

for j in range(0,N_seq):
    #seqs1[:,:,j] = utils.seq2mat(wt_pdgor_mut_array[j][3:33])
    # seqs2[:,:,j] = utils.seq2mat(wt_pdgor_mut_array[j][13:43])
    seqs3[:,:,j] = utils.seq2mat(wt_pdgor_mut_array[j][34:-1])


# dot1 = np.array(emat_RNAP)[:,:,sp.newaxis]*seqs1
# energies_RNAP1 = dot1.sum(0).sum(0)

# dot2 = np.array(emat_RNAP)[:,:,sp.newaxis]*seqs2
# energies_RNAP2 = dot2.sum(0).sum(0)

dot3 = np.array(emat_RNAP)[:,:,sp.newaxis]*seqs3
energies_RNAP3 = dot3.sum(0).sum(0)


# calculate expression (i.e. pbound \propto P_1/N exp(energies_1) +
# P_2/N exp(energies_2) + ... )

#Lets place these arrays in a Pandas dataframe for ease of handling
df =  pd.DataFrame(wt_pdgor_mut_array)
df.columns = ['RNAP full Sequence']
# df['RNAP site 1 energy prediction'] = pd.DataFrame(energies_RNAP1)
# df['RNAP site 2 energy prediction'] = pd.DataFrame(energies_RNAP2)
df['RNAP site 3 energy prediction'] = pd.DataFrame(energies_RNAP3)


#lets add some noise to make it more realizstic
P = 3000 # approx number for RNAP
df['P'] = np.random.normal(P, P*0.25,len(df['RNAP full Sequence']))

# calculate Pbound
# p1_state = np.exp(-df['RNAP site 1 energy prediction'])
# p2_state = np.exp(-df['RNAP site 2 energy prediction'])
p3_state = np.exp(-df['RNAP site 3 energy prediction'])


N_ns = 4.6E6
df['tau'] = (df['P']/N_ns) * (p3_state)/ (1 +  (df['P']/N_ns) *( p3_state ))

#------------------------------------------------------------------------------#
# sort cells into four bins of 15% windows
#------------------------------------------------------------------------------#

# identify the log(expression) percentile values for each of our 'FACS' sorts.
tau = 'tau'
p1 = np.percentile(df[tau], 15)
p2 = np.percentile(df[tau], 25)
p3 = np.percentile(df[tau], 40)
p4 = np.percentile(df[tau], 60)
p5 = np.percentile(df[tau], 75)
p6 = np.percentile(df[tau], 85)

df_bin1 = df[df[tau] <= p1]
df_bin2 = df[(df[tau] >= p2) & (df[tau] <= p3)]
df_bin3 = df[(df[tau] >= p4) & (df[tau] <= p5)]
df_bin4 = df[df[tau] >= p6]


#now we need to save our binned data as fasta files for later analysis
#I might be making more work for myself by converting to _scalefactor3.fasta

fname = '20170822_pdgor_1RNAP_' + emat_promoters[i] + 'bin'
#since the library sequences were randomly generated, I'm just going to
#take the first 200,000 sequences in eachbin.
count = 0
file = open(datadir + 'sortedseq/' + fname + '1.fasta', "w")
while count <= 2*10**5:
    file.write(">")
    file.write(str(count))
    file.write("\n")
    file.write(df_bin1['RNAP full Sequence'].values[count][10:])
    file.write("\n")
    count += 1
file.close()

count = 0
file = open(datadir + 'sortedseq/' + fname + '2.fasta', "w")
while count <= 2*10**5:
    file.write(">")
    file.write(str(count))
    file.write("\n")
    file.write(df_bin2['RNAP full Sequence'].values[count][10:])
    file.write("\n")
    count += 1
file.close()

count = 0
file = open(datadir + 'sortedseq/' + fname + '3.fasta', "w")
while count <= 2*10**5:
    file.write(">")
    file.write(str(count))
    file.write("\n")
    file.write(df_bin3['RNAP full Sequence'].values[count][10:])
    file.write("\n")
    count += 1
file.close()

count = 0
file = open(datadir + 'sortedseq/' + fname + '4.fasta', "w")
while count <= 2*10**5:
    file.write(">")
    file.write(str(count))
    file.write("\n")
    file.write(df_bin4['RNAP full Sequence'].values[count][10:])
    file.write("\n")
    count += 1
file.close()
