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

#------------------------------------------------------------------------------#
# make mutated library
#------------------------------------------------------------------------------#
# region that contains upstream RNAP site(s) on dogR promoter
RNAP_full = 'TTTTTATGCATTGTTCTTTTTGTGATCTAAATTGTAGTACAACAATATAAGTTTGTACTACATTA'

# the first 10bp were not mutated in experiments
RNAP_WT = 'TTTTTATGCA'
RNAP_mut_region = 'TTGTTCTTTTTGTGATCTAAATTGTAGTACAACAATATAAGTTTGTACTACATTA'

#Lets generate our library of mutated sequences ; of size N_seq
N_seq = 5*10**6
wt_pdgor_mut_array = ["" for x in range(N_seq)]
for i,bp in enumerate(wt_pdgor_mut_array):
    wt_pdgor_mut_array[i] = RNAP_WT + seq_mut(RNAP_mut_region)

#Save to txt so that we don't have to rerun this
np.savetxt(datadir + '20170814_dgoR_sequences_10percentMut.txt',
           wt_pdgor_mut_array , fmt="%s")
