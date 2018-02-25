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

# Seaborn, useful for graphics
import seaborn as sns

sns.set_palette("deep", color_codes=True)
utils.set_plotting_style4()

#===============================================================================
# Set output directory
#===============================================================================
output = 'output_figs/'

#------------------------------------------------------------------------------#
# Load in the summary emat csv file
#------------------------------------------------------------------------------#
# 20170818 xylE MCMC fit using MPAthic
datadir = '../sortseq/20160707_purT_xylE_dgoR/'

emat_xylE = pd.read_csv(datadir + '20160710_xylE_MG1655_M9xylose_na_mut2_4bins_XylR_emat_mean.csv')


emat_xylE = emat_xylE[['A','C','G','T']]
emat_xylE = np.array(emat_xylE.T)#[:,5:]
#emat_xylE = utils.fix_matrix_gauge(emat_xylE)
len_emat = len(emat_xylE[0,:])
#------------------------------------------------------------------------------#
# make shared axis plots of information footprint and delta bin shift
#------------------------------------------------------------------------------#

# xylFG promoter
seq = 'ATTGAACTCCATAATCAGGTAATGCCGCGGGTGATGGATGATGTCGTAATATTGGGCACTCCCTTTCAGTTGCTCAATTATGTTATTTCACACTGCTATTGAGATAATTCACAAGTGTGCGCTCGCTCGCAAAATAAAATGGAATGATGAAACTGGGTAATTCCTCGAAGAGAAAAATGCAATAAGTACAATTGCGCAACAAAAGTAAGATCTCGGTCATAAATCAAGAAATAAACCAAAAATCGTAATCGAAAGATAAAAATCTGTAATTGTTTTCCCCTGTTTAGTTGCTAAAAATTGGTTACGTTTATCGCGGTGATTGTTACTTATTAAAACTGTCCTCTAACTACAGAAGGCCCTACACC'
seq_ind = np.arange(len(seq)-len_emat) - len(seq)
energies_scan_xylE = np.zeros(len(seq_ind))
fig, (ax1, ax2, ax3) = plt.subplots(3,figsize=(2,1.5), sharey=True, sharex=True)
#figsize=utils.cm2inch(4,5)

for i in range(0,len(seq)-len_emat):
    # mat = utils.genomescan(seq[i:(i+len_emat)],emat_xylE)
    mat = utils.seq2mat(seq[i:(i+len_emat)])
    energies_scan_xylE[i] = np.sum(np.multiply(mat, emat_xylE))
    if energies_scan_xylE[i] < -2.3:
        print(seq[i:(i+len_emat)])
        print(energies_scan_xylE[i])
ax1.plot(seq_ind,energies_scan_xylE)
# ax1.set_xlabel('position from start codon (bp)', fontsize = 5)
ax1.set_ylabel('predicted\nenergy (a.u.)', fontsize = 5)
ax1.tick_params(labelsize=5)
ax1.yaxis.grid(False)
ax1.set_xlim(seq_ind[0],0)

ax1.set_yticks([-3.0,-2.0, -1.0, 0.0, 1.0])
# ax1.set_xticks(np.arange(min(energies_scan_xylE), max(energies_scan_xylE)+1, 0.1))
# ax1.tick_params(fontsize=5)

# xylAB promoter (reverse complement of xylFG)
seq = 'GGTGTAGGGCCTTCTGTAGTTAGAGGACAGTTTTAATAAGTAACAATCACCGCGATAAACGTAACCAATTTTTAGCAACTAAACAGGGGAAAACAATTACAGATTTTTATCTTTCGATTACGATTTTTGGTTTATTTCTTGATTTATGACCGAGATCTTACTTTTGTTGCGCAATTGTACTTATTGCATTTTTCTCTTCGAGGAATTACCCAGTTTCATCATTCCATTTTATTTTGCGAGCGAGCGCACACTTGTGAATTATCTCAATAGCAGTGTGAAATAACATAATTGAGCAACTGAAAGGGAGTGCCCAATATTACGACATCATCCATCACCCGCGGCATTACCTGATTATGGAGTTCAAT'
seq_ind = np.arange(len(seq)-len_emat) - len(seq)
energies_scan_xylE = np.zeros(len(seq_ind))
for i in range(0,len(seq)-len_emat):
    # mat = utils.genomescan(seq[i:(i+len_emat)],emat_xylE)
    mat = utils.seq2mat(seq[i:(i+len_emat)])
    energies_scan_xylE[i] = np.sum(np.multiply(mat, emat_xylE))
    if energies_scan_xylE[i] < -2.3:
        print(seq[i:(i+len_emat)])
        print(energies_scan_xylE[i])
ax2.plot(seq_ind,energies_scan_xylE)
# ax2.set_xlabel('position from start codon (bp)', fontsize = 5)
ax2.set_ylabel('predicted\nenergy (a.u.)', fontsize = 5)
ax2.tick_params(labelsize=5)
ax2.yaxis.grid(False)
ax2.set_xlim(seq_ind[0],0)


# xylE promoter
seq = 'TCAACGGCCTGACGGCAGGTGGTGTGAAAGGTTAAAGATGTTGTTCTGCCAATGTTATGCCGCTGCACCCTCAACTTACGTTATCCCAACTTGTGACTGTTATTCGGCGCTCCACGGAGCGCCTTTTTTTCTTTCGTCTGCAATCTGAATCGTTCGCCGGTTAATATTTCCATCATAGAGCTTATTATTTTTACGTTATTTGTTTTCCCACTTACGATAATTCTCTTTCGTGCTCTGAGTCACGGCAATAGTATTGTTTTTATCAATTTTGGATAATTATCACAATTAAGATCACAGAAAAGACATTACGTAAACGCATTGTAAAAAATGATAATTGCCTTAACTGCCTGACAATTCCAACATCAATGCACTGATAAAAGATCAGAATGGTCTAAGGCAGGTCTGA'
seq_ind = np.arange(len(seq)-len_emat) - len(seq)
energies_scan_xylE = np.zeros(len(seq_ind))
for i in range(0,len(seq)-len_emat):
    # mat = utils.genomescan(seq[i:(i+len_emat)],emat_xylE)
    mat = utils.seq2mat(seq[i:(i+len_emat)])
    energies_scan_xylE[i] = np.sum(np.multiply(mat, emat_xylE))
    if energies_scan_xylE[i] < -2.3:
        print(seq[i:(i+len_emat)])
        print(energies_scan_xylE[i])
ax3.plot(seq_ind,energies_scan_xylE)
ax3.set_xlabel('position from start codon (bp)', fontsize = 5)
ax3.set_ylabel('predicted\nenergy (a.u.)', fontsize = 5)
ax3.tick_params(labelsize=5)
ax3.yaxis.grid(False)
ax3.set_xlim(seq_ind[0],0)
# ax3.text(seq_ind[10],0.55, r'$xylE$ promoter',
#         horizontalalignment='center',
#         verticalalignment='center', alpha = 0.8, fontsize='8')
#
# # dgoR promoter
# seq = 'GCAGCCTTAAAAAATCCCTTAAAATTCAACTTCATGATTAATCCTTTAAACGAGTAGCTAAGCCGTTGAAAGCATACTCTTCTTTCCCTCCTGTTTAGCGCCTGATATCCCTTTTCAGATTTCTGCCCGACGCATGTCATTTTTTTATGCATTGTTCTTTTTGTGATCTAAATTGTAGTACAACAATATAAGTTTGTACTACATTACACGCACGGCAAACGCGAACGTCATCACGCTGGTACTACAAAGTTGCCGCGTTATGCATCGATCGGGGTAAAGTAGAGAAGAACATACAGAGCACAAGGACTCTCC'
# seq_ind = np.arange(len(seq)) - len(seq)
# energies_scan_xylE = np.zeros(len(seq))
# for i in range(0,len(seq)):
#     mat = utils.genomescan(seq[i:(i+len_emat)],emat_xylE)
#     energies_scan_xylE[i] = np.sum(np.multiply(mat, emat_xylE))
#     if energies_scan_xylE[i] < -2.0:
#         print(seq[i:(i+len_emat)])
# ax4.plot(seq_ind,energies_scan_xylE)
# # purT promoter
# seq = 'GAATGACTACGTATTTAACTTCAACCGCCATTTGCAGCCTCTCATAATAACTGTGATTTTATACAGTATATTTCTTTTCGGTTGAGAAATCAACATCAGCAATAAAGACACACGCAAACGTTTTCGTTTATACTGCGCGCGGAATTAATCAGGGGATATTCGTTATGACGTTATTAGGCACTGCGCTGCGTCCGGCAGC'
# seq_ind = np.arange(len(seq)) - len(seq)
# energies_scan_xylE = np.zeros(len(seq))
# for i in range(0,len(seq)):
#     mat = utils.genomescan(seq[i:(i+len_emat)],emat_xylE)
#     energies_scan_xylE[i] = np.sum(np.multiply(mat, emat_xylE))
#     if energies_scan_xylE[i] < -2.0:
#         print(seq[i:(i+len_emat)])
# ax4.plot(seq_ind,energies_scan_xylE)
# #lacZ promoter
# seq = 'AGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGAAACAGCT'
# seq_ind = np.arange(len(seq)) - len(seq)
# energies_scan_xylE = np.zeros(len(seq))
# for i in range(0,len(seq)):
#     mat = utils.genomescan(seq[i:(i+len_emat)],emat_xylE)
#     energies_scan_xylE[i] = np.sum(np.multiply(mat, emat_xylE))
#     if energies_scan_xylE[i] < -2.0:
#         print(seq[i:(i+len_emat)])
# ax4.plot(seq_ind,energies_scan_xylE)
#
# ax4.set_xlabel('position from start codon (bp)', fontsize = 5)
# ax4.set_ylabel('predicted binding\nenergy (a.u.)', fontsize = 5)
# ax4.tick_params(labelsize=5)
# ax4.yaxis.grid(False)

fig.subplots_adjust(hspace=0.1)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

# plt.tight_layout()
figname_out = 'figS8_xylR_matrix_scans.pdf'
fig.savefig(output + figname_out, format='pdf')
