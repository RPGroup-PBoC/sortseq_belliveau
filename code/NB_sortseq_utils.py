"""
Title:
    NB_sortseq_utils.py
Creation Date:
    2017-07-01
Author(s):
    Nathan Belliveau
Purpose:
    This file contains a variety of functions used for data processing, analysis,
    and inferences. The functions are split into different categories and we
    list them below.

        Sort-Seq processing (of .fasta or .fastq files)
        -------------------------------
        freqcount:
            Returns the number of A, C, G, Ts at each positions in the provided
            sequence file.

        freqcount_cond:
            Returns the number of A, C, G, Ts at each positions in the provided
            sequence file - for a subset of sequences that have letter 'Let'
            at the specified position.

        footprint_MI_calc:
            Returns the information footprint mutual information values based
            on the frequencies of each A, C, G, T along the mutated region
            considered.

        footprint_KLD_calc:
            Returns the Kullback-Leibler divergence between the nucleotide
            frequencies (A, C, G, T ) of two growth conditions.

        calc_mutrate:
            Returns the mutation rate at each position along the mutated region.

        calc_deltabin:
            Returns the average expression shift at each position (i.e.
            the average bin of mutated sequences relative to the average of the
            wild-type sequence)

        calc_deltabin_3bpavg:
            Returns the average expression shift as above, but by considering
            three positions at a time.

        Energy matrix and sequence logo tools
        -------------------------------

        seq2mat:
            Returns a 4xL array representation of a sequences (1=WT, 0 otherwise)

        zero_matrix:
            Returns a 4xL matrix model where the energies in each column have
            been shifted such that smallest element in each column has zero
            energy.

        zero_matrix_WT:
            Returns a 4xL matrix model where the energies in each column have
            been shifted such that the wild-type nucleotide has zero energy.

        compute_energies:
            Multiplies energy matrix with sequence matrix (array of 1 or more)
            to calculate the energies.

        fix_matrix_gauge:
            Returns a 4xL matrix model into the canonical gauge.

        estimate_scalefactor:
            Returns scaled matrix model in which each column's minimmum entry
            shifted to zero minimum energy equal to zero and then scaled
            such that the average energetic cost of a mutation is 2.5 k_B T.

        emattoIUPAC:
            Returns IUPAC sequence as determined from a matrix model

        hpd:
            Returns highest probability density region given by a set of
            samples.

        MI tools for MCMC matrix processing
        -------------------------------

        load_seqs_batches:
            Loads in sequences in data, to be used to check correlation
            between energy predictions and expected sign (e.g. negative repressor
            binding energy should correspond to low fluorescence bin), and
            to provide estimation of mutual information using compute_MI().

        load_unique_seqs_batches:
            Loads in unique sequences in data, to be used to check correlation
            between energy predictions and expected sign (e.g. negative repressor
            binding energy should correspond to low fluorescence bin), and
            to provide estimation of mutual information using compute_MI().

        compute_MI:
            Returns estimate of mutual information between matrix model
            predictions and observed fluoresence binning. Also returns the
            estimated joint probability distribution.

        Plotting styles and other plotting tools
        ----------------------
        set_plotting_style:
            Formats plotting enviroment to similar settings as used in Physical
            Biology of the Cell, 2nd edition. To format all plots within a
            script, execute `NB_sortseq_utils.set_plotting_style() in the
            preamble.

        plot_data:
            To be used by a Jupyter Notebook to easily plot expression shift,
            information footprint, and mutation rate for a single set of
            Sort-Seq experiments

        logo_on_matrix:
            To be used by a Jupyter Notebook to easily plot energy matrix and
            sequence logo for any inferred matrix.

License: MIT
    Copyright (c) 2017 Rob Phillips group @ California Institute of Technology

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.MIT

"""
import numpy as np
import pandas as pd
import scipy as sp
import scipy.ndimage
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as colors
from IPython.core.pylabtools import figsize
from Bio import SeqIO
from numpy import linalg as LA

#==============================================================================#
#        Sort-Seq processing (of .fasta or .fastq files)
#        -------------------------------
#==============================================================================#

def freqcount(records):
    """counts the frequency of each nucliotide at each position along the
    sequence considered.

    Parameter
    ---------
    records: this is an entry provided by the BioPython iterator,
        SeqIO.parse(seqname, seqtype) for a sequence file called seqname
        with sequence type equal to seqtype (e.g. .fastq).

    Returns
    -------
    countletter : float.
        Array of counts of A, C, G, Ts in record.seq, at each position of
        the sequence.
    """
    count = 0

    #check length of a sequence
    len_seq = 0
    for record in records:
        len_seq = len(record.seq)
        if len_seq != 0:
            break

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

    return countletter

#------------------------------------------------------------------------------#

def freqcount_cond(records, Let, position):
    """counts the frequency of each nucliotide at each position along the
    sequence considered. Only counts letters if sequence has letter 'Let'
    at specified 'position'.

    Parameter
    ---------
    records: this is an entry provided by the BioPython iterator,
        SeqIO.parse(seqname, seqtype) for a sequence file called seqname
        with sequence type equal to seqtype (e.g. .fastq).

    Let: letter ('A', 'C', 'G', 'T')

    position: position along sequence that contains Let

    Returns
    -------
    countletter : float.
        Array of counts of A, C, G, Ts in record.seq, at each position of
        the sequence.
    """
    count = 0

    #check length of a sequence
    len_seq = 0
    for record in records:
        len_seq = len(record.seq)
        if len_seq != 0:
            break

    countletter = np.zeros((4, len_seq))
    for record in records:
        count += 1
        if (record.seq[position] == Let):
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

    return countletter

#------------------------------------------------------------------------------#

def footprint_MI_calc(bin_freqs, finite_corr = True):
    """Calculates the information footprint mutual information values at
    each position along the sequence, using the nucleotide frequencies from each
    bin as calculated using freqcount().

    Parameter
    ---------
    bin_freqs: a numpy array (np.zeros([# bins, # letters (i.e. 4),
            length sequence]) that contains the letter frequences from each
            bin.

    finite_corr: specify True if finite data MI correction from Treves and
            Panzeri, 1995 is to be applied.

    Returns
    -------
    bin_freqs.sum(axis=1)[:,0]: total number of reads

    seqLength: length of sequence considered

    MI_mu: array 1*seqLength of mutual information values

    """

    # sequence length
    seqLength = len(bin_freqs[0][1,:])

    #MI calculation correction parameters
    nb = 4 #number of bases
    nmu = len(bin_freqs[:][:,1]) #number of bins
    #Calculate mutual information

    MI_mu = np.zeros(seqLength)
    bintotal = bin_freqs.sum(axis=1)[:,0].sum()


    #calculate MI (base; bin)
    for i in range(0,seqLength):
        # check that position has been mutated and skip position otherwise
        if np.any(bin_freqs[:,:,i] == 0):
            continue
        for j in range(0,nb):
            for k in range(0,nmu):
                Totalsum_bin = 0.0
                for x in range(0,nmu):
                    Totalsum_bin += np.sum(bin_freqs[x][j,i])
                MI_mu[i] += (bin_freqs[k][j,i]/bintotal)* \
                            np.log2(((bin_freqs[k][j,i])/bintotal))
                MI_mu[i] += -(bin_freqs[k][j,i]/bintotal)* \
                            np.log2(np.sum([bin_freqs[k][:,i]])/bintotal)
                MI_mu[i] += -(bin_freqs[k][j,i]/bintotal)* \
                            np.log2(Totalsum_bin/bintotal)

        # apply finite data correction if specified
        if finite_corr == True:
            MI_mu[i] += - ((nb-1)*(nmu-1)*np.log2(np.exp(1)))/(2*bintotal)

    return bin_freqs.sum(axis=1)[:,0], seqLength, MI_mu

def footprint_KLD_calc(bin_freqs1, bin_freqs2):
    """Calculates the Returns the Kullback-Leibler divergence values at
    each position along the sequence, using the nucleotide frequencies from each
    bin as calculated using freqcount() for two different datasets (e.g. two
    different growth conditions, or between wild-type and TF deletion strain).

    Parameter
    ---------
    bin_freqs1: a numpy array (np.zeros([# bins, # letters (i.e. 4),
            length sequence]) that contains the letter frequences from each
            bin.

    bin_freqs2: a numpy array (np.zeros([# bins, # letters (i.e. 4),
            length sequence]) that contains the letter frequences from each
            bin.

    Returns
    -------
    bin_freqs.sum(axis=1)[:,0]: total number of reads

    seqLength: length of sequence considered

    KLD_: array 1*seqLength of mutual information values

    """

    # sequence length
    seqLength = len(bin_freqs1[0][1,:])

    # #MI calculation correction parameters
    nb = 4 #number of bases
    nmu = len(bin_freqs1[:][:,1]) #number of bins

    #Calculate mutual information
    KL_calc = np.zeros(seqLength)
    bintotal1 = bin_freqs1.sum(axis=1)[:,0].sum()
    bintotal2 = bin_freqs2.sum(axis=1)[:,0].sum()

    #calculate KL
    # I will collapse 2d (bin,nucleotide) into single array
    # and calculate P(X) = P(bin,nucleotide)_dataset1 and
    # P(Y) = P(bin,nucleotide)_dataset2

    for i in range(0,seqLength):
        # check that position has been mutated and skip position otherwise
        if np.any(bin_freqs1[:,:,i] == 0):
            continue
        for j in range(0,nb):
            for k in range(0,nmu):
                Totalsum_bin1 = 0.0
                Totalsum_bin2 = 0.0
                for x in range(0,nmu):
                    Totalsum_bin1 += np.sum(bin_freqs1[x][j,i])
                    Totalsum_bin2 += np.sum(bin_freqs2[x][j,i])
                # P_X = (bin_freqs1[k][j,i]/bintotal1)
                P_X = (bin_freqs1[k][j,i]/Totalsum_bin1)

                # P_Y = (bin_freqs2[k][j,i]/bintotal2)
                P_Y = (bin_freqs2[k][j,i]/Totalsum_bin2)
                KL_calc[i] += -P_Y * np.log2(P_X/P_Y)


    return KL_calc

#------------------------------------------------------------------------------#

def calc_mutrate(bin_freqs):
    """Estimates the mutation rate at each position along the mutated sequence
    using the sorted sequences.

    Parameter
    ---------
    bin_freqs: a numpy array (np.zeros([# bins, # letters (i.e. 4),
            length sequence]) that contained the letter frequences from each
            bin.

    Returns
    -------

    mutrate: array 1*seqLength of mutation rate at each position in the library.

    """

    # sequence length
    seqLength = len(bin_freqs[0][1,:])

    #MI calculation correction parameters
    nb = 4 #number of bases
    nmu = len(bin_freqs[:][:,1]) #number of bins

    bintotal = bin_freqs.sum(axis=1)[:,0].sum()

    mutrate = np.zeros(seqLength)
    mutrate_bp = np.zeros([seqLength,4])

    for i in range(0,seqLength):
        for k in range(0,nmu):
            mutrate[i] += (np.sum([bin_freqs[k][:,i]]) - np.amax([bin_freqs[k][:,i]]))/bintotal
            for j in range(0,4):
                mutrate_bp[i,j] += bin_freqs[k][j,i]/bintotal
        for j in range(0,4):
            if mutrate_bp[i,j] >= 0.5:
                mutrate_bp[i,j] = 0
    return mutrate

#------------------------------------------------------------------------------#

def calc_deltabin(seq, files, bin_freqs, seqtype = "fastq"):
    """ At each position (starting at i), count number of sequences where region
    i is mutated. We are assuming that on average a mutation messes up binding,
    however this is not always the case. For example, especially with RNAP,
    there might be a couple positions that are not-at-all optimal for DNA
    binding.

    Parameter
    ---------
    seq: wild-type sequence of library region

    files: filenames (used to identify bin number, '...bin*.fastq')

    bin_freqs: numpy array (np.zeros([# bins, # letters (i.e. 4),
            length sequence]) that contained the letter frequences from each
            bin.

    seqtype: sequence file type (i.e. '.fastq' or '.fasta')

    Returns
    -------

    avgBin_counts: array 1*seqLength; contains counts used to calculate average
            of mutated nucleotides at each position.

    avgBin-avgbin_WT: average bin of mutated nucleotides at each position
            relative to wild-type average bin.

    avgbin_WT: average bin of wild-type nucleotides

    """

    seqLength = len(seq)

    avgBin_counts = np.zeros([len(files),seqLength])
    avgBin = np.zeros(seqLength)


    avgbin_WT = 0
    for j in range(0,len(files)):
        avgbin_WT += ( (j+1)*bin_freqs[j,:,0].sum() )/ bin_freqs[:,:,0].sum()

    for i in range(0,seqLength):

        for j, fname in enumerate(files):
            count = 0
            binnumber = int(fname[-7]) - 1

            for rec in SeqIO.parse(fname, seqtype):
                if (rec.seq[i:(i+1)] != seq[i:(i+1)]):
                    count += 1

            avgBin_counts[binnumber,i] = count


    for i in range(0,seqLength):
        for j in range(0,len(files)):
            avgBin[i] += ( (j+1)*avgBin_counts[j,i] )/avgBin_counts[:,i].sum()
    print(avgbin_WT)

    return avgBin_counts, (avgBin-avgbin_WT), avgbin_WT

#------------------------------------------------------------------------------#

def calc_deltabin_compare(seq, files, bin_freqs, seqtype = "fastq"):
    """ At each position (starting at i), count number of sequences where region
    i is mutated. We are assuming that on average a mutation messes up binding,
    however this is not always the case. For example, especially with RNAP,
    there might be a couple positions that are not-at-all optimal for DNA
    binding.

    Parameter
    ---------
    seq: wild-type sequence of library region

    files: filenames (used to identify bin number, '...bin*.fastq')

    bin_freqs: numpy array (np.zeros([# bins, # letters (i.e. 4),
            length sequence]) that contained the letter frequences from each
            bin.

    seqtype: sequence file type (i.e. '.fastq' or '.fasta')

    Returns
    -------

    avgBin_counts: array 1*seqLength; contains counts used to calculate average
            of mutated nucleotides at each position.

    avgBin-avgbin_WT: average bin of mutated nucleotides at each position
            relative to wild-type average bin.

    avgbin_WT: average bin of wild-type nucleotides

    """

    seqLength = len(seq)

    avgBin_counts = np.zeros([4,len(files),seqLength])
    avgBin = np.zeros(seqLength)

    avgBin_shift_param = np.zeros([4,seqLength])

    avgbin_WT = 0
    for j in range(0,len(files)):
        avgbin_WT += ( (j+1)*bin_freqs[j,:,0].sum() )/ bin_freqs[:,:,0].sum()

    # now, for each position, and each bp (except WT),
    # compare the average bin for that nucleotide relative
    # to the average above. If it is +, give value of +1;
    # if it is a negative difference, asign value of -1.
    # asign the WT bp as zero (might reconsider this step)

    for i in range(0,seqLength):
        for j, fname in enumerate(files):
            count = 0
            binnumber = int(fname[-7]) - 1

            for rec in SeqIO.parse(fname, seqtype):
                if (rec.seq[i:(i+1)] != seq[i:(i+1)]):
                    count += 1

            avgBin_counts[binnumber,i] = count


    for i in range(0,seqLength):
        for j in range(0,len(files)):
            avgBin[i] += ( (j+1)*avgBin_counts[j,i] )/avgBin_counts[:,i].sum()
    print(avgbin_WT)

    return avgBin_counts, (avgBin-avgbin_WT), avgbin_WT

#------------------------------------------------------------------------------#

def calc_deltabin_3bpavg(seq, files, bin_freqs, seqtype = "fastq"):
    """
    At each position (starting at i), count number of sequences where
    region (i):(i+3) is mutated. This is sort of a rolling average and not critical
    to the result. It just ends up a bit cleaner than if we looked at a single
    base pair since. We are assuming that on average a mutation messes up binding,
    however this is not always the case. For example, especially with RNAP, there might
    be a couple positions that are not-at-all optimal for DNA binding.

    Parameter
    ---------
    seq: wild-type sequence of library region

    files: filenames (used to identify bin number, '...bin*.fastq')

    bin_freqs: numpy array (np.zeros([# bins, # letters (i.e. 4),
            length sequence]) that contained the letter frequences from each
            bin.

    seqtype: sequence file type (i.e. '.fastq' or '.fasta')

    Returns
    -------

    avgBin_counts: array 1*seqLength; contains counts used to calculate average
            of mutated nucleotides at each position.

    avgBin-avgbin_WT: average bin of mutated nucleotides at each position
            relative to wild-type average bin.

    """
    seqLength = len(seq)

    avgBin_counts = np.zeros([len(files),seqLength])
    avgBin = np.zeros(seqLength)

    #filecount = 0

    avgbin_WT = 0
    for j in range(0,len(files)):
        avgbin_WT += ( (j+1)*bin_freqs[j,:,0].sum() )/ bin_freqs[:,:,0].sum()
    print('average_bin_WT', avgbin_WT)

    for i in range(0,seqLength-2):

        for j, fname in enumerate(files):
            count = 0
            binnumber = int(fname[-7]) - 1

            for rec in SeqIO.parse(fname, seqtype):
                if (rec.seq[i:(i+2)] != seq[i:(i+2)]):
                    count += 1

            avgBin_counts[binnumber,i] = count

    for i in range(0,seqLength-2):
        for j in range(0,len(files)):
            avgBin[i] += ( (j+1)*avgBin_counts[j,i] )/avgBin_counts[:,i].sum()


    return avgBin_counts, (avgBin-avgbin_WT)

#------------------------------------------------------------------------------#
#
# def genomescan(record, EnergyMat):
#     """counts the frequency of each nucliotide at each position along
#     the seuqence considered
#
#     Parameter
#     ---------
#
#
#     Returns
#     -------
#
#     """
#
#     GenomeMatscan = 0*EnergyMat
#     for i in range(0,len(record)):
#         if record[i] == 'A':
#             GenomeMatscan[0,i]= 1
#         elif record[i] == 'C':
#             GenomeMatscan[1,i]= 1;
#         elif record[i] == 'G':
#             GenomeMatscan[2,i]= 1;
#         elif record[i] == 'T':
#             GenomeMatscan[3,i]= 1;
#         else:
#             print("error: unexpected letter")
#     return GenomeMatscan
#



#==============================================================================#
#        Energy matrix and sequence logo tools
#        -------------------------------
#==============================================================================#

def seq2mat(seq):
    """
    Takes input nucleotide sequence and returns in matrix representation.

    Parameter
    ---------
    seq: sequence string of A,C,G,Ts

    Returns
    -------
    mat:
        4xlen(seq) array of ones and zeros (WT position = 1, zero otherwise)

    """
    # # check that emat matrix is 4xL
    # if matrix_model.shape[0] != 4:
    #     emat = emat.T

    seq_dict = {'A':0,'C':1,'G':2,'T':3}
    mat = sp.zeros((4,len(seq)),dtype=int)
    for i,bp in enumerate(seq):
        mat[seq_dict[bp],i] = 1
    return mat

#------------------------------------------------------------------------------#

def zero_matrix(emat):
    """
    Takes in matrix model and set the smallest element of each column (for
    each position) to zero energy.

    Parameter
    ---------
    emat: 4xlen(binding site) array of matrix model

    Returns
    -------
    emat:
        Returns 4xlen(binding site) array of matrix model in which each column
        has been shifted such that the smallest element has zero energy.

    """
    # check that emat matrix is 4xL
    if emat.shape[0] != 4:
        emat = emat.T

    for j in range(emat.shape[1]):
        emat[:,j] = emat[:,j] - emat[:,j].min()

    return emat


def zero_matrix_WT(emat, seq):
    """
    Takes in matrix model and set the smallest element of each column (for
    each position) to zero energy.

    Parameter
    ---------
    emat: 4xlen(binding site) array of matrix model

    seq: wild-type sequence string of A,C,G,Ts

    Returns
    -------
    emat:
        Returns 4xlen(binding site) array of matrix model in which each column
        has been shifted such that the smallest element has zero energy.

    """
    # check that emat matrix is 4xL
    if emat.shape[0] != 4:
        emat = emat.T

    # make matrix of WT sequence
    seqmat = seq2mat(seq)
    # set the smallest element of each column to zero
    for j in range(emat.shape[1]):
        emat[:,j] = emat[:,j] - np.dot(emat[:,j],seqmat[:,j])

    return emat

def compute_energies(seqs,batches,emat):
    """seqs: matrix of sequences, should be 4xLxN
    batches: vector of batches
    emat: energy matrix, 4xL"""
    dot = emat[:,:,sp.newaxis]*seqs
    energies = dot.sum(0).sum(0)
    return energies

# def genomescan(seqs,emat):
#     """seqs: matrix of sequences, should be 4xLxN
#     batches: vector of batches
#     emat: energy matrix, 4xL"""
#     dot = emat[:,:,sp.newaxis]*seqs
#     energies = dot.sum(0).sum(0)
#     return energies

#------------------------------------------------------------------------------#


def fix_matrix_gauge(emat):

    """Fix gauge of an energy matrix such that the average value
    of each column is zero (columns correspond to positions), and
    overall matrix norm is equal to 1."""
    # fix mean
    for j in range(emat.shape[1]):
        emat[:,j] = emat[:,j] -sp.mean(emat[:,j])
    # fix inner product
    # emat = emat/sp.sqrt(sp.sum(emat*emat))
    emat = emat / LA.norm(emat)

    return emat


#------------------------------------------------------------------------------#

def estimate_scalefactor(emat):
    """This function estimates the expected scale factor necessary to
    convert an energy weight matrix that is in arbitrary units. for
    a matrix entry, \theta_ij, and with kBT units, A*\theta_ij + B,
    this provides  the A factor when the energy matrix has average value
    equal to zero and overall matrix norm is equal to 1. We will assume
    that the average change in binding energy for a single point mutation
    is equal to 2 kBT and determine the scale factor A to satisfy this.


    """
    # I'm going to calculate the average shift in binding energy for every
    # single point mutation and then say this should be equal to 2-3 KbT.

    # check that emat matrix is 4xL
    if emat.shape[0] != 4:
        emat = emat.T

    # # re-fix gauge just to be safe
    # emat = fix_matrix(emat)

    # set the smallest element of each column to zero
    for j in range(emat.shape[1]):
        emat[:,j] = emat[:,j] - emat[:,j].min()

    # aim for an average mutation shift of 2 kBT
    return 2/(emat.sum()/(4*emat.shape[1]))

#------------------------------------------------------------------------------#


def emattoIUPAC(emat, scalefactor = 15):
    """
    This function will take an energy matrix dataframe and
    return the IUPAC sequence.
    """

    seq = {0:'A', 1:'C', 2:'G', 3:'T'}

    # # Load in pandas dataframe energy matrix
    # emat_RNAP = pd.read_csv(ematdir)

    # make approx conversion to kBT with scale factor equal to 3.5
    emat_= np.array(scalefactor * emat[['A','C','G','T']].T)

    # Exponentiate energy parameters
    emat_exp = np.exp(-emat_)

    # Determine p_ij
    p_ij_emat = emat_exp / emat_exp.sum(axis=0)

    IUPAC = ''
    for i in range(0,len(emat_exp.sum(axis=0))):
        if p_ij_emat[0,i] >= 0.6:
            IUPAC += seq[0]
        elif p_ij_emat[1,i] >= 0.6:
            IUPAC += seq[1]
        elif p_ij_emat[2,i] >= 0.6:
            IUPAC += seq[2]
        elif p_ij_emat[3,i] >= 0.6:
            IUPAC += seq[3]
            continue
        elif p_ij_emat[0,i]+ p_ij_emat[1,i] >= 0.7:
            IUPAC += 'M'
            continue
        elif p_ij_emat[0,i]+ p_ij_emat[2,i] >= 0.7:
            IUPAC += 'R'
            continue
        elif p_ij_emat[0,i]+ p_ij_emat[3,i] >= 0.7:
            IUPAC += 'W'
            continue
        elif p_ij_emat[1,i]+ p_ij_emat[2,i] >= 0.7:
            IUPAC += 'S'
            continue
        elif p_ij_emat[1,i]+ p_ij_emat[3,i] >= 0.7:
            IUPAC += 'Y'
            continue
        elif p_ij_emat[2,i]+ p_ij_emat[3,i] >= 0.7:
            IUPAC += 'K'
            continue
        elif p_ij_emat[1,i]+ p_ij_emat[2,i] + p_ij_emat[3,i] >= 0.9:
            IUPAC += 'B'
            continue
        elif p_ij_emat[0,i]+ p_ij_emat[2,i] + p_ij_emat[3,i] >= 0.9:
            IUPAC += 'D'
            continue
        elif p_ij_emat[0,i]+ p_ij_emat[1,i] + p_ij_emat[3,i] >= 0.9:
            IUPAC += 'H'
            continue
        elif p_ij_emat[0,i]+ p_ij_emat[1,i] + p_ij_emat[2,i] >= 0.9:
            IUPAC += 'V'
            continue

        else:
            IUPAC += 'N'

    return IUPAC


#==============================================================================#


def hpd(trace, mass_frac):
    """
    Returns highest probability density region given by
    a set of samples.
    Parameters
    ----------
    trace : array
        1D array of MCMC samples for a single variable
    mass_frac : float with 0 < mass_frac <= 1
        The fraction of the probability to be included in
        the HPD.  For hreple, `massfrac` = 0.95 gives a
        95% HPD.

    Returns
    -------
    output : array, shape (2,)
        The bounds of the HPD

    Notes
    -----
    We thank Justin Bois (BBE, Caltech) for developing this function.
    http://bebi103.caltech.edu/2015/tutorials/l06_credible_regions.html
    """
    # Get sorted list
    d = np.sort(np.copy(trace))

    # Number of total samples taken
    n = len(trace)

    # Get number of samples that should be included in HPD
    n_samples = np.floor(mass_frac * n).astype(int)

    # Get width (in units of data) of all intervals with n_samples samples
    int_width = d[n_samples:] - d[:n - n_samples]

    # Pick out minimal interval
    min_int = np.argmin(int_width)

    # Return interval
    return np.array([d[min_int], d[min_int + n_samples]])



#==============================================================================#
#        MI tools for MCMC matrix processing
#        -------------------------------
#==============================================================================#

def load_unique_seqs_batches(data_fn,mut_region_start,mut_region_length):
    """Load in unique sequence-batche pairs from data file.

    Parameter
    ---------
    data_fn: csv file containing sequence, batch
    mut_region_start: sequence index corresponding to start of ROI
    mut_region_length: self-evident

    Returns
    -------
    seq_mat: 4xmut_region_lengthxN matrix containing sequence information
    batch_vec: N length vector containing batches
    """


    f = open(data_fn)
    # read lines into one big list and transform into a set. This
    # automatically gets rid of duplicates
    # lines with region of interest selected
    roi_list = [(line.split(',')[0][mut_region_start:mut_region_start+mut_region_length],
            line.split(',')[1].strip()) for line in f if line.strip()]
    f.close()

    lines_unique = list(set(roi_list))
    N = len(lines_unique)

    # instantiate batch and sequence variables
    batch_vec = sp.zeros(N,dtype=int)
    seq_mat = sp.zeros((4,mut_region_length,N),dtype=int)

    for i, line in enumerate(lines_unique):
        batch_vec[i] = int(line[1])
        seq_mat[:,:,i] = seq2mat(line[0])
    batch_vec = batch_vec-batch_vec.min()

    return seq_mat,batch_vec


def load_seqs_batches(data_fn,mut_region_start,mut_region_length):
    """Load in unique sequence-batche pairs from data file.

    Parameter
    ---------
    data_fn: csv file containing sequence, batch
    mut_region_start: sequence index corresponding to start of ROI
    mut_region_length: self-evident

    Returns
    -------
    seq_mat: 4xmut_region_lengthxN matrix containing sequence information
    batch_vec: N length vector containing batches
    """

    N = 0
    f = open(data_fn)
    for line in f:
        if line.strip():
            N = N + 1
    f.close()
    print(N)

    # instantiate batch and sequence variables
    batch_vec = sp.zeros(N,dtype=int)
    seq_mat = sp.zeros((4,mut_region_length,N),dtype=int)

    f = open(data_fn)
    for i, line in enumerate(f):
        if line.strip():
            sb = line.split(',')
            batch_vec[i] = int(sb[1])
            seq_mat[:,:,i] = seq2mat(sb[0][mut_region_start:mut_region_start+mut_region_length])
    f.close()
    batch_vec = batch_vec-batch_vec.min()

    return seq_mat,batch_vec

#------------------------------------------------------------------------------#

def compute_MI(seqs,batches,emat):
    '''
    Note that this is not exactly the same code used to calculate MI
    that is used in MPAthic, but provide a similar estimate and calculation
    of the joint distribution that is useful for visualizing the result.

    Parameters
    -------

    Returns
    -------

    '''

    # preliminaries
    n_seqs = len(batches)
    n_batches = batches.max() + 1 # assumes zero indexed batches
    n_bins = 1000

    #energies = sp.zeros(n_seqs)
    f = sp.zeros((n_batches,n_seqs))

    # compute energies
    #for i in range(n_seqs):
    #    energies[i] = sp.sum(seqs[:,:,i]*emat)
    # alternate way
    dot = emat[:,:,sp.newaxis]*seqs
    energies = dot.sum(0).sum(0)


    # sort energies
    inds = sp.argsort(energies)
    for i,ind in enumerate(inds):
        f[batches[ind],i] = 1.0/n_seqs # batches aren't zero indexed


    # bin and convolve with Gaussian
    f_binned = sp.zeros((n_batches,n_bins))

    for i in range(n_batches):
        f_binned[i,:] = sp.histogram(f[i,:].nonzero()[0],bins=n_bins,range=(0,n_seqs))[0]
    #f_binned = f_binned/f_binned.sum()
    f_reg = scipy.ndimage.gaussian_filter1d(f_binned,0.04*n_bins,axis=1)
    f_reg = f_reg/f_reg.sum()

    # compute marginal probabilities
    p_b = sp.sum(f_reg,axis=1)
    p_s = sp.sum(f_reg,axis=0)

    # finally sum to compute the MI
    MI = 0
    for i in range(n_batches):
        for j in range(n_bins):
            if f_reg[i,j] != 0:
                MI = MI + f_reg[i,j]*sp.log2(f_reg[i,j]/(p_b[i]*p_s[j]))
    print(MI)
    return MI,f_reg




#==============================================================================#
#       Plotting Configuration
#        -------------------------------
#==============================================================================#

def cm2inch(*tupl):
    """Converts inches into mm for figures.
    """

    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def set_plotting_style1():
    """
    Formats plotting enviroment to that used in Physical Biology of the Cell,
    2nd edition. To format all plots within a script, simply execute
    `plotting_defaults.set_plotting_style() in the preamble.

    Parameters
    -------

    Returns
    -------

    """
    rc = {'lines.linewidth': 2,
          'axes.labelsize': 22,
          'axes.titlesize': 22,
          'axes.facecolor': '#E3DCD0',
          'xtick.major' : 22,
          'ytick.major' : 22,
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': ':',
          'grid.linewidth': 1.5,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.fontsize': 22}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)
    sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)


def set_plotting_style2():
    """
    Formats plotting enviroment to that used in Physical Biology of the Cell,
    2nd edition. To format all plots within a script, simply execute
    `plotting_defaults.set_plotting_style() in the preamble.
    """
    rc = {'lines.linewidth': 2,
          'axes.labelsize': 20,
          'axes.titlesize': 20,
          'axes.facecolor': '#E3DCD0',
          'xtick.major' : 20,
          'xtick.labelsize': 14,
          'ytick.labelsize': 20,
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': ':',
          'grid.linewidth': 1.5,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.fontsize': 13}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    #sns.set_style('darkgrid', rc=rc)
    sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)

def set_plotting_style3():
    """
    for mass spec plots: Formats plotting enviroment to that used in
    Physical Biology of the Cell,
    2nd edition. To format all plots within a script, simply execute
    `plotting_defaults.set_plotting_style() in the preamble.
    """
    rc = {'lines.linewidth': 2,
          'axes.labelsize': 18,
          'axes.titlesize': 20,
          'axes.facecolor': '#E3DCD0',
          'xtick.major' : 13,
          'xtick.labelsize': 'large',
          'ytick.labelsize': 13,
          'ytick.linewidth': 1.5,
          'ytick.color': '#ffffff',
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': ':',
          'grid.linewidth': 1.5,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.fontsize': 13}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)
    sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)


def set_plotting_style4():
    """
    Formats plotting enviroment to that used in Physical Biology of the Cell,
    2nd edition. To format all plots within a script, simply execute
    `plotting_defaults.set_plotting_style() in the preamble.
    """
    rc = {'lines.linewidth': 0.3,
          'axes.labelsize': 8,
          'axes.titlesize': 8,
          'axes.facecolor': '#E3DCD0',
          'xtick.major' : 8,
          'xtick.labelsize': 6,
          'ytick.labelsize': 8,
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': ':',
          'grid.linewidth': 0.3,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.fontsize': 13}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)
    sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)


def set_plotting_style_emat():
    """
    Formats plotting enviroment to that used in Physical Biology of the Cell,
    2nd edition. To format all plots within a script, simply execute
    `plotting_defaults.set_plotting_style() in the preamble.
    """
    rc = {'lines.linewidth': 0.3,
          'axes.labelsize': 8,
          'axes.titlesize': 8,
          'axes.facecolor': '#E3DCD0',
          'xtick.major' : 8,
          'xtick.labelsize': 6,
          'ytick.labelsize': 8,
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': ':',
          'grid.linewidth': 0.3,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.fontsize': 13}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    # sns.set_style('darkgrid', rc=rc)
    sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)


def set_plotting_style_MS():
    """
    Formats plotting enviroment to that used in Physical Biology of the Cell,
    2nd edition. To format all plots within a script, simply execute
    `plotting_defaults.set_plotting_style() in the preamble.
    """
    rc = {'lines.linewidth': 0.8,
          'axes.labelsize': 8,
          'axes.titlesize': 8,
          'axes.facecolor': '#E3DCD0',
          'xtick.major' : 8,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,
          'ytick.minor.linewidth': 0.3,
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': ':',
          'grid.linewidth': 0.5,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.fontsize': 13}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)
    sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)

def set_plotting_style_scan():
    """
    Formats plotting enviroment to that used in Physical Biology of the Cell,
    2nd edition. To format all plots within a script, simply execute
    `plotting_defaults.set_plotting_style() in the preamble.
    """
    rc = {'lines.linewidth': 0.8,
          'axes.labelsize': 8,
          'axes.titlesize': 8,
          'axes.facecolor': '#E3DCD0',
          'xtick.major' : 8,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,
          'ytick.minor.linewidth': 0.3,
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': ':',
          'grid.linewidth': 0.5,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.fontsize': 13}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    #sns.set_style('darkgrid', rc=rc)
    sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)

def plot_data(df, plottype = 'expshift'):

    if (len(df.promoter.unique())!=1) or (len(df.strain.unique())!=1) or (len(df.media.unique())!=1):
        raise RuntimeError('Looks like there is data from multiple promoters or different strains.')

    colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]

    if plottype == 'expshift':

        fig = plt.figure(figsize=(12,3.5))
        ax1 = fig.add_subplot(111)

        seqLength = len(df.position.unique())
        ind = df.groupby(['position'])['position'].mean().values

        ax1.bar(ind, df.groupby(['position'])['expshift'].mean(), width=0.8, \
                linewidth = 0, color = 'k')
        ax1.set_xlabel('position')
        ax1.grid(b=False)
        ax1.set_facecolor('white')

        ax1.set_xlim([ind[0]-.5, ind[-1]+.5])

        # Annotate y axis
        ylim = ax1.get_ylim()
        yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
        ax1.set_yticks(yticks)
        ax1.set_yticklabels(['-','+'], fontname='Arial')
        ax1.set_ylabel('expression\nshift')

        # Make grid lines
        xticks = [x for x in ind if x%20 == 0]
        for n, x in enumerate(xticks):
                ax1.axvline(x, color='lightgray', linewidth=0.36, zorder=-1)

    if plottype == 'infofootprint':
        fig = plt.figure(figsize=(12,3.5))
        ax1 = fig.add_subplot(111)

        ind = df.position.unique()
        ax1.bar(ind, df.groupby(['position'])['MI'].mean(), width=0.8, \
                linewidth = 0,color=colours[2])
        ax1.set_xlim(ind[0],ind[-1])
        ax1.set_ylabel('mutual information\n(bits)')
        ax1.set_xlabel(r'position (relative to start of $' + df.promoter.unique() + '$ start codon)')
        ax1.grid(b=False)

    if plottype == 'mutrate':

        fig = plt.figure(figsize=(12,3.5))
        ax1 = plt.subplot(111)

        ind = df.position.unique()
        ax1.bar(ind, df.groupby(['position'])['mutation_rate'].mean(), width=0.8, \
                linewidth = 0)

        ax1.set_xlim(ind[0],ind[-1])
        ax1.set_ylabel('mutation rate')
        ax1.set_xlabel(r'position (relative to start of $' + df.promoter.unique() + '$ start codon)')



# Plot logo on top of energy matrix
def logo_on_matrix(ax, energy_df, relative_scale=1, relative_spacing=.1,
                   fontsize=9, show_positions=False, wt_seq=None, acgt_pad = 7):
    import anylogo
    if wt_seq==None:
        wt_seq = ''.join(energy_df.sort_values(by='position')['WT_sequence'].values)


    energy_df = energy_df[['A','C','G','T']].copy()
    energy_df_scaled = estimate_scalefactor(np.array(energy_df))*energy_df.copy()
    energy_df_scaled = energy_df_scaled[['A','C','G','T']]
    # Create background array
    gc = .508
    background_array =pd.DataFrame( [[(1-gc)/2,gc/2,gc/2,(1-gc)/2]])
    # create background nucleotide frequencies dataframe
    background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df), 1)), columns=['A','C','G','T'])

    # Set color scale - I want the colorbar to be symmetric and will pick values#
    # that seem appropriate for all matrices.

    emat_min=-np.max([-energy_df.min().min(),energy_df.max().max()])#-0.4
    print(emat_min)
    emat_max=np.max([-energy_df.min().min(),energy_df.max().max()])#0.4
    print(emat_max)
    mid_val=0.0

    emat_ymin = -2 * (relative_scale + relative_spacing)
    emat_ymax = -2 * relative_spacing
    yticks = np.linspace(emat_ymin, emat_ymax, 9)[[1, 3, 5, 7]]
    yticklabels = list('TGCA')
    cmap = plt.get_cmap('RdBu_r')
    anylogo.draw(ax, effect_df=energy_df_scaled, logo_type='information', background = background_df,
            find_beta=False)
    L = len(energy_df)
    # im = ax.imshow(utils.zero_matrix_WT(energy_df.values.T,wt_seq), aspect='auto',
    #           extent=(-.5, L - .5, emat_ymin, emat_ymax), zorder=100, cmap=cmap)
    im = ax.imshow(zero_matrix_WT(energy_df.values.T,wt_seq),
                interpolation='none',
                cmap='RdBu_r',
                clim=(emat_min, emat_max),
                norm = MidpointNormalize(midpoint = mid_val,
                        vmin = emat_min, vmax = emat_max),
                extent=(-.5, L - .5, emat_ymin, emat_ymax),
                zorder=100,
                aspect='auto')
    ax.set_ylim([emat_ymin, 2])
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=fontsize, horizontalalignment='center')
    ax.set_ylabel('')
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_tick_params(length=0)
    if not show_positions:
        ax.set_xticks([])
    y = .5*emat_ymax
    if not wt_seq is None:
        assert len(wt_seq) == L, \
            'Error! len(wt_seq)=%d does not match len(energy_df)=%d.' % (len(wt_seq), L)
        for i in range(L):
            ax.text(i, y, wt_seq[i], horizontalalignment='center', verticalalignment='center')
    ax.tick_params(axis='y', pad=acgt_pad)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # create an axes on the right side of ax. The width of cax will be 3%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)

    cbar = plt.colorbar(im, cax=cax, ticks=[emat_min, 0, emat_max], label='energy (a.u.)')
    cbar.ax.set_yticklabels([('%.2f' % emat_min), '0', ('%.2f' % emat_max)], fontname='Arial')
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(axis=u'both', which=u'both',length=0)


# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
     http://chris35wills.github.io/matplotlib_diverging_colorbar/
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
