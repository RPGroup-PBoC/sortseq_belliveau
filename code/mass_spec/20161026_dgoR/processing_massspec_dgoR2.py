# Purpose of script:
# For each DNA affinity purification dataset, this script will extract the
# relavant information from the 'proteingroups.txt' output after running the
# mass spec .raw data on MaxQuant (version 1.5.XX). It will also add in
# the additional experimental details (DNA sequence used, strain used for lysate
# media cells were grown in, etc. that are noted below). Lastely, it will
# calculate median shifted enrichment ratios using all available ratio
# measurements.

# Import dependencies.
import sys
import glob
import numpy as np
import pandas as pd

promoter = 'dgoR'
date = 20161026
WT_sequence = 'CAAAACTGCAG TTTTATGCATTGTTC TTTTTGTGATCTAAATTGTAGTACAACAATATAAGTTTGTACTAC ATTACACGCACGGCA'
control_sequence = 'CAAAACTGCAG TTTTATGCATTGTTC TTAACAATTTAAAAACGTTAAGGTACATACGGTGATGACTTCACC ATTACACGCACGGCA'
strain = 'MG1655deltalysA'
media = 'M9glucose'
binding_additions = 'na'
replicate = 2
spec_machine = 'orbitrap_elite'
labeling = 'forward'
df_details = pd.DataFrame([[promoter, date, WT_sequence, control_sequence,
strain, media, binding_additions, replicate, spec_machine, labeling]])

#==============================================================================
# Load in the data files
#==============================================================================
# Load the output of protein groups containing SILAC ratios into a DataFrame
df = pd.read_csv('analysis/20161115_DNAAP_dgoR_three_trials/data/NBM1001_dgoR_scramble_proteinGroups_dgoR_threetrials.txt',
                 delimiter='	')
# remove non-e.coli items detected
df = df[~df['Protein IDs'].str.contains("REV")]
df = df[~df['Protein IDs'].str.contains("CON")]

# from MaxQuant output, I want to collect the following for each replicate:
df = df[['Protein IDs','Gene names','Peptides', 'Intensity 2', 'Intensity L 2',
         'Intensity H 2', 'Ratio H/L type', 'Ratio H/L 2']]
# Load the list of proteins with predicted or confirmed DNA binding motifs
# in E. coli
EcoliTF = pd.read_csv('../20150712_EcoCyc_TF_query-results.txt',delimiter='	')

# compare list of proteins with DNA binding motifs to those detected
df_temp = np.zeros(len(df))
for i, gene in enumerate(df['Gene names']):
    test = pd.Series(gene).isin(EcoliTF['Protein'])

    if test[0] == True:
        df_temp[i] = 1

# add to dataframe
df['TF_check'] = df_temp

df = df.reset_index(drop=True)

# Tile experimental details for number of protein entries in df
# i.e. make an array containing the same experimental details in
# each row, which will be merged with the experimental results.
df_details = pd.DataFrame(pd.np.tile(df_details, (len(df), 1)))

#==============================================================================
# Combine experimental data with experimental details and rename columns.
#==============================================================================

# combine mass spec data with experimental details into tidy format
df = pd.concat([df, df_details], axis=1, ignore_index=True)

# Rename the columns of the data_frame
df.columns = ['ProteinIDs','gene','peptides', 'intensity_total', 'intensity_L',
              'intensity_H',  'maxquant_ratio_type', 'maxquant_ratio',
              'TF_check', 'promoter', 'date',
              'wt_sequence', 'control_sequence', 'strain', 'media',
              'binding_additions', 'replicate', 'spec_machine', 'labeling']

#==============================================================================
# Renormalization of SILAC ratio data correction.
# Due to potential experimental variability, the log ratios
# are median shifted so that the entire sample population is centered
# at 1:1 (heavy: light).
#==============================================================================

# Note: Following initial analysis, it was found that that DgoR was the TF with
# greatest enrichment by both intensity and the 'MaxQuant' reported ratios.
# However, there appeared to be quite a bit of difference, with the 'MaxQuant'
# ratio shifted significantly down. The data was re-analyzed using MaxQuant
# without 'requantify' and 'match between runs' turned on. This analysis did not
# identify any 'light' DgoR, suggesting there may be limited accuracy in the
# magnitude of the enirchment ratio and it may be larger than reported. The
# enirchment value based on the intensity values was used.

# df['maxquant_ratio'][df.gene == 'dgoR'] = \
#     df['intensity_H'][df.gene == 'dgoR']/ df['intensity_L'][df.gene == 'dgoR']

# Let's median shift the maxquant_ratio
df['maxquant_logratio'] = np.log(df['maxquant_ratio'])
# Calculate the median of log ratios
median_shift = df['maxquant_logratio'].replace([-np.inf,np.inf], np.nan).dropna().median()
# Calculate the median shifted log ratios and non-log ratios for all proteins
# in the sample.
df['maxquant_logratio_medianshift'] = df['maxquant_logratio'] - median_shift
df['maxquant_ratio_medianshift'] = np.exp(df['maxquant_logratio_medianshift'])

#==============================================================================
# Generate the summary file
#==============================================================================

df.to_csv('output/' + str(date) + '_' + promoter + '_' + str(replicate) +
        '_mass_spec_summary.csv', index=False)
