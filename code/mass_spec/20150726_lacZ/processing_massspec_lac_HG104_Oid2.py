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

promoter = 'lacZ_Oid'
date = 20150726
WT_sequence = 'CAAAACTGCAG CTCGTATGTTGTGTGG AATTGTGAGCGCTCACAATT TCACACAGGAAA'
control_sequence = 'CAAAACTGCAG CTCGTATGTTGTGTGG AATCCCCCTAGGTGAACCATT TCACACAGGAAA'
strain = 'HG104deltalysA'
media = 'M9glucose'
binding_additions = 'na'
replicate = 2
spec_machine = 'orbitrap_classic'
labeling = 'reverse'
df_details = pd.DataFrame([[promoter, date, WT_sequence, control_sequence,
strain, media, binding_additions, replicate, spec_machine, labeling]])

#==============================================================================
# Load in the data files
#==============================================================================
# Load the output of protein groups containing SILAC ratios into a DataFrame
df = pd.read_csv('../../data/mass_spec/20150726_HG104/proteinGroups.txt',
                 delimiter='	')
# remove non-e.coli items detected
df = df[~df['Protein IDs'].str.contains("REV")]
df = df[~df['Protein IDs'].str.contains("CON")]
# Ratio H/L variability [%] = standard deviation of the natural logs of the constituent peptide ratios multiplied by 100.
df['std_dev'] = np.exp((df['Ratio H/L variability [%] HG104_Oid_L_H'] / 100))

# from MaxQuant output, I want to collect the following for each replicate:
df = df[['Protein IDs','Gene names','Peptides', 'Intensity HG104_Oid_L_H',
         'Intensity L HG104_Oid_L_H', 'Intensity H HG104_Oid_L_H',
         'Ratio H/L type', 'Ratio H/L HG104_Oid_L_H', 'std_dev']]

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
              'MaxQuant_ratio_std_dev', 'TF_check',
              'promoter', 'date', 'wt_sequence', 'control_sequence', 'strain',
              'media', 'binding_additions', 'replicate', 'spec_machine',
              'labeling']

#==============================================================================
# Renormalization of SILAC ratio data correction.
# Due to potential experimental variability, the log ratios
# are median shifted so that the entire sample population is centered
# at 1:1 (heavy: light).
#==============================================================================

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

df.to_csv(+ str(date) + '_' + promoter + '_' + str(replicate) +
        '_mass_spec_summary.csv', index=False)
