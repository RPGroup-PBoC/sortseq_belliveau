# Purpose of script:
# Grab the proteinGroups txt file and add in column that identifies whether
# a gene/protein has expected DNA binding domain.
# remove xylE 'short DNA' columns since data was not used in paper.

# Import dependencies.
import sys
import glob
import numpy as np
import pandas as pd

#==============================================================================
# Load in the data files
#==============================================================================
# Load the output of protein groups containing SILAC ratios into a DataFrame
df = pd.read_csv('../../data/mass_spec/20170925_all_requantify_ON/proteinGroups.txt',
                 delimiter='	')
# # remove non-e.coli items detected
# df = df[~df['Protein IDs'].str.contains("REV")]
# df = df[~df['Protein IDs'].str.contains("CON")]
#
# # from MaxQuant output, I want to collect the following for each replicate:
# df = df[['Protein IDs','Gene names','Peptides', 'Intensity 1', 'Intensity L 1',
#          'Intensity H 1', 'Ratio H/L type', 'Ratio H/L 1']]

# Load the list of proteins with predicted or confirmed DNA binding motifs
# in E. coli
EcoliTF = pd.read_csv('../20150712_EcoCyc_TF_query-results.txt',delimiter='	')

# compare list of proteins with DNA binding motifs to those detected
df_schmidt['TF_check'] = np.where(df_schmidt.Gene.isin(EcoliTF['Protein']), 'yes', 'no')

# df['TF_check'] = ''
# # compare list of proteins with DNA binding motifs to those detected
# # df_temp = np.zeros(len(df))
# for i, gene in enumerate(df['Gene names']):
#     # test = pd.Series(gene).isin(EcoliTF['Protein'])
#     print(pd.Series(gene).isin(EcoliTF['Protein'])[0] )
#     print(gene)
#     if pd.Series(gene).isin(EcoliTF['Protein'])[0] == True:
#         df['TF_check'][df['Gene names']==gene] = 'yes'
#     else:
#         df['TF_check'][df['Gene names']==gene] = 'no'


#==============================================================================
# Generate the summary file
#==============================================================================

df.to_csv('output/20170929_proteinGroups_all.txt', index=False,sep='	')
