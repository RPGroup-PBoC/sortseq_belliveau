# Purpose of script:
# For each flow cytometry dataset, I have several files that I want to use to
# generate compiled tidy datasets. This includes: 1) DATE_hist.txt -> histogram
# generated in the Summit flow cytometry software, 2) DATE-stats.csv ->
# statistics based on this histogram, 3) notes from the measurements that day,
# and 4) a tidy csv file that contains the details related to the promoter,
# growth conditions, flow cytometry measurement FL1 voltage and gain, strain,
# etc. The first file should combine the histogram data with the measurement
# notes and the details in the tidy csv. The second file should combine the
# stats with  the measurement notes and the details in the tidy csv.

# Import dependencies.
import os
import sys
import itertools
import glob
import re
import numpy as np
import pandas as pd

# Set the plotting style.
import sys
sys.path.insert(0, '../../')

promoter = 'titration'
date = 20170506
WT_sequence = 'na'
control_sequence = 'na'
strain = 'MG1655deltalysA'
media = 'M9glucose'
binding_additions = 'na'
replicate = 3
spec_machine = 'orbitrap_elite'
labeling = 'na'
df_details = pd.DataFrame([[promoter, date, WT_sequence, control_sequence,
strain, media, binding_additions, replicate, spec_machine, labeling]])

#==============================================================================
# Load in the data files
#==============================================================================
# Load the output of protein groups containing SILAC ratios into a DataFrame
df = pd.read_csv('../../data/mass_spec/20170506_titration/proteinGroups.txt',
                 delimiter='	')
# remove non-e.coli items detected
df = df[~df['Protein IDs'].str.contains("REV")]
df = df[~df['Protein IDs'].str.contains("CON")]

# from MaxQuant output, I want to collect the following for each replicate:
df = df[['Protein IDs','Gene names','Peptides',
    'Intensity titration_trial3_01', 'Intensity L titration_trial3_01',
    'Intensity H titration_trial3_01', 'Intensity titration_trial3_1',
    'Intensity L titration_trial3_1', 'Intensity H titration_trial3_1',
    'Intensity titration_trial3_10', 'Intensity L titration_trial3_10',
    'Intensity H titration_trial3_10', 'Intensity titration_trial3_100',
    'Intensity L titration_trial3_100', 'Intensity H titration_trial3_100',
    'Intensity titration_trial3_1000', 'Intensity L titration_trial3_1000',
    'Intensity H titration_trial3_1000',
    'Ratio H/L type', 'Ratio H/L titration_trial3_01',
    'Ratio H/L titration_trial3_1', 'Ratio H/L titration_trial3_10',
    'Ratio H/L titration_trial3_100','Ratio H/L titration_trial3_1000',
    'Ratio H/L normalized titration_trial3_01',
    'Ratio H/L normalized titration_trial3_1',
    'Ratio H/L normalized titration_trial3_10',
    'Ratio H/L normalized titration_trial3_100',
    'Ratio H/L normalized titration_trial3_1000']]

df = df.reset_index(drop=True)

# tile experimental details for number of protein entries in df
df_details = pd.DataFrame(pd.np.tile(df_details, (len(df), 1)))

#==============================================================================
# analysis of SILAC ratio and data correction
#==============================================================================


#==============================================================================
# Generate the summary file
#==============================================================================

# combine mass spec data with experimental details into tidy format
df = pd.concat([df, df_details], axis=1, ignore_index=True)

# Rename the columns of the data_frame
df.columns = ['ProteinIDs','gene','peptides',
            'intensity_total_01', 'intensity_L_01','intensity_H_01',
            'intensity_total_1', 'intensity_L_1','intensity_H_1',
            'intensity_total_10', 'intensity_L_10','intensity_H_10',
            'intensity_total_100', 'intensity_L_100','intensity_H_100',
            'intensity_total_1000', 'intensity_L_1000','intensity_H_1000',
            'MaxQuant_ratio_type', 'maxquant_ratio_01','maxquant_ratio_1',
            'maxquant_ratio_10', 'maxquant_ratio_100', 'maxquant_ratio_1000',
            'maxquant_ratio_normalized_01', 'maxquant_ratio_normalized_1',
            'maxquant_ratio_normalized_10', 'maxquant_ratio_normalized_100',
            'maxquant_ratio_normalized_1000',
            'promoter', 'date', 'wt_sequence', 'control_sequence', 'strain',
            'media', 'binding_additions', 'replicate', 'spec_machine',
            'labeling']


# save to file
df.to_csv(+ str(date) + '_' + promoter + '_' + str(replicate) +
        '_mass_spec_summary.csv', index=False)
