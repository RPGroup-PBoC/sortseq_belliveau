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

date = 20170229

#==============================================================================
# Load in the data files
#==============================================================================

# List the target directory.
datadir = '/Users/Nathan/Desktop/20170509_FACS_organization_temp/'
files = np.array(os.listdir(datadir))
# Find file containing strain and measurement details
csv_bool_details = np.array([str(date) in f and ('csv' in f or 'txt' in f)\
                    and 'details' in f for f in files])
files_details = files[np.array(csv_bool_details)]
# Find file containing measurement statistics
csv_bool_stats = np.array([str(date) in f and ('csv' in f or 'txt' in f)\
                    and 'stats' in f for f in files])
files_stats = files[np.array(csv_bool_stats)]
# Find file containing measurement histgram counts
csv_bool_hist = np.array([str(date) in f and ('csv' in f or 'txt' in f)\
                    and 'hist' in f for f in files])
files_hist = files[np.array(csv_bool_hist)]
# Find comments file
csv_bool_comments = np.array([str(date) in f and ('csv' in f or 'txt' in f)\
                    and 'comments' in f for f in files])
files_comments = files[np.array(csv_bool_comments)]

print(files_details[0])

#==============================================================================
# Generate the summary stats file
#==============================================================================

# Initialize the DataFrame to hold the tidy formated measurement details.
df = pd.read_csv(files_details[0])
# Add in the histogram stats from the flow cytometry measurements
dataframe_stats = pd.read_csv(files_stats[0],sep='\t',skiprows=1)
df[['Region','Count', '%Hist', '%All', 'Bounds', 'ModeCount',
       'Mode', 'Mean', 'Median', 'StdDev.', 'CV', 'CV(hm)', 'Skew']] = \
       dataframe_stats[['Region','Count', '%Hist', '%All', 'Bounds',
         'ModeCount', 'Mode', 'Mean', 'Median', 'StdDev.', 'CV', 'CV(hm)',
         'Skew']]

# write
df.to_csv('output/tmp.csv', index=False)

# Add the comments to the header of the data file
filenames = [files_comments[0], 'output/tmp.csv']

with open('output/' + str(date) + '_flow_stats.csv', 'w') as output:
    for fname in filenames:
        with open(fname) as infile:
            output.write(infile.read())

# Remove temporary file
os.remove(filenames[1])

#==============================================================================
# Generate the histogram counts file
#==============================================================================

hist = np.loadtxt(files_hist[0])
x = np.logspace(0.1,5,num=256)

# Initialize the DataFrame to save the mean expression levels
df2 = pd.DataFrame()
for index, row in df.iterrows():
    for i, counts in enumerate(hist[index,:]):
        # df2 = df2.append(row[['date', 'promoter', 'strain', 'region', 'bin',
        #  'sequence', 'voltage', 'gain', 'media', 'condition']],ignore_index=True)
        df2 = df2.append([[row['date'], row['promoter'], row['strain'], row['region'], row['bin'],
         row['sequence'], row['voltage'], row['gain'], row['media'], row['condition'],x[i],counts, counts/hist[index,:].sum()]],ignore_index=True)
# Rename the columns of the data_frame
df2.columns = ['date', 'promoter', 'strain', 'region', 'bin',
              'sequence', 'voltage', 'gain', 'media', 'condition',
              'fluorescence', 'count', 'fraction']

print(df2.columns.values)
df2.to_csv('output/tmp_2.csv', index=False)

# Add the comments to the header of the data file
filenames = [files_comments[0], 'output/tmp_2.csv']

with open('output/' + str(date) + '_flow_hist.csv', 'w') as output:
    for fname in filenames:
        with open(fname) as infile:
            output.write(infile.read())
# Remove temporary file
os.remove(filenames[1])
