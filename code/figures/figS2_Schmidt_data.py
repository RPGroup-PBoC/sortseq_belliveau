import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.special
import scipy.optimize
import scipy.io
import glob

# Import the project utils
import sys
sys.path.insert(0, '../')
import NB_sortseq_utils as utils

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Seaborn, useful for graphics
import seaborn as sns

sns.set_palette("deep", color_codes=True)
utils.set_plotting_style1()


#===============================================================================
# Set output directory based on the graphicspath.tex file to print in dropbox
#===============================================================================
output = 'output_figs/'

#==============================================================================#
# plot bar plot with a subset of the protein measurements for galE and  dgoD
#==============================================================================#
# load in the Schmidt data (partial data from SI)
data = pd.read_csv('input_data/Nature_Bio_2015_Schmidt_supS7_2.csv')
# a subset of the protein measurements from Schmidt et al. 2016
df = pd.read_csv('input_data/schmidt_data_trimmed.csv')

fig1 = plt.figure(figsize=(6,4))
ax1 = plt.subplot(111)

# plt.rcParams['xtick.major.pad']='20'

objects = ['Glucose', 'Xylose', 'Acetate', 'Galactose', 'Glycerol']
y_pos = np.arange(len(objects))
protein_copy = np.array([data[data.Gene == 'galE'].Glucose.values[0],
                data[data.Gene == 'galE'].Xylose.values[0],
                data[data.Gene == 'galE'].Acetate.values[0],
                data[data.Gene == 'galE'].Galactose.values[0],
                data[data.Gene == 'galE'].Glycerol.values[0]])

ax1.bar(y_pos, protein_copy, align='center')
ax1.set_xticks(y_pos, objects)#, fontsize=12)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_yticks(np.arange(0, 6000, 1000))
ax1.set_ylabel('GalE copy\nnumber / cell')

plt.tight_layout()

figname_out = output + 'figS2_Schmidt_galE.pdf'
fig1.savefig(figname_out, format='pdf')


fig2 = plt.figure(figsize=(6,4))
ax2 = plt.subplot(111)

objects = ['Glucose', 'Xylose', 'Acetate', 'Galactose', 'Glycerol']
y_pos = np.arange(len(objects))
protein_copy = [data[data.Gene == 'dgoD'].Glucose.values[0],
                data[data.Gene == 'dgoD'].Xylose.values[0],
                data[data.Gene == 'dgoD'].Acetate.values[0],
                data[data.Gene == 'dgoD'].Galactose.values[0],
                data[data.Gene == 'dgoD'].Glycerol.values[0]]

ax2.bar(y_pos, protein_copy, align='center')
ax2.set_xticks(y_pos, objects)#, fontsize=12)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.set_yticks(np.arange(0, 700, 100))
ax2.set_ylabel('DgoD copy\nnumber / cell')

# plt.title('Programming language usage')
plt.tight_layout()
fig2.savefig(output + 'figS2_Schmidt_dgoD.pdf', format='pdf')

plt.close()
#==============================================================================#
# plot coefficient of variation across all the Schmidt et al. data
#==============================================================================#


# now I need to load the RegulonDB file that lists which genes are regulated (and by what)
regDB_TF_genes = pd.read_csv('input_data/network_tf_gene_.txt',sep='\t')

# determine which genes in the Schmidt data are regulated.
data['TF Annotated'] = np.where(data['Gene'].isin(regDB_TF_genes['Gene']), 'Yes', 'No')

#determine mean values
data['Copy number mean'] = data.values[:,7:29].mean(dtype=np.float64,axis=1)
data['Copy number std dev'] = data.values[:,7:29].std(dtype=np.float64,axis=1)
data['Copy number max'] = data.values[:,7:29].max(axis=1).astype(float)
data['Copy number min'] = data.values[:,7:29].min(axis=1).astype(float)
data['range'] = data['Copy number max'] - data['Copy number min']

# calculate coefficient of variation (std. dev. / mean)
data['cov'] = data['Copy number std dev'] / data['Copy number mean']

conditions=[ 'Glucose', 'LB', 'Glycerol + AA', 'Acetate' , 'Fumarate', 'Glucosamine', 'Glycerol',
 'Pyruvate', 'Chemostat mu=0.5', 'Chemostat mu=0.35', 'Chemostat mu=0.20', 'Chemostat mu=0.12',
'Stationary phase 1 day', 'Stationary phase 3 days', 'Osmotic-stress glucose', '42C glucose',
'pH6 glucose', 'Xylose', 'Mannose', 'Galactose', 'Succinate', 'Fructose']

#randomize data since higher copynumber/ higher confidence measurements are at top of spreadsheet.
data2 = data.reindex(np.random.permutation(data[data['cov']>0].index))


f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 4))
f.subplots_adjust(wspace=0.1)


yvals = data2[data2['TF Annotated']=='Yes']
vals = np.random.uniform(0.98,1.02,size=len(yvals))
ax1.errorbar(yvals['cov'].values, vals,
            linestyle='none',fmt='o', alpha = 0.6, markersize=8, color='#3953a4')
ax1.set_ylim(0.95,1.05)
ax1.set_yticklabels([''])
ax1.yaxis.grid(False)
#plt.xscale('log')


#xylE
count = 0
for i in yvals.Gene:
    if i == 'galE':
        break
    count += 1
yvals_dgoR = data2[(data2['Gene']=='galE')]
print(yvals_dgoR['cov'].values)
val_dgoR = count
print(val_dgoR)
ax1.errorbar(yvals_dgoR['cov'].values, vals[val_dgoR],
            linestyle='none',fmt='o', alpha = 1, markersize=8, color = 'k')
ax1.text(yvals_dgoR['cov'].values, vals[val_dgoR]+0.005, 'galE')#, fontsize=14)

ax1.set_ylim(0.95,1.05)
ax1.set_yticklabels([''])
ax1.yaxis.grid(False)
#plt.xscale('log')


yvals = data2[data2['TF Annotated']=='No']
vals = np.random.uniform(0.98,1.02,size=len(yvals))
ax2.errorbar(yvals['cov'].values, vals,
            linestyle='none',fmt='o', alpha = 0.6, markersize=8, color = '#ee332a')
#dgoD
count = 0
for i in yvals.Gene:
    if i == 'dgoD':
        break
    count += 1
yvals_dgoR = data2[(data2['Gene']=='dgoD')]
print(yvals_dgoR['cov'].values)
val_dgoR = count
print(val_dgoR)
ax2.errorbar(yvals_dgoR['cov'].values, vals[val_dgoR],
            linestyle='none',fmt='o', alpha = 1, markersize=8, color = 'k')
ax2.text(yvals_dgoR['cov'].values, vals[val_dgoR]+0.005, 'dgoD', fontsize=14)
#dgoK
count = 0
for i in yvals.Gene:
    if i == 'dgoK':
        break
    count += 1
yvals_dgoR = data2[(data2['Gene']=='dgoK')]
print(yvals_dgoR['cov'].values)
val_dgoR = count
print(val_dgoR)
ax2.errorbar(yvals_dgoR['cov'].values, vals[val_dgoR],
            linestyle='none',fmt='o', alpha = 1, markersize=8, color = 'k')
ax2.text(yvals_dgoR['cov'].values, vals[val_dgoR]+0.005, 'dgoK', fontsize=14)
#dgoA
count = 0
for i in yvals.Gene:
    if i == 'dgoA':
        break
    count += 1
yvals_dgoR = data2[(data2['Gene']=='dgoA')]
print(yvals_dgoR['cov'].values)
val_dgoR = count
print(val_dgoR)
ax2.errorbar(yvals_dgoR['cov'].values, vals[val_dgoR],
            linestyle='none',fmt='o', alpha = 1, markersize=8, color = 'k')
ax2.text(yvals_dgoR['cov'].values, vals[val_dgoR]+0.005, 'dgoA', fontsize=14)
#dgoR
count = 0
for i in yvals.Gene:
    if i == 'dgoR':
        break
    count += 1
yvals_dgoR = data2[(data2['Gene']=='dgoR')]
print(yvals_dgoR['cov'].values)
val_dgoR = count
print(val_dgoR)
ax2.errorbar(yvals_dgoR['cov'].values, vals[val_dgoR],
            linestyle='none',fmt='o', alpha = 1, markersize=8, color = 'k')
ax2.text(yvals_dgoR['cov'].values, vals[val_dgoR]+0.005, 'dgoR', fontsize=14)

#xylE
count = 0
for i in yvals.Gene:
    if i == 'xylE':
        break
    count += 1
yvals_dgoR = data2[(data2['Gene']=='xylE')]
print(yvals_dgoR['cov'].values)
val_dgoR = count
print(val_dgoR)
ax2.errorbar(yvals_dgoR['cov'].values, vals[val_dgoR],
            linestyle='none',fmt='o', alpha = 1, markersize=8, color = 'k')
ax2.text(yvals_dgoR['cov'].values, vals[val_dgoR]+0.005, 'xylE', fontsize=14)


ax2.set_ylim(0.95,1.05)
ax2.set_yticklabels([''])
ax2.yaxis.grid(False)
ax2.tick_params(axis='both', which='major', labelsize=20)

ax2.set_xlim(0,4.5)
#plt.xscale('log')
plt.tight_layout()
plt.savefig(output + 'figS2_Schmidt_stripplot_cov.pdf', format='pdf')
