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

#===============================================================================
# Read the data
#===============================================================================

datadir = '../mass_spec/*/'
files = glob.glob(datadir+'*.csv')

df = pd.DataFrame()

for f in enumerate(files):
    if "titration" not in f[1]:
        continue
    df_temp = pd.DataFrame()
    df_temp = pd.read_csv(f[1])

    # let's median shift the maxquant_ratio for sample that is 1:1

    df_temp['maxquant_logratio_1'] = np.log(df_temp['maxquant_ratio_1'])
    median_shift = df_temp['maxquant_logratio_1'].replace([-np.inf,np.inf], np.nan).dropna().median()
    print(np.exp(median_shift))
    df_temp['maxquant_ratio_medianshift_01'] = df_temp['maxquant_ratio_01']/np.exp(median_shift)
    df_temp['maxquant_ratio_medianshift_1'] = df_temp['maxquant_ratio_1']/np.exp(median_shift)
    df_temp['maxquant_ratio_medianshift_10'] = df_temp['maxquant_ratio_10']/np.exp(median_shift)
    df_temp['maxquant_ratio_medianshift_100'] = df_temp['maxquant_ratio_100']/np.exp(median_shift)
    df_temp['maxquant_ratio_medianshift_1000'] = df_temp['maxquant_ratio_1000']/np.exp(median_shift)

    # append data to df
    df = df.append(df_temp)

#===============================================================================
# Calculate the mean and standard deviation of heavy/light ratio for all
# proteins that were detected.
#===============================================================================


df = df.join(df.groupby(['promoter', 'gene', 'strain'])['maxquant_ratio_medianshift_01'].mean(),\
             on=['promoter', 'gene', 'strain'],rsuffix='_avg')
# df = df.join(df.groupby(['promoter', 'gene', 'strain'])['maxquant_ratio_medianshift_01'].std(),\
#              on=['promoter', 'gene', 'strain'],rsuffix='_std')
df = df.join(df.groupby(['promoter', 'gene', 'strain'])['maxquant_ratio_medianshift_1'].mean(),\
             on=['promoter', 'gene', 'strain'],rsuffix='_avg')
# df = df.join(df.groupby(['promoter', 'gene', 'strain'])['maxquant_ratio_medianshift_1'].std(),\
#              on=['promoter', 'gene', 'strain'],rsuffix='_std')
df = df.join(df.groupby(['promoter', 'gene', 'strain'])['maxquant_ratio_medianshift_10'].mean(),\
             on=['promoter', 'gene', 'strain'],rsuffix='_avg')
# df = df.join(df.groupby(['promoter', 'gene', 'strain'])['maxquant_ratio_medianshift_10'].std(),\
#              on=['promoter', 'gene', 'strain'],rsuffix='_std')
df = df.join(df.groupby(['promoter', 'gene', 'strain'])['maxquant_ratio_medianshift_100'].mean(),\
             on=['promoter', 'gene', 'strain'],rsuffix='_avg')
# df = df.join(df.groupby(['promoter', 'gene', 'strain'])['maxquant_ratio_medianshift_100'].std(),\
#              on=['promoter', 'gene', 'strain'],rsuffix='_std')
df = df.join(df.groupby(['promoter', 'gene', 'strain'])['maxquant_ratio_medianshift_1000'].mean(),\
             on=['promoter', 'gene', 'strain'],rsuffix='_avg')
# df = df.join(df.groupby(['promoter', 'gene', 'strain'])['maxquant_ratio_medianshift_1000'].std(),\
#              on=['promoter', 'gene', 'strain'],rsuffix='_std')

#===============================================================================
# lets now take the means above and make a new DataFrame with the expected and
# measured values
#===============================================================================
df_tes = pd.DataFrame()
for i, let in enumerate(df['maxquant_ratio_medianshift_01_avg']):
    try:
        df_tes = df_tes.append([[0.1, let]],ignore_index=True)
    except:
        pass
for i, let in enumerate(df['maxquant_ratio_medianshift_1_avg']):
    try:
        df_tes = df_tes.append([[1, let]], ignore_index=True)
    except:
        pass
for i, let in enumerate(df['maxquant_ratio_medianshift_10_avg']):
    try:
        df_tes = df_tes.append([[10, let]], ignore_index=True)
    except:
        pass
for i, let in enumerate(df['maxquant_ratio_medianshift_100_avg']):
    try:
        df_tes = df_tes.append([[100, let]], ignore_index=True)
    except:
        pass
for i, let in enumerate(df['maxquant_ratio_medianshift_1000_avg']):
    try:
        df_tes = df_tes.append([[1000, let]], ignore_index=True)
    except:
        pass

df_tes = df_tes.rename(columns={0: 'ratio_expect', 1: 'ratio_measured'})

#===============================================================================
# Lets calculate an effecting labeling efficiency to plot for comparison.
# Here we will assume that the heavy label isn't perfectly heavy. This may arise
# due to e.g. the heavy lysine not being 100%.
#===============================================================================

# Lets also plot with a 99.1 percent incorporation (effective)
H = np.logspace(-2,4,50)
#H = np.array([0.1, 1.0, 10.0, 100.0, 1000.0])
H_corr = H * 0.991
L = np.ones(len(H))
#L = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
ratio = H/L
ratio_corr = H_corr/(L+(H-H_corr))

#===============================================================================
# lets plot the ratios
#===============================================================================

fig = plt.figure(figsize=(6.5,4))
ax = fig.add_subplot(111)

sns.stripplot(x='ratio_expect', y='ratio_measured', data=df_tes, jitter=0.1, size=2,
              alpha=0.6, color = "#3498db")

ax.plot(np.linspace(-1,5,7),np.logspace(-2,4,7, base=10), linestyle ='--', alpha = 0.8,
         color = 'grey')
ax.plot(np.linspace(-1,5,50), ratio_corr, linestyle ='-.', alpha = 0.5, color = 'red')

ax.set_ylabel('measured\n(heavy / light)')
ax.set_xlabel('expected\n(heavy / light)')
# ax.ylabel('')
# ax.xlabel('')

ax.set_ylim(ymax =1E4, ymin = 1E-2)
ax.set_yscale('log')

plt.tight_layout()

fig.savefig(output + 'figS5_mass_spec_titration.pdf', format='pdf')
