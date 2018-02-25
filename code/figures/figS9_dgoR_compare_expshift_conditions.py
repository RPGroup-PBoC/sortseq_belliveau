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
utils.set_plotting_style1()

#===============================================================================
# Set output directory
#===============================================================================
output = 'output_figs/'
#------------------------------------------------------------------------------#
# Load in the summary csv file
#------------------------------------------------------------------------------#
fname_1 = 'input_data/fig7_dgoR_expshfit_glucose.csv'
df1 = pd.read_csv(fname_1)

fname_2 = 'input_data/fig7_dgoR_expshfit_galactonate.csv'
df2 = pd.read_csv(fname_2)

fname_3 = 'input_data/fig7_dgoR_deltadgoR_expshift_glucose.csv'
df3 = pd.read_csv(fname_3)

seqLength = 147

#------------------------------------------------------------------------------#
# define linear regression function
#------------------------------------------------------------------------------#

import scipy.optimize

# Define linear function
def linear_fun(x, slope, intercept):
    return slope * x + intercept

#------------------------------------------------------------------------------#
# make plots comparing information footprint and delta bin shift
# across different conditiosn
#------------------------------------------------------------------------------#

colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]

fig, ((ax1, ax2)) = plt.subplots(nrows=1, ncols=2,figsize=(12,4))

#------------------------------------------------------------------------------#
# galactonate vs glucose - MG1655 - bin shift values
#------------------------------------------------------------------------------#
X = df2[df2.position >=-121]['expshift'].values
Y = df1[df1.position >=-121]['expshift'].values

# Compute the curve fit (Guess is unit slope and zero intercept)
popt, _ = scipy.optimize.curve_fit(linear_fun, X, Y,
                                   p0=[1,0])
# Parse the results
slope, intercept = popt

residuals = Y- linear_fun(X, slope, intercept)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((Y-np.mean(Y))**2)
r_squared_1 = 1 - (ss_res / ss_tot)


# Plot data points
ax1.plot(X, Y, '.', alpha=1, color=colours[2])

# Plot best fit line
x = np.linspace(-0.6, 0.6, 200)
y = linear_fun(x, slope, intercept)
ax1.plot(x, y, '-', linewidth=2, zorder=1, color='k', alpha=0.6, linestyle='--')

# # Hide the right and top spines
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
# # Only show ticks on the left and bottom spines
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.set_xlabel('fluorescence bin shift')
ax1.set_xlim(-0.6,0.6)
ax1.set_ylim(-0.6,0.6)
ax1.grid(False)

ax1.set_yticks([-0.6,0,0.6], minor=False)
ax1.set_xticks([-0.6,0,0.6], minor=False)
ax1.yaxis.grid(True, which='major')
ax1.xaxis.grid(True, which='major')
ax1.text(0.3,0.4, r'$R^2$' + ' = %0.2f' %r_squared_1,
        horizontalalignment='center',
        verticalalignment='center', alpha = 0.8, fontsize='20')

#------------------------------------------------------------------------------#
# galactonate vs deltadgoR glucose - MG1655 - bin shift values
#------------------------------------------------------------------------------#
X = df2[df2.position >=-121]['expshift'].values
Y = df3[df3.position >=-121]['expshift'].values

# Compute the curve fit (Guess is unit slope and zero intercept)
popt, _ = scipy.optimize.curve_fit(linear_fun, X, Y,
                                   p0=[1,0])
# Parse the results
slope, intercept = popt

residuals = Y- linear_fun(X, slope, intercept)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((Y-np.mean(Y))**2)
r_squared_2 = 1 - (ss_res / ss_tot)

# Plot data points
ax2.plot(X, Y, '.', alpha=1, color=colours[2])

# Plot best fit line
x = np.linspace(-0.6, 0.6, 200)
y = linear_fun(x, slope, intercept)
ax2.plot(x, y, '-', linewidth=2, zorder=1, color='k', alpha=0.6, linestyle='--')

# ax2.scatter(df2[df2.position >=40]['expshift'].values, df3[df3.position >=40]['expshift'].values, alpha=1, color=colours[2])
# # Hide the right and top spines
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
# # Only show ticks on the left and bottom spines
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')
ax2.set_xlabel('fluorescence bin shift')
ax2.set_xlim(-0.6,0.6)
ax2.set_ylim(-0.6,0.6)

ax2.set_yticks([-0.6,0,0.6], minor=False)
ax2.set_xticks([-0.6,0,0.6], minor=False)
ax2.yaxis.grid(True, which='major')
ax2.xaxis.grid(True, which='major')
ax2.text(0.3,0.4, r'$R^2$' + ' = %0.2f' %r_squared_2,
        horizontalalignment='center',
        verticalalignment='center', alpha = 0.8, fontsize='20')
plt.tight_layout()

figname_out = 'figS9_dgoR_compare_conditions_trimmed.pdf'
fig.savefig(output + figname_out, format='pdf')

print(r_squared_1)
print(r_squared_2)
