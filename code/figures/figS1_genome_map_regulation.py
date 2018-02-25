import numpy as np
import pandas as pd

# plotting stuff
import matplotlib.pyplot as plt
# import pylab as pl
from IPython.core.pylabtools import figsize


#===============================================================================
# Set output directory based on the graphicspath.tex file to print in dropbox
#===============================================================================
output = 'output_figs/'

#==============================================================================#
# Load in the operons and annotated TF binding sites on RegulonDB
#==============================================================================#

# load in data of known TF binding sites
df_bindingsites = pd.read_csv('input_data/network_tf_operon.txt',  delimiter='	', skiprows=16)
# split operon column into two columns; one with operon, the other with associated genes
df_bindingsites[['operon','operon_genes']] = df_bindingsites['operon'].str[:-1].str.split('[', expand=True).astype(str)

# load in list of operons
df_operons = pd.read_csv('input_data/OperonSet.txt',  delimiter='	', skiprows=36)

#==============================================================================#
# Identify operons that have any binding site associated with them.
#==============================================================================#

# groupby operons
operons = df_operons.groupby('operon')

# Go through list of operons and check whether a binding site
# is listed for it in the df_bindingsites dataframe.
# If so, append to operon_reg array and otherwise append
# to operon_noreg array.
operon_reg = []
operon_noreg = []
for op in operons:
    # ignore phage insertions or genes labeled as phantom genes
    if ('ins' in op[0]) or ('C0'  in op[0]) or ('IS128' in op[0]):
        continue
    # check whether binding site is listed for operon op
    d = df_bindingsites[df_bindingsites.operon == op[0]].operon
    if len(d) >= 1:
        operon_reg.append(op[0])
    if len(d) == 0:
        operon_noreg.append(op[0])

#==============================================================================#
# Now we need to make an array that notes the start position
# for each operon - for both operon_reg and operon_noreg arrays.
#==============================================================================#

# generate start position array for operon_reg
operon_reg_pos = []
for op in operon_reg:
    if df_operons[df_operons.operon==op].strand.values == 'forward':
        operon_reg_pos.append(df_operons[df_operons.operon==str(op)].pos_left.values[0])
    else:
        operon_reg_pos.append(df_operons[df_operons.operon==str(op)].pos_right.values[0])

# generate start position array for operon_reg
operon_noreg_pos = []
for op in operon_noreg:
    if df_operons[df_operons.operon==op].strand.values == 'forward':
        operon_noreg_pos.append(df_operons[df_operons.operon==str(op)].pos_left.values[0])
    else:
        operon_noreg_pos.append(df_operons[df_operons.operon==str(op)].pos_right.values[0])

#==============================================================================#
# Calculate fraction of annotated and unannotated operons based on
# data from RegulonDB
#==============================================================================#

# fraction of un-regulated operons:
fraction_reg = float(len(operon_reg))/float(len(operon_reg) + len(operon_noreg))
print('fraction of regulated operons: %.2f' % fraction_reg)

# fraction of regulated operons:
fraction_noreg = float(len(operon_noreg))/float(len(operon_reg) + len(operon_noreg))
print('fraction of unannoated operons: %.2f' % fraction_noreg)#',fraction_noreg)

#==============================================================================#
# convert start positions to degrees for circular map
# E. coli GenBank; U00096.3   4641652 bp
#==============================================================================#

# regulated operon positions
reg_pos_degree = 360.0*(np.sort(operon_reg_pos)/4641652.0)
# unannotated operon positions
noreg_pos_degree = 360.0*(np.sort(operon_noreg_pos)/4641652.0)



# make circular plot of result

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111, projection='polar')

# place line at positions using the axis ticks in a
# polar plot.

# this just sets the length of each tick line
tick1 = [ax.get_rmax()*0.86,ax.get_rmax()*1.14] # marks oriC
tick2 = [ax.get_rmax()*0.94,ax.get_rmax()*1.06] # marks regulated and unannotated operons
tick = [ax.get_rmax()*0.9,ax.get_rmax()*1.1] # marks operons considered in this work.

# note that position of promoters are based on
# first gene; gene and oriC degree obtained on EcoCyc
oriC_pos = 56  # is at 304 degree, shift to make 0
for t  in np.deg2rad([304+oriC_pos]):
    ax.plot([t,t], tick1, lw=4, color="k")


# regulated and non-regulated genes
for t  in np.deg2rad(-(noreg_pos_degree+oriC_pos)):
    ax.plot([t,t], tick2, lw=0.6, color="#EE352B", alpha = 1)
for t  in np.deg2rad(-(reg_pos_degree+oriC_pos)):
    ax.plot([t,t], tick2, lw=0.6, color="#3954A4", alpha = 1)


# tick = [ax.get_rmax()*0.85,ax.get_rmax()*1.15]
# lacZ
for t  in np.deg2rad([-(28+oriC_pos)]):
    ax.plot([t,t], tick, lw=10, color="#C24E9D")

# relB
for t  in np.deg2rad([-(128+oriC_pos)]):
    ax.plot([t,t], tick, lw=10, color="#E55E68")

# marR
for t  in np.deg2rad([-(126+oriC_pos)]):
    ax.plot([t,t], tick, lw=10, color="#F47530")

# purT (purT and yebG at 150)
for t  in np.deg2rad([-(151+oriC_pos)]):
    ax.plot([t,t], tick, lw=10, color="#19733A")

# yebG
for t  in np.deg2rad([-(149+oriC_pos)]):
    ax.plot([t,t], tick, lw=10, color="#43669A")

# xylE
for t  in np.deg2rad([-(329+oriC_pos)]):
    ax.plot([t,t], tick, lw=10, color="#3C8DA0")

# dgoR
for t  in np.deg2rad([-(300+oriC_pos)]):
    ax.plot([t,t], tick, lw=10, color="#E8B770")


ax.grid(False)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.spines['polar'].set_visible(False)
# ax.set_theta_offset(np.deg2rad(oriC_pos))

figname_out = output + 'figS1_genome_regulation_summary.pdf'
plt.savefig(figname_out, format='pdf')
