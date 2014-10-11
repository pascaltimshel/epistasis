#!/usr/bin/env python

import os
import sys
import numpy as np
import pandas as pd

# this module contains Matlab-like commands
import matplotlib.pyplot as plt

import argparse

###################################### USAGE ######################################

#python gen_interacting_SNP_sets.py --width 500000

###################################### FILE snippets ######################################
### file_bim
# 1       rs11497407      0       558390  0       2
# 1       rs12565286      0       711153  0       2
# 1       rs11804171      0       713682  0       2
# 1       rs2977670       0       713754  0       1


### file_interactions
### IMPORTANT information: 
#1) the first column listing the position will be used (e.g. 67280788 and not 67280789)
#2) the left part of the file (column 1-4) will be denoted interaction_XX_A. 
#3) the right part of the file (column 5-8) wil be denoted interaction_XX_B.
#4) THE last column contains "updated" interaction IDs. There are 92932 of these updated (there where 96137 of the original)
# chr12   67280788        67280789        interaction1    chr12   66561751        66561752        interaction1    interaction1
# chr3    150566213       150566214       interaction2    chr3    150061015       150061016       interaction2    interaction2
# chr10   28973034        28973035        interaction3    chr10   28917902        28917903        interaction3    interaction3
# chr10   28978176        28978177        interaction4    chr10   29702662        29702663        interaction4    interaction4
# chr17   52725358        52725359        interaction5    chr17   53563372        53563373        interaction5    interaction5


### OSX
file_bim = "/Users/pascaltimshel/Dropbox/EGCUT_DATA/geno/Prote_370k_251011.no_mixup.chr_infered.bim"
file_interactions = "/Users/pascaltimshel/p_HiC/Lan_et_al_chromosomal_interactions/lift_findItersection.intersection.paste.updatedIDs"

### Broad
#XX

###################################### PARAMETERs ###################################### 
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--width", help="INTERGER. Width in bp. Width is the distance to boundary from the focal/center point. x - width <= snp_pos <= x + width", default=500*1000)
args = arg_parser.parse_args()


interaction_width = int(args.width)
#interaction_width = 500*1000 # in bp. try [20kb, 500kb, 2Mb and 10Mb] or [10kb, 250kb, 1Mb and 5Mb], where the numbers are TWICE the distance to the region boundary from the focal point.

path_base_out = "/Users/pascaltimshel/p_HiC/Lan_et_al_interaction_SNP_sets/{}".format(interaction_width) # OSX
path_snp_sets = path_base_out+"/snp_sets"
path_figs = path_base_out+"/figs"
path_stats = path_base_out+"/stats"

paths = [path_base_out, path_snp_sets, path_figs, path_stats]
for path in paths:
	if not os.path.exists(path):
		os.makedirs(path)

## figs
## stats
## snpsets

###################################### READ bim file ######################################
header_bim = ["chr", "snp", "cm", "pos", "dummy1", "dummy2"]
df_bim = pd.read_csv(file_bim, sep="\t", header=None, names=header_bim)


###################################### Processing ######################################
df_set_stats = pd.DataFrame() #columns=['setA_size', 'setA_range' 'setB_size', 'setB_range']
with open(file_interactions, 'r') as f:
	for line_no, line in enumerate(f):
		if line_no == 100: break
		
		fields = line.strip().split()
		(chr_A, pos_A, chr_B, pos_B, interaction_ID) = (fields[0], fields[1], fields[4], fields[5], fields[8])
		
		## Turn off X chromosomes
		#if (chr_A.upper() == 'CHRX') or (chr_B.upper() == 'CHRX'):
		#    continue
		## Modifying chromosomes
		if 'X' in chr_A.upper(): chr_A = chr_A.replace('X', '23')
		if 'X' in chr_B.upper(): chr_B = chr_B.replace('X', '23')
		chr_A = int(chr_A.lstrip('chr'))
		chr_B = int(chr_B.lstrip('chr'))
		pos_A = int(pos_A)
		pos_B = int(pos_B)
		if line_no % 50 == 0: print line_no, chr_A, pos_A, chr_B, pos_B, interaction_ID
		
		### A
		df_A_extract = df_bim[df_bim["chr"]==chr_A]
		df_A_extract = df_A_extract[(df_A_extract["pos"] >= pos_A-interaction_width) & (df_bim["pos"] <= pos_A+interaction_width)]
		### B
		df_B_extract = df_bim[df_bim["chr"]==chr_B]
		df_B_extract = df_B_extract[(df_B_extract["pos"] >= pos_B-interaction_width) & (df_bim["pos"] <= pos_B+interaction_width)]       
		
		## set calculations
		tmp_intersection = len(set(df_A_extract['snp']).intersection(set(df_B_extract['snp'])))
		tmp_union = len(set(df_A_extract['snp']).union(set(df_B_extract['snp'])))
		## number of tests
		n_tests_plink = len(df_A_extract)*len(df_B_extract)-tmp_intersection
		n_tests_non_redundant = len(df_A_extract)*len(df_B_extract)-tmp_intersection**2+(tmp_intersection-1)*tmp_intersection/2 # int division ok! # N_A*N_B-N_AB^2+(N_AB-1)*N_AB/2
		## df_stats
		df_stats = pd.DataFrame({'set_interchromosomal':1 if chr_A != chr_B else 0,
								 'set_distance': abs(pos_A - pos_B) if chr_A == chr_B else None, # or np.nan
								'set_intersect':tmp_intersection,
								'set_union':tmp_union,
								'set_self_interaction':1 if (chr_A == chr_B) and (pos_A == pos_B) else 0,
								'set_percentage_shared':tmp_intersection/float(tmp_union)*100 if tmp_union!=0 else None, # float divison
								'set_n_tests_plink':n_tests_plink,
								'set_n_tests_non_redundant': n_tests_non_redundant,
								 'setA_size':len(df_A_extract), 
								 'setA_range':df_A_extract["pos"].max()-df_A_extract["pos"].min(), 
								 'setB_size':len(df_B_extract), 
								 'setB_range':df_B_extract["pos"].max()-df_B_extract["pos"].min()}, 
								index=[interaction_ID])
		df_set_stats = df_set_stats.append(df_stats) # index and columns are appended to the end of the dataframe
		# If no columns are passed, the columns will be the sorted list of dict keys.
		
		file_setA = path_snp_sets + "/{}_A.txt".format(interaction_ID)
		file_setB = path_snp_sets + "/{}_B.txt".format(interaction_ID)
		with open(file_setA, 'w') as fA:
			for snp in df_A_extract["snp"]:
				fA.write(snp+"\n")
		with open(file_setB, 'w') as fB:
			for snp in df_B_extract["snp"]:
				fB.write(snp+"\n")

###################################### Creating summary  ######################################
df_set_summay = pd.Series(
	{
	'setA_size_mean': df_set_stats['setA_size'].mean(),
	'setB_size_mean': df_set_stats['setB_size'].mean(),
	'set_n_interchromosomal_sum': df_set_stats['set_interchromosomal'].sum(),
	'set_n_tests_non_redundant_sum':df_set_stats['set_n_tests_non_redundant'].sum(),
	'set_n_tests_plink_sum':df_set_stats['set_n_tests_plink'].sum(),
	'set_intersect_mean':df_set_stats['set_intersect'].mean(),
	'set_percentage_shared_sum':df_set_stats['set_percentage_shared'].mean()
	})
###################################### Write stats to csv ######################################
df_set_stats.to_csv(path_stats+"/df_set_stats.txt", sep="\t") # Series
df_set_summay.to_csv(path_stats+"/df_set_summay.txt", sep="\t", index=True) # DataFrame

###################################### Plotting ######################################

for column_name, series in df_set_stats.iteritems(): # Iterator over (column, series) pairs
	print column_name
	file_figure = "{}/fig_{}.pdf".format(path_figs, column_name)

	plt.figure()
	_ = plt.hist(series.dropna()) # needed for NaN values in distance
	plt.title(column_name, fontsize=18)
	plt.xlabel(column_name)
	plt.ylabel('counts')
	#plt.savefig(file_figure)
	plt.savefig(file_figure)
	plt.close() # or plt.close(fig)

