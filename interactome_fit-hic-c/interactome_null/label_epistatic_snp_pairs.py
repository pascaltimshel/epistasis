#!/usr/bin/env python

import os
import sys

import matplotlib
matplotlib.use('Agg') #Agg backend and not an X-using backend that required an X11 connection. Call use BEFORE importing pyplot!
# REF: http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pyplot as plt

import argparse

import re
import collections
import pandas as pd
import numpy as np # for plotting an array

###################################### USAGE ######################################

#python gen_SNP2interaction_map.py --path_interaction_table XXX --file_null_table XXX --q_threshold XXX --path_main_out XXX


### Examples 1 - fastEpi_compiled/hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8
#python gen_SNP2interaction_map.py --path_interaction_table /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_snpsets/maf_5_sets/500_snppool_hIMR90_q_1e-08 --file_null_table /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/null_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-08.nperm_1000.txt --q_threshold 1e-08 --path_main_out /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8

### Examples 2 - hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10 (not complete yet)
#python gen_SNP2interaction_map.py --path_interaction_table /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_snpsets/maf_5_sets/50000_snppool_hIMR90_q_1e-09 --file_null_table /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/null_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-09.nperm_1000.txt --q_threshold 1e-09 --path_main_out /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10
	# at null_192_1017 --> Segmentation fault: 11 (likely because it was using 6 GB memory)



###################################### SYNOPSIS ######################################
# ???

### Output files
# ??

###################################### TODO ######################################

###################################### TESTs to make ######################################


###################################### FILE snippets ######################################
################## file_fastepistasis_lm_combined ##################
### Example ### 
# CHR	SNP_A	CHR	SNP_B	BETA	CHISQ	PVALUE	PHENOTYPE
# 11	rs1508531	11	rs4939298	-0.13990	33.32995	7.77755E-09	ILMN_1716816
# 1	rs4650623	6	rs3798983	-0.22814	33.00383	9.19773E-09	ILMN_1716816
# 1	rs4650623	6	rs3798982	-0.22814	33.00383	9.19773E-09	ILMN_1716816
# 4	rs12501969	21	rs762173	-0.04258	32.96917	9.36321E-09	ILMN_1666935
### Expected format
# a tab file with the *COLUMN NAMES* in the format: {experiment_type}_{experiment_no}
# experiment_type:	{hic, null}
# experiment_no:	{1, 2, .., N_perm}
# e.g. hic_1 or null_1, null_2, ..

################## file_interaction_table_snps ##################
### Note that columns "snps_A" and "snps_B" can be empty
# interaction_ID	chr_A	pos_A	chr_B	pos_B	setA_size	setB_size	snps_A	snps_B
# interaction_1	10	100515000	15	100515000	24	11	rs7912656;rs7913008;rs11189802;rs11189803;rs7474758;rs10786469;rs10883191;rs11189804;rs2065877;rs4582899;rs11189811;rs7902691;rs7906326;rs2801410;rs11189812;rs1336503;rs1336505;rs2785218;rs2801392;rs10883197;rs10883198;rs2785219;rs2801393;rs11189814	rs3784812;rs7172796;rs7172808;rs8039168;rs8040533;rs4489968;rs7172038;rs7178084;rs11634550;rs3940418;rs16957846
# interaction_2	3	52755000	10	52755000	6	7	rs11130323;rs2268026;rs2072390;rs13082208;rs13082960;rs2336545	rs7919460;rs10883534;rs10883535;rs12254876;rs7071427;rs7911997;rs7900763
# interaction_3	3	50495000	10	50495000	9	26	rs916288;rs736471;rs12492113;rs1467914;rs6786523;rs12494849;rs1467913;rs9867588;rs763030	rs12777801;rs7094916;rs17221456;rs17221463;rs10737004;rs7070724;rs7095853;rs7900801;rs7071250;rs7918118;rs2242271;rs7076107;rs10752022;rs11253558;rs11253559;rs2306409;rs11253560;rs7099607;rs2620942;rs7917055;rs1871621;rs11253564;rs12415759;rs2387213;rs2306408;rs11253567
# interaction_4	10	103175000	17	103175000	16	4	rs12355803;rs12258171;rs11191002;rs10786635;rs10159775;rs10883648;rs10883649;rs4917940;rs7082055;rs7082129;rs11191009;rs10883650;rs4919545;rs7900797;rs10786636;rs12261987	rs12944420;rs12603358;rs11654121;rs4792930
# interaction_5	10	103255000	15	103255000	10	4	rs9420832;rs11191030;rs4244345;rs9420838;rs9420839;rs9420843;rs9420847;rs9419921;rs9325503;rs12769629	rs963375;rs4777503;rs10438367;rs17823279
# interaction_6	10	10345000	22	10345000	2	14	rs1413288;rs7905520	rs139513;rs139515;rs5751071;rs5751072;rs139516;rs2284079;rs139519;rs139520;rs2235852;rs139525;rs139526;rs139528;rs8139515;rs60171
# interaction_7	3	110875000	10	110875000	7	0	rs1732207;rs12492614;rs1614658;rs1729585;rs1406493;rs1113042;rs1358591	


###################################### PARAMETERs ###################################### 
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--path_interaction_table", required=True, help="""
	*Main* path where the interaction table sits. E.g. XXX/10000_snppool_hIMR90_q_1e-08 and NOT XXXX/10000_snppool_hIMR90_q_1e-08/stats/df_interaction_table_snps.txt
	Expected folder organization in path_interaction_table:
	<path_interaction_table>/
		snp_sets/ [optional]
		figs/ [optional]
		stats/ [required!]
			df_interaction_table_snps.txt
			MORE FILES...
		errors/ [optional]
	""")
arg_parser.add_argument("--file_null_table", required=True, help="A full filename (path+filename) to the null table")
arg_parser.add_argument("--q_threshold", required=True, help="The q-value threshold used. Formatting matters! Example: '1e-08'. This value is only used for validity check and forming output filenames. [PLEASE NOTE THAT THE CHECK IS NOT PERFECT!]")
arg_parser.add_argument("--path_main_out", required=True, help="Main path to write output files. Path will be created if it does not exists.")
args = arg_parser.parse_args()

path_interaction_table = os.path.abspath(args.path_interaction_table)
file_null_table = os.path.abspath(args.file_null_table)
q_threshold = args.q_threshold
path_main_out = os.path.abspath(args.path_main_out) + "/SNP2interaction_map"

file_interaction_table_snps = path_interaction_table+"/stats/df_interaction_table_snps.txt"

###################################### INITIAL Error checks ######################################

################## Check that files exists ##################
if not os.path.exists(file_interaction_table_snps):
	raise Exception("file_interaction_table_snps does not exists: %s" % file_interaction_table_snps)
if not os.path.exists(file_null_table):
	raise Exception("file_null_table does not exists: %s" % file_null_table)

################## Make sure that q_threshold and the input files/paths match ##################
# This script is only sensitive to the correct q_threshold because we need the null and df_interaction_table_snps to match up. The length of the files must be equal
if not re.search(q_threshold, path_interaction_table, flags=re.IGNORECASE): # re.search(pattern, string, flags=0). The re.search function returns a search object on success, None on failure.
	print "q_threshold", q_threshold 
	print "path_interaction_table", path_interaction_table
	raise Exception("Could regex match q_threshold pattern in *path_interaction_table*. Check that the q_threshold is correct")
if not re.search(q_threshold, file_null_table, flags=re.IGNORECASE): # re.search(pattern, string, flags=0). The re.search function returns a search object on success, None on failure.
	raise Exception("Could regex match q_threshold pattern in *file_null_table*. Check that the q_threshold is correct")



###################################### OUTPUT ######################################
if not os.path.exists(path_main_out):
	print "path_main_out did not exists. Will make new dir"
	os.makedirs(path_main_out)

file_SNP2interaction_map = path_main_out + "/SNP2interaction_map.txt" # add information about origin of file
file_bonferroni_correction = path_main_out + "/bonferroni_correction.txt"
#file_interactions_per_snp = path_main_out + "/interactions_per_snp.txt"

file_fig_n_tests_barplot = path_main_out + "/fig_n_tests_barplot.pdf"
file_fig_interactions_per_snp_histogram = path_main_out + "/fig_interactions_per_snp_histogram.pdf"



###################################### READ FILES ######################################
### Interaction table with SNPs
df_interaction_table_snps = pd.read_csv(file_interaction_table_snps, sep="\t")

### Null table
df_null_table = pd.read_csv(file_null_table, sep="\t")

###################################### Initialize container ######################################
### Types of container:
# 1) Defaultdict + set
# 2) Defaultdict + OrderedDict
# 3) Deque

SNP2interaction_dict = collections.defaultdict(set) # using SET()
#SNP2interaction_dict = collections.defaultdict(collections.OrderedDict) # OrderedDict is not a lot slower than a plain dict, but at least doubles the memory.
#SNP2interaction_dict = {}

bonferroni_correction_dict = collections.defaultdict(int)
#SNP_interaction_count_dict = {}

###################################### MAIN LOOP ######################################
for column_name, series in df_null_table.iteritems(): # Iterator over (column, series) pairs
	print column_name
	series_corrected_offset = series-1 # We need to subtract one from the index because Pandas/Python is zero-based and R is one-based
	for interaction_no, elem in enumerate(series_corrected_offset, start=1):
		## Note that interaction_no starts at one! This is because we want the interaction IDs to run from 1...N_interactions.
		
		### This is the index to access rows/interaction in Pandas Dataframes
		interaction_idx = interaction_no - 1

		snps = list() # Default value - empty list. Will be evaluated as False in boolean context.
		
		#if interaction_no==100: break

		try:
			snps_A = df_interaction_table_snps["snps_A"][interaction_idx].split(";")
			snps_B = df_interaction_table_snps["snps_B"][elem].split(";")
			snps = snps_A + snps_B # similar to .extend(). OBS: this list could potentially contain DUPLICATES. However, that should not be a problem
		except AttributeError: # E.g. "AttributeError 'float' object has no attribute 'split'" when splitting "nan" ("nan" is type "float")")
			pass
			# --> Got AttributeError. Value is likely 'nan'. Will do NOTHING because this is ok! 

		### Constructing interaction_identifier
		# column_name example: hic_1, null_1, null_2, ...
		experiment_type = column_name.split("_")[0] # {hic, null}
		experiment_no = column_name.split("_")[1] # {1, 2, .., N_perm}
		interaction_identifier = "{experiment_type}_{experiment_no}_{interaction_no}".format(experiment_type=experiment_type, experiment_no=experiment_no, interaction_no=interaction_no)
		print interaction_identifier
		
		for snp in snps: 
			### Remarks:
			# "snps" may contain duplicates if snps_A and snps_B overlap because we are USING LIST PLUS (+) OPERATOR: snps_A + snps_B
			# duplicates is not a problem for this code. Just think about it.
			# duplicates should not exists when using INTER-CHROMOSOMAL interactions

			SNP2interaction_dict[snp].add(interaction_identifier) # using set()
			#SNP2interaction_dict[snp][interaction_identifier] = 1 # using OrderedDict
			#SNP2interaction_dict[snp] = 1

		################## Bonferroni correction calculation ##################
		# Purpose: we multiply the number of SNPs in setA and setB to get the number of tests.
		# Because I choose to generate the the "null" by shuffling B and KEEPING A FIXED, both "hic" and "null" will have the same setA_size for a given interaction number.
		# setB is a ordered sequence, 1..N_interactions, for "hic". setB is shuffled for the "null".
		# IMPORTANT: *elem* is the INDEX pointing to a interaction in the df_interaction_table_snps.
		# *elem* comes from series_corrected_offset where the numbers are corrected to be zero-based
		tmp_correction = df_interaction_table_snps.ix[interaction_idx,"setA_size"] * df_interaction_table_snps.ix[elem,"setB_size"]
		bonferroni_correction_dict[column_name] += tmp_correction
			
###################################### Writing files ######################################

# print "Writing interactions_per_snp"
# ################## file_interactions_per_snp ##################
# with open(file_interactions_per_snp, 'w') as fh:
# 	for key in sorted(SNP2interaction_dict, key=lambda x: x.split("_")[-1]):
# 		fh.write( "{}\t{}\n".format(key, len(SNP2interaction_dict[key])) )


print "Writing bonferroni_correction"
################## file_bonferroni_correction ##################
with open(file_bonferroni_correction, 'w') as fh_bon:
	for key in sorted(bonferroni_correction_dict, key=lambda x: x.split("_")[-1]):
		fh_bon.write( "{}\t{}\n".format(key, bonferroni_correction_dict[key]) )
		#print key, bonferroni_correction_dict[key]


print "Writing file_SNP2interaction_map"
################## file_SNP2interaction_map ##################
with open(file_SNP2interaction_map, 'w') as fh_map:
	for snp in sorted(SNP2interaction_dict): # .items() does not work because we are using sorted() 
		set_string = ";".join(sorted(SNP2interaction_dict[snp])) # USING set() sorted set
		#set_string = ";".join(SNP2interaction_dict[snp]) # # USING OrderedDict()
		fh_map.write( "{}\t{}\t{}\n".format(snp, len(SNP2interaction_dict[snp]), set_string) )


###################################### Plotting ######################################

################## Histogram ##################
print "will plot histogram now..."

x = np.array( [len(value) for value in SNP2interaction_dict.itervalues()] ) # note the use of .itervalues()
plt.figure()
_ = plt.hist(x, bins=500) # 100 looks ok
plt.title('interactions_per_snp', fontsize=18)
plt.xlabel('interactions_per_snp')
plt.ylabel('Counts')
plt.savefig(file_fig_interactions_per_snp_histogram)
plt.close() # or plt.close(fig)


################## Barplot ##################
### Simple bar-chart
# http://cs.smith.edu/dftwiki/index.php/MatPlotLib_Tutorial_1

print "will plot barplot now..."

width = 1
N = len( bonferroni_correction_dict )
x = np.arange(1, N+1)
y = np.array(bonferroni_correction_dict.values()) # can also be list?
labels = bonferroni_correction_dict.keys()

print "number of observations for plotting: %s" % N

try:
	plt.figure()
	_ = plt.bar( x, y, width, color="y" )
	plt.title('bonferroni_correction_dict', fontsize=18)
	plt.xlabel('Experiment type/no')
	plt.ylabel('Counts')

	plt.xticks(x + width/2.0, labels )
	
	plt.savefig(file_fig_n_tests_barplot)
	plt.close() # or plt.close(fig)
except Exception, e:
	print "Exception occurred during plotting: %s" % e # 






print "The script is complete!"

