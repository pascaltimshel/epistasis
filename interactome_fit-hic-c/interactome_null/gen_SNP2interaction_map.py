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

try:
	import cPickle as pickle # NO subclass - cPickle is written in C.
except:
	print "loaded pickle - could not load cPickle..."
	import pickle

import time

import memory_profiler; import psutil # OBS: the program RUN EXTREMELY SLOW when memory_profiling WITHOUT psutil. To install pip in homedir on a server: "pip install -user psutil"
	### See warning
	#/home/unix/ptimshel/.local/lib/python2.7/site-packages/memory_profiler.py:62: UserWarning: psutil module not found. memory_profiler will be slow
	#  warnings.warn("psutil module not found. memory_profiler will be slow")

import pdb

###################################### USAGE ######################################

#python gen_SNP2interaction_map.py --path_interaction_table XXX --file_null_table XXX --q_threshold XXX --path_main_out XXX

################## OSX ##################
### fastEpi_compiled/hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8
#python gen_SNP2interaction_map.py --path_interaction_table /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_snpsets/maf_5_sets/500_snppool_hIMR90_q_1e-08 --file_null_table /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/null_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-08.nperm_1000.txt --q_threshold 1e-08 --path_main_out /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8
# mprof run gen_SNP2interaction_map.py XXX

### fastEpi_compiled/hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8
#python gen_SNP2interaction_map.py --path_interaction_table /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_snpsets/maf_5_sets/500_snppool_hIMR90_q_1e-06 --file_null_table /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/null_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-06.nperm_1000.txt --q_threshold 1e-06 --path_main_out /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8


### hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10 (could not complete - TOO MUCH MEMORY?)
#python gen_SNP2interaction_map.py --path_interaction_table /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_snpsets/maf_5_sets/50000_snppool_hIMR90_q_1e-09 --file_null_table /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/null_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-09.nperm_1000.txt --q_threshold 1e-09 --path_main_out /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10
	# at null_192_1017 --> Segmentation fault: 11 (likely because it was using 6 GB memory)


################## Broad ##################

####### hIMR90 #######
### hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10
#python gen_SNP2interaction_map.py --path_interaction_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/50000_snppool_hIMR90_q_1e-09 --file_null_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/null_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-09.nperm_1000.txt --q_threshold 1e-09 --path_main_out /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10/fastEpi_compiled

### hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8
#python gen_SNP2interaction_map.py --path_interaction_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/500_snppool_hIMR90_q_1e-08 --file_null_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/null_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-08.nperm_1000.txt --q_threshold 1e-08 --path_main_out /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8/fastEpi_compiled

### hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8
#python gen_SNP2interaction_map.py --path_interaction_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/500_snppool_hIMR90_q_1e-06 --file_null_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/null_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-06.nperm_1000.txt --q_threshold 1e-06 --path_main_out /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8/fastEpi_compiled

### hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8
#python gen_SNP2interaction_map.py --path_interaction_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/2500_snppool_hIMR90_q_1e-07 --file_null_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/null_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-07.nperm_1000.txt --q_threshold 1e-07 --path_main_out /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8/fastEpi_compiled

####### hESC #######
### hESC_width_1000_maf_5_q_1e-12_epi1_1e-10
#python gen_SNP2interaction_map.py --path_interaction_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/1000_snppool_hESC_q_1e-12 --file_null_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/null_table.fit-hi-c.nosex.interchromosomal.hESC.q_1e-12.nperm_1000.txt --q_threshold 1e-12 --path_main_out /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hESC_width_1000_maf_5_q_1e-12_epi1_1e-10/cat_epistasis_150210_125525
#bsub -J gen_map_40GB_MEDPOP -q MEDPOP -R 'rusage[mem=40]' -o fastEpi_compiled_MEDPOP_40GB.bsub.out python gen_SNP2interaction_map.py --path_interaction_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/1000_snppool_hESC_q_1e-12 --file_null_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/null_table.fit-hi-c.nosex.interchromosomal.hESC.q_1e-12.nperm_1000.txt --q_threshold 1e-12 --path_main_out /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hESC_width_1000_maf_5_q_1e-12_epi1_1e-10/fastEpi_compiled_MEDPOP_40GB
#bsub -J gen_map_40GB_week -q week -R 'rusage[mem=40]' -o fastEpi_compiled_week_40GB.bsub.out python gen_SNP2interaction_map.py --path_interaction_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/1000_snppool_hESC_q_1e-12 --file_null_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/null_table.fit-hi-c.nosex.interchromosomal.hESC.q_1e-12.nperm_1000.txt --q_threshold 1e-12 --path_main_out /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hESC_width_1000_maf_5_q_1e-12_epi1_1e-10/fastEpi_compiled_week_40GB
#bsub -J gen_map_20GB_week -q week -R 'rusage[mem=20]' -o fastEpi_compiled_week_20GB.bsub.out python gen_SNP2interaction_map.py --path_interaction_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/1000_snppool_hESC_q_1e-12 --file_null_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/null_table.fit-hi-c.nosex.interchromosomal.hESC.q_1e-12.nperm_1000.txt --q_threshold 1e-12 --path_main_out /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hESC_width_1000_maf_5_q_1e-12_epi1_1e-10/fastEpi_compiled_week_20GB

#bsub -J gen_map_40GB -q priority -R 'rusage[mem=40]' -o fastEpi_compiled_priority_40GB.bsub.out python gen_SNP2interaction_map.py --path_interaction_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/1000_snppool_hESC_q_1e-12 --file_null_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/null_table.fit-hi-c.nosex.interchromosomal.hESC.q_1e-12.nperm_1000.txt --q_threshold 1e-12 --path_main_out /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hESC_width_1000_maf_5_q_1e-12_epi1_1e-10/fastEpi_compiled_priority_40GB
#bsub -J gen_map_10GB -q priority -R 'rusage[mem=10]' -o /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hESC_width_1000_maf_5_q_1e-12_epi1_1e-10/fastEpi_compiled_priority_40GB.bsub.out python gen_SNP2interaction_map.py --path_interaction_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/1000_snppool_hESC_q_1e-12 --file_null_table /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/null_table.fit-hi-c.nosex.interchromosomal.hESC.q_1e-12.nperm_1000.txt --q_threshold 1e-12 --path_main_out /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hESC_width_1000_maf_5_q_1e-12_epi1_1e-10/fastEpi_compiled_priority_10GB

###################################### SYNOPSIS ######################################
# ???

### Output files
# ??

###################################### TODO ######################################
# 1) Count the number of interactions for each SNP. (Plot histogram?)
# 2) Think about an improved path to write files to: where should the files go?
# 3) 
# (consider writing out only interaction_identifiers were BOTH snps_A and snps_B are non-empty)

### What we know:
# 1) The histogram over experiments are uniform: this means that "hic" and "null" have the same number of tests, and thus SNPs, in them

###################################### TESTs to make ######################################
# OK: check q_threshold match
# OK: check equal length of null and df_interaction_table
# check header of df_interaction_table: interaction_ID	chr_A	pos_A	chr_B	pos_B	setA_size	setB_size	snps_A	snps_B


###################################### FILE snippets ######################################
################## file_null_table ##################
### Example ### 
# null_1	null_2	null_3	null_4	null_5	null_6	null_7	null_8	... LINES CONTINUE ...
# 921	3140	55	1961	579	2484	1644	2225	2871	511	... LINES CONTINUE ...
# 1291	1236	1770	2322	2281	72	1997	2259	2147	... LINES CONTINUE ...
# 1986	1975	2361	1574	1998	636	1189	3392	3069	... LINES CONTINUE ...
# 3147	401	725	650	607	814	1618	2510	2293	1726	3468	... LINES CONTINUE ...
# 699	838	2342	1261	2277	2453	2378	767	2130	1340	... LINES CONTINUE ...
# 3112	496	367	94	1944	2589	2829	1363	2430	1357	... LINES CONTINUE ...
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
arg_parser.add_argument("--path_main_out", required=True, help="Main path to write output files. Path will be created if it does not exists. PASCAL RECOMMENDS USING the respective 'fastEpi_compiled' path for this purpose.")
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
else:
	print "path_main_out={} exists".format(path_main_out)
	ans = ""
	while ans != "yes":
		ans = raw_input("Do you want overwrite the content (type 'yes' to continue)? ")
	print # add extra newline

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

################## Print messages ##################
print "len(df_interaction_table_snps): {}".format(len(df_interaction_table_snps))
print "len(df_null_table): {}".format(len(df_null_table))

################## Tests ##################
### Check Equal Length Of Data Frames
# The data frames must/should have equal length (rows) because they contain the same number of interactions (rows)
assert len(df_interaction_table_snps) == len(df_null_table)

assert list(df_interaction_table_snps.columns.values) == ["interaction_ID", "chr_A", "pos_A", "chr_B", "pos_B", "setA_size", "setB_size", "snps_A", "snps_B"]

###################################### Initialize container ######################################
### Types of container:
# 1) Defaultdict + set
# 2) Defaultdict + OrderedDict
# 3) Deque

###################################### MAIN LOOP ######################################

#@memory_profiler.profile
def populate_snp2interaction_dict():
	SNP2interaction_dict = collections.defaultdict(set) # using SET()
	#SNP2interaction_dict = collections.defaultdict(collections.OrderedDict) # OrderedDict is not a lot slower than a plain dict, but at least doubles the memory.
	#SNP2interaction_dict = {}

	bonferroni_correction_dict = collections.defaultdict(int)
	#SNP_interaction_count_dict = {}
	
	for column_no, (column_name, series) in enumerate(df_null_table.iteritems(), start=1): # Iterator over (column, series) pairs
		#print "COLUMN_NAME: {}".format(column_name)

		time_start = time.time()
		series_corrected_offset = series-1 # We need to subtract one from the index because Pandas/Python is zero-based and R is one-based
		for interaction_no, elem in enumerate(series_corrected_offset, start=1):
			## Note that interaction_no starts at one! This is because we want the interaction IDs to run from 1...N_interactions.
			
			### This is the index to access rows/interaction in Pandas Dataframes
			interaction_idx = interaction_no - 1

			snps = list() # Default value - empty list. Will be evaluated as False in boolean context.
			
			#if interaction_no==100: break

			### KEEP THIS - THERE ARE SOME USEFUL INFORMATION IN THE CODE ### 
			### OLD METHOD WHICH RESULTED IN UNEQUAL REPRESENTATION OF EXPERIMENTS_IDENTIFERS ACROSS SNPS
			### That is, some EXPERIMENTS_IDENTIFERS (e.g. null_617) would not be represented 3 times in a SNP (rs1322573) like the majority of EXPERIMENTS_IDENTIFERS 
			# try:
			# 	snps_A = df_interaction_table_snps["snps_A"][interaction_idx].split(";")
			# 	snps_B = df_interaction_table_snps["snps_B"][elem].split(";")
			# 	snps = snps_A + snps_B # similar to .extend(). OBS: this list could potentially contain DUPLICATES. However, that should not be a problem
			# 	#snps = snps_B
			# except AttributeError: # E.g. "AttributeError 'float' object has no attribute 'split'" when splitting "nan" ("nan" is type "float")")
			# 	pass
			# 	# --> Got AttributeError. Value is likely 'nan'. Will do NOTHING because this is ok! 


			try:
				snps_A = df_interaction_table_snps["snps_A"][interaction_idx].split(";")
			except AttributeError: # E.g. "AttributeError 'float' object has no attribute 'split'" when splitting "nan" ("nan" is type "float")")
				snps_A = []

			try:
				snps_B = df_interaction_table_snps["snps_B"][elem].split(";")
			except AttributeError: # E.g. "AttributeError 'float' object has no attribute 'split'" when splitting "nan" ("nan" is type "float")")
				snps_B = []

			snps = snps_A + snps_B # similar to .extend(). OBS: this list could potentially contain DUPLICATES. However, that should not be a problem

			#print "Got {} SNPs in 'snps'".format(len(snps))
			
			### Constructing interaction_identifier
			# column_name example: hic_1, null_1, null_2, ...
			experiment_type = column_name.split("_")[0] # {hic, null}
			experiment_no = column_name.split("_")[1] # {1, 2, .., N_perm}
			interaction_identifier = "{experiment_type}_{experiment_no}_{interaction_no}".format(experiment_type=experiment_type, experiment_no=experiment_no, interaction_no=interaction_no)
			#print interaction_identifier

			### DEBUG			
			# if interaction_identifier.endswith("_394"):
			# 	print interaction_identifier
			# 	pdb.set_trace()


			for snp in snps: 
				### Remarks:
				# "snps" may contain duplicates if snps_A and snps_B overlap because we are USING LIST PLUS (+) OPERATOR: snps_A + snps_B
				# duplicates is not a problem for this code. Just think about it.
				# duplicates should not exists when using INTER-CHROMOSOMAL interactions

				#if snp == 'rs1322573': pdb.set_trace() # ***DEBUGGING
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

		time_elapsed = time.time() - time_start # time nested loop
		time_elapsed_populate_dict = time.time() - time_start_populate_dict # time main loop
		print "Elapsed time for one experiment: {} s ({} min)".format( time_elapsed, time_elapsed/60 )
		print "Main loop | COLUMN_NAME: {identifer} | Experiment: #{status_current}/#{status_total} | {pct_complete:.2f} % done | {sec:.2f} sec [{min:.2f} min]".format(identifer=column_name, status_current=column_no, status_total=n_experiments, pct_complete=(column_no/float(n_experiments))*100, sec=time_elapsed_populate_dict, min=time_elapsed_populate_dict/60)
		
		#return (SNP2interaction_dict, bonferroni_correction_dict) # FOR TESTING QUICKLY!
	return (SNP2interaction_dict, bonferroni_correction_dict)

time_start_populate_dict = time.time()
n_experiments = len(df_null_table.columns) # *OBS*
(SNP2interaction_dict, bonferroni_correction_dict) = populate_snp2interaction_dict()
time_elapsed_populate_dict = time.time() - time_start_populate_dict
print "Elapsed time populate_snp2interaction_dict(): {} s ({} min)".format( time_elapsed_populate_dict, time_elapsed_populate_dict/60 )

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
print "Done"


print "Writing file_SNP2interaction_map. FORMAT=txt"
################## file_SNP2interaction_map ##################
with open(file_SNP2interaction_map, 'w') as fh_map:
	for snp in sorted(SNP2interaction_dict): # .items() does not work because we are using sorted() 
		set_string = ";".join(sorted(SNP2interaction_dict[snp])) # USING set() sorted set
		#set_string = ";".join(SNP2interaction_dict[snp]) # # USING OrderedDict()
		fh_map.write( "{}\t{}\t{}\n".format(snp, len(SNP2interaction_dict[snp]), set_string) )
print "Done"


###################################### Alternative methods for writing files ######################################

################## Pickle ##################
#ans = raw_input("'DEBUG mode': the program is about to write pickle file. Pres <Enter> to continue.")
#print "Got answer: [{}]".format(ans)

print "Writing file_SNP2interaction_map. FORMAT=pickle"
file_SNP2interaction_map_pickle = path_main_out + "/SNP2interaction_map.pickle" # add information about origin of file
with open(file_SNP2interaction_map_pickle, 'wb') as f:
	pickle.dump(SNP2interaction_dict, f, protocol=2) # pickle.dump(obj, file[, protocol]). pickle.HIGHEST_PROTOCOL
print "Done"


### File size
#SNP2interaction_map_pickle: 74 KB
#SNP2interaction_map.txt: 42 KB

################## JSON ##################

# import json
# OBS: Python sets are not json serializable: http://stackoverflow.com/questions/8230315/python-sets-are-not-json-serializable

# file_SNP2interaction_map_json = path_main_out + "/SNP2interaction_map.json" # add information about origin of file

# with open(file_SNP2interaction_map_json, 'w') as f:
# 	json.dump(SNP2interaction_dict, f)
# 	#Alternative: json.dump(obj, fp, indent=None, separators=2, sort_keys=True)


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

