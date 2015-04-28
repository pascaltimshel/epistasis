#!/usr/bin/env python

import os
import sys
import time

import argparse

import re
import collections

import glob

import resource # to allow enough filehandles to be opened

###################################### Library ######################################

import matplotlib
matplotlib.use('Agg') #Agg backend and not an X-using backend that required an X11 connection. Call use BEFORE importing pyplot!
# REF: http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np # for plotting an array

import epistasis_table_library
###################################### USAGE ######################################

#python epistasis_table_processing.py --path_main_input XXX


################## Broad ##################
##### hIMR90 #####

##### hESC #####

################## OSX ##################
##### TEST #####
# python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8_test_case/assigned

####### hIMR90 [5] #######
# JOB_DIR_NAME="hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME
# JOB_DIR_NAME="hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME
# JOB_DIR_NAME="hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME
# JOB_DIR_NAME="hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME
# JOB_DIR_NAME="hIMR90_width_1000_maf_5_q_1e-06_epi1_1e-10"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME

####### hESC [4] #######
# JOB_DIR_NAME="hESC_width_1000_maf_5_q_1e-12_epi1_1e-10"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME
# JOB_DIR_NAME="hESC_width_500_maf_5_q_1e-14_epi1_1e-8"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME
# JOB_DIR_NAME="hESC_width_500_maf_5_q_1e-16_epi1_1e-8"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME
# JOB_DIR_NAME="hESC_width_2500_maf_5_q_1e-13_epi1_1e-10"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME

####### lan-et-al_K562 [2] #######
# JOB_DIR_NAME="lan-et-al_K562_width_1000_maf_5_q_OUTLIER_RM_epi1_1e-8"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME
# JOB_DIR_NAME="lan-et-al_K562_width_5000_maf_5_q_OUTLIER_RM_epi1_1e-8"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME

####### contactCount_1 #######
# JOB_DIR_NAME="hESC-contactCount_1_width_1000_maf_5_q_1_epi1_1e-8"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME
# JOB_DIR_NAME="hIMR90-contactCount_1_width_1000_maf_5_q_1_epi1_1e-8"; python epistasis_table_processing.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/$JOB_DIR_NAME


###################################### SYNOPSIS ######################################


###################################### TODO ######################################


###################################### TESTs to make ######################################


###################################### FILE snippets ######################################
################## file_epistasis_table ##################
### Example ### 
# CHR_1	SNP_1	CHR_2	SNP_2	BETA	CHISQ	PVALUE	PHENOTYPE	EIID	EID
# 1	rs12564773	3	rs1514423	+0.10238	41.25541	1.33581E-10	ILMN_1653904	hic_1_1197	hic_1
# 1	rs12564755	3	rs1514423	+0.10238	41.25541	1.33581E-10	ILMN_1653904	hic_1_1197	hic_1
# 3	rs17009271	4	rs17676993	-0.41600	57.89248	2.76844E-14	ILMN_1659762	hic_1_1512	hic_1
### Expected format
# THIS IS A TAB SEPERATED FILE

################## file_experiments_list ##################
### Example
# hic_1	3.57003651219e-09
# null_1	3.22434260007e-09
# null_2	3.30130059359e-09
# null_3	3.39907068048e-09
# null_4	3.16754267171e-09
# null_5	3.49826535014e-09
### FORMAT: two column, tab seperated.


###################################### Batch job vs interactive ######################################
### QUESTION: "find out if running in shell or not (e.g. sun grid engine queue)"
# SOURCE #1: http://stackoverflow.com/questions/967369/python-find-out-if-running-in-shell-or-not-e-g-sun-grid-engine-queue
# SOURCE #2: http://stackoverflow.com/questions/1077113/how-do-i-detect-whether-sys-stdout-is-attached-to-terminal-or-not
if sys.stdout.isatty():
    print("Script is running in 'Interactive' mode, i.e. from a shell/terminal")
    flag_mode_interactive = True
else:
    print("Script is running in 'Non-interactive' mode, i.e. from batch queueing system")
    flag_mode_interactive = False

### Code dependent on the interactive mode ###
# 1) IN 'Interactive' mode: the prompt is enabled if "path_main_out" exists.


###################################### PARAMETERs ###################################### 
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--path_main_input", required=True, help="""
	*Main* path is where the *epistasis_table* lives.
	Expected folder organization in path_main_input:
	<path_main_input>/ (e.g. fastEpi_compiled/assigned)
			epistasis_table.txt												[*required*]
			probe_list.txt													[*required*]
			experiments_list.txt			[*required*]

	The program relies on the *specific name* in path_main_input to match the specified pattern.
	""")
args = arg_parser.parse_args()

path_main_input = os.path.abspath(args.path_main_input)

################## Set output filenames ##################
path_main_out = path_main_input + "/epistasis_table_processing" # dir
# extension_XXX = ".txt" # file extension name

### Epistasis: epistasis tables
file_epistasis_table_processed = path_main_out + "/epistasis_table_pruned_processed.txt"
file_epistasis_table_pruned_EIID = path_main_out + "/epistasis_table_pruned_EIID.txt"
file_epistasis_table_pruned_hemani = path_main_out + "/epistasis_table_pruned_hemani.txt"

### Epistasis: experiment counts and stats file
file_epistatic_counts = path_main_out + "/epistatic_counts.txt"
file_epistatic_counts_csv = path_main_out + "/epistatic_counts.csv"
file_epistatic_stats = path_main_out + "/stats_processing.txt" # CHANGE THE NAME

### Probes
file_probe_counts_csv = path_main_out + "/probe_counts.csv"


#file_epistatic_counts_significant = path_main_out + "/epistatic_counts_significant.txt"

################## Get input files ##################
### file_epistasis_table
file_epistasis_table = path_main_input + "/epistasis_table.txt"
if not os.path.exists(file_epistasis_table):
	raise Exception("file_epistasis_table does not exists: %s" % file_epistasis_table)

### Probe list
file_probe_list = path_main_input+"/probe_list.txt"
if not os.path.exists(file_probe_list):
	raise Exception("file_probe_list does not exists: %s" % file_probe_list)

### Experiments list
file_experiments_list = path_main_input+"/experiments_list.txt"
if not os.path.exists(file_experiments_list):
	raise Exception("file_experiments_list does not exists: %s" % file_experiments_list)



###################################### INITIAL Error checks ######################################

###################################### OUTPUT ######################################
if os.path.exists(path_main_out):
	print "WARNING: path_main_out={} exists".format(path_main_out)
	
	# if flag_mode_interactive: # only ask for user input if the script is running from command line
	# 	ans = ""
	# 	while ans != "yes":
	# 		ans = raw_input("Do you want overwrite the content (type 'yes' to continue)? ")
	# 	print # add extra newline

################## *Make OUTPUT dirs* ##################
for path in [path_main_out]:
	if not os.path.exists(path):
		os.makedirs(path)

###################################### Read file with experiments ######################################

experiments_dict = {}
print "Reading list of experiments: ".format(file_experiments_list)
with open(file_experiments_list, 'r') as fh:
	for line in fh:
		fields = line.strip().split()
		### Get and save identifier
		EID = fields[0]
		experiments_dict[EID] = 1
print "Read n={} experiments".format(len(experiments_dict))

###################################### Read file_probe_list ######################################
probes_dict = dict()
with open(file_probe_list, 'r') as fh:
	# DESCRIPTION of file_probe_list
	# with NO header
	# One illumina_probe_id (e.g. ILMN_2221) per line
	for line in fh:
		illumina_probe_id = line.strip()
		probes_dict[illumina_probe_id] = 1

### Checks and variable definition
n_probes = len(probes_dict.keys())
assert(n_probes == epistasis_table_library.N_PHENOTYPES_TESTED) # 9269

###################################### Get probe annotation ######################################

print "Will call function to get probe annotation"
df_probe_annotation = epistasis_table_library.get_probe_annotation() # index of df is Probe_Id (illumina_probe_id)


###################################### Main Loop ######################################

################## Dict for epistatic SNP-pairs counts ##################
#epistatic_counts_dict = collections.defaultdict(lambda: collections.defaultdict(int)) # e.g. epistatic_counts_dict['MASTER_KEY']['EID'] = INTEGER
	# the expression inside the parenthesis must be a 'callable'. That is why I use a 'lambda function'
epistatic_counts_dict_master_keys = ["count_significant", "count_significant_pruned_EIID", "count_significant_pruned_hemani"] # *OBS* this is used for defining the order of the subplot histogram
epistatic_counts_dict = dict() #collections.defaultdict(dict)
epistatic_counts_dict["count_significant"] = collections.defaultdict(int)
epistatic_counts_dict["count_significant_pruned_EIID"] = collections.defaultdict(set) # *OBS*: set()
epistatic_counts_dict["count_significant_pruned_hemani"] = collections.defaultdict(int)
	### MASTER KEYS IN DICT:
	# count_significant
	# count_significant_pruned_EIID
	# count_significant_pruned_hemani


################## Dict for epistasis_table_*.txt ##################
# EID = Experimental IDentifier. E.g. hic_1, null_1
# EIID = Experimental Interaction IDentifier. E.g. hic_1_1923, null_1_119
epistasis_table_dict = collections.defaultdict(list) # *OBS*: list
# 	# D[EID].append(line)
epistasis_table_pruned_EIID_dict = epistasis_table_library.makehash() # "hash"/auto-vivification
	# D[EID][EIID]['best_pvalue'] = pvalue
	# D[EID][EIID]['best_line'] = line
epistasis_table_pruned_hemani_dict = epistasis_table_library.makehash() # "hash"/auto-vivification
	# D[EID][illumina_probe_id][chr_pair]['best_pvalue'] = pvalue
	# D[EID][illumina_probe_id][chr_pair]['best_line'] = line
	# *OBS*: 
		# chr_pair is formed by sorting chrA and chrB. Example: "2:10" or "10:19"
		# this is IMPORTANT to represent a chromosme pair uniquely and unambiguously

################## Dict for distribution of significant probe counts ##################
### DataFrame for distribution of significant probe counts
df_probe_counts_master_keys = ["count_assigned_significant"]
df_probe_counts = pd.DataFrame(columns=df_probe_counts_master_keys, index=probes_dict.keys()) # initialize data frame with NaNs.
	# INDEX: illumina_probe_id, e.g. ILMN_2221
	# COLUMNS: df_probe_counts_master_keys
df_probe_counts = df_probe_counts.fillna(0) # fill the data with 0s rather than NaNs



### Count the number of lines in the file
n_lines_file_epistasis_table = sum(1 for line in open(file_epistasis_table)) 

time_start_loop = time.time()

### READ and WRITE line-by-line
with open(file_epistasis_table, 'r') as fh_compiled:
	next(fh_compiled) # SKIPPING THE FIRST LINE!
	for line_no, line in enumerate(fh_compiled, start=1):

		time_elapsed_loop = time.time() - time_start_loop # <type 'float'>
		if line_no % 5000 == 0: 
			print "Main loop | #{line_no}/#{n_lines} | {pct_complete:.2f} % done | {sec:.2f} sec [{min:.2f} min]".format(line_no=line_no, n_lines=n_lines_file_epistasis_table, pct_complete=(line_no/float(n_lines_file_epistasis_table))*100, sec=time_elapsed_loop, min=time_elapsed_loop/60)
			sys.stdout.flush()


		# CHR_1	SNP_1	CHR_2	SNP_2	BETA	CHISQ	PVALUE	PHENOTYPE	EIID	EID
		# 1	rs12564773	3	rs1514423	+0.10238	41.25541	1.33581E-10	ILMN_1653904	hic_1_1197	hic_1
		line = line.strip() # strip() without argument will remove BOTH leading and trailing whitespace.
		fields = line.split() 
		chr_A, snp_A, chr_B, snp_B, pvalue, illumina_probe_id, EIID, EID = fields[0], fields[1], fields[2], fields[3], fields[6], fields[7], fields[8], fields[9]


		################## Defining/converting variables ##################
		### Chromosome pair
		try:
			chr_A = int(chr_A)
			chr_B = int(chr_B)
		except ValueError as e:
			print "ERROR: could not convert chr_A [{}] or chr_B [{}] to int.".format(chr_A, chr_B)
			print "Exception: [{exception}]".format(exception=e)
			print "Will re-raise exception..."
			raise
		
		chr_pair = ":".join(map(str, sorted([chr_A, chr_B]))) # NUMERICALLY sorted, e.g. "9:21".
			#^^ REMEMBER: .join() ONLY works on a list with string elements.

		### pvalue
		pvalue = float(pvalue)


		### probes counting
		df_probe_counts.ix[illumina_probe_id, "count_assigned_significant"] += 1
		
		### get probe annoation ###
		#columns of interest:
			# (Probe_Id --> illumina_probe_id)
			# RefSeq_ID --> e.g. NM_017583.3, XM_373469.3 [OBS: Trim the dot (".") and the trailing number - this will allow to map genes via DAVID]
			# Symbol --> e.g. TRIM44, LOC387701
			# Chromosome --> 11, 10
			# Probe_Chr_Orientation --> +, -
			# Probe_Coordinates --> [35786070-35786119], [92811754-92811767:92811768-92811803]
		phenotype_annotation_dict = dict()
		phenotype_annotation_dict["RefSeq_ID"] = df_probe_annotation.ix[illumina_probe_id, "RefSeq_ID"]
		phenotype_annotation_dict["Symbol"] = df_probe_annotation.ix[illumina_probe_id, "Symbol"]
		phenotype_annotation_dict["Chromosome"] = int(df_probe_annotation.ix[illumina_probe_id, "Chromosome"]) # OBS: int()
		phenotype_annotation_dict["Probe_Chr_Orientation"] = df_probe_annotation.ix[illumina_probe_id, "Probe_Chr_Orientation"]
		phenotype_annotation_dict["Probe_Coordinates_start"] = int(df_probe_annotation.ix[illumina_probe_id, "Probe_Coordinates"].split("-")[0]) # OBS: int()
			# e.g. 92811754-92811767:92811768-92811803 --> 92811754
		phenotype_annotation_dict["Probe_Coordinates"] = df_probe_annotation.ix[illumina_probe_id, "Probe_Coordinates"]


		### epistasis
		epistatic_counts_dict['count_significant'][EID] += 1
		epistatic_counts_dict['count_significant_pruned_EIID'][EID].add(EIID)


		### Line output
		line_out = "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
			line,
			phenotype_annotation_dict["RefSeq_ID"],
			phenotype_annotation_dict["Symbol"],
			phenotype_annotation_dict["Chromosome"],
			phenotype_annotation_dict["Probe_Chr_Orientation"],
			phenotype_annotation_dict["Probe_Coordinates_start"],
			phenotype_annotation_dict["Probe_Coordinates"]
			)
		#line_out = line

		# epistasis_table_dict
		epistasis_table_dict[EID].append( line_out )
		# epistasis_table_pruned_EIID_dict
		if pvalue < epistasis_table_pruned_EIID_dict[EID][EIID]['best_pvalue']:
			# D[EID][EIID]['best_pvalue'] = pvalue
			epistasis_table_pruned_EIID_dict[EID][EIID]['best_pvalue'] = pvalue
			epistasis_table_pruned_EIID_dict[EID][EIID]['best_line'] = line_out
		# epistasis_table_pruned_hemani_dict
		if pvalue < epistasis_table_pruned_hemani_dict[EID][illumina_probe_id][chr_pair]['best_pvalue']:
			# D[EID][illumina_probe_id][chr_pair]['best_pvalue'] = pvalue
			epistasis_table_pruned_hemani_dict[EID][illumina_probe_id][chr_pair]['best_pvalue'] = pvalue
			epistasis_table_pruned_hemani_dict[EID][illumina_probe_id][chr_pair]['best_line'] = line_out

print "Done with main loop!"

###################################### PROCESS and WRITE epistasis_table_* content ######################################

################## POPULATE *epistatic_counts_dict* AND WRITE epistasis_table_* ##################
### epistasis_table_dict ###
with open(file_epistasis_table_processed, 'w') as fh:
	for EID in sorted(epistasis_table_dict, key=epistasis_table_library.function_sort_EID):
		list_of_line_out = epistasis_table_dict[EID]
		for line_out in list_of_line_out:
			fh.write( line_out + "\n" )

### epistasis_table_pruned_EIID_dict ###
## epistasis_table_pruned_EIID_dict[EID][EIID]['best_pvalue'] = pvalue
with open(file_epistasis_table_pruned_EIID, 'w') as fh:
	for EID in sorted(epistasis_table_pruned_EIID_dict, key=epistasis_table_library.function_sort_EID):
		EID_dict = epistasis_table_pruned_EIID_dict[EID]
			# EID_dict.keys() is EIID
			# EXAMPLE --> ['hic_1_2013', 'hic_1_2574']

		### The below can be used for validation/extra check ###
		##n_epistatic_counts = len(EID_dict[EID].keys()) # this is the number of EIID (EIID) that is significant for a given EID
		##epistatic_counts_dict["count_significant_pruned_EIID__alternative_counting_method"][EID] = n_epistatic_counts # no need for incrementing ("+=")

		for EIID in sorted(EID_dict, key=epistasis_table_library.function_sort_EIID):
			EIID_dict = EID_dict[EIID]
				# EIID_dict.keys()
				# EXAMPLE [*FIXED*] --> ['best_pvalue', 'best_line']
			line_out = EIID_dict["best_line"]
			fh.write( line_out + "\n" )


### epistasis_table_pruned_hemani_dict ###
## epistasis_table_pruned_hemani_dict[EID][illumina_probe_id][chr_pair]['best_pvalue'] = pvalue
#for EID, EID_dict in epistasis_table_pruned_hemani_dict.items():
with open(file_epistasis_table_pruned_hemani, 'w') as fh:
	for EID in sorted(epistasis_table_pruned_hemani_dict, key=epistasis_table_library.function_sort_EID):
		EID_dict = epistasis_table_pruned_hemani_dict[EID]
			# EID_dict.keys() is illumina_probe_id
			# EXAMPLE --> ['ILMN_1684445', 'ILMN_1789639', 'ILMN_1701551', 'ILMN_1663793', 'ILMN_1814688', 'ILMN_1775441']
		
		# for illumina_probe_id, illumina_probe_id_dict in EID_dict.items():
		for illumina_probe_id in EID_dict:
			illumina_probe_id_dict = EID_dict[illumina_probe_id]
				# illumina_probe_id_dict.keys() is chr_pair
				# EXAMPLE --> ['4:14', '4:6']
			n_chr_pair_epistatic_count = len(illumina_probe_id_dict.keys()) # this is the number of (unique) chromosome pairs with significant epistatic SNP-pairs
			epistatic_counts_dict["count_significant_pruned_hemani"][EID] += n_chr_pair_epistatic_count # populating master key in epistatic_counts_dict
			# ^OBS: we INCREMENT the count
			
			### DEBUGGGING ###
			# if n_chr_pair_epistatic_count > 1:
			# 	pdb.set_trace()
			# from pprint import pprint
			# illumina_probe_id_dict.keys()

			for chr_pair in sorted(illumina_probe_id_dict, key=lambda x: x.split(":")[0]): #sorted on first chromosome. OBS: we are not using .items().
				#POTENTIAL: you could also sort by best p-value for a given probe.
				#pdb.set_trace()
				line_out = illumina_probe_id_dict[chr_pair]["best_line"]
				fh.write( line_out + "\n" )


###################################### Writing MORE files ######################################

###################################### CREATING DataFrame: Calculating Empirical P-value ######################################

#df_experiment_identifier_counts = pd.DataFrame()
df_experiment_identifier_counts = pd.DataFrame(columns=epistatic_counts_dict.keys(), index=experiments_dict.keys()) # initialize data frame with NaNs.
df_experiment_identifier_counts = df_experiment_identifier_counts.fillna(0) # fill the data with 0s rather than NaNs
### COLUMNS IN DATAFRAME: "keys of epistatic_counts_dict"
	# count_significant
	# count_significant_pruned_EIID
	# count_significant_pruned_hemani


for key_master in epistatic_counts_dict:
	for EID in epistatic_counts_dict[key_master]:
	#for EID in experiments_dict:
		#df_experiment_identifier_counts[key_master][]
		if key_master == 'count_significant_pruned_EIID':
			n_epistatic_counts = len(epistatic_counts_dict[key_master][EID]) # length of set()
			### Validation: the two methods of counting should give the same
			assert( n_epistatic_counts == len(epistasis_table_pruned_EIID_dict[EID].keys()) )
			#if not n_epistatic_counts == len(epistasis_table_pruned_EIID_dict[EID].keys()):
				#pdb.set_trace()
			
			### assignment to data frame
			df_experiment_identifier_counts.ix[EID, key_master] = n_epistatic_counts
		else:
			df_experiment_identifier_counts.ix[EID, key_master] = epistatic_counts_dict[key_master][EID]
		#print key_master, EID
		#if 'c' in df_experiment_identifier_counts.index: pdb.set_trace()

################## Creating P-value dict ##################
### About p-values: "obtaining a result EQUAL TO or MORE EXTREME than what was actually observed" --> p_val = sum(X >= X_OBS)
p_value_dict = {} # KEYS will be: count_significant_pruned_EIID, count_significant, count_all
for key_master in epistatic_counts_dict:
	p_value_dict[key_master] = sum(df_experiment_identifier_counts.ix[:, key_master] >= df_experiment_identifier_counts.ix['hic_1', key_master])/float(len(df_experiment_identifier_counts))
### OLD:
# p_value_count_significant_pruned = sum(df_experiment_identifier_counts.ix[:, 'count_significant_pruned_EIID'] >= df_experiment_identifier_counts.ix['hic_1', 'count_significant_pruned_EIID'])/float(len(df_experiment_identifier_counts))
# p_value_count_significant = sum(df_experiment_identifier_counts.ix[:, 'count_significant'] >= df_experiment_identifier_counts.ix['hic_1', 'count_significant'])/float(len(df_experiment_identifier_counts))
# p_value_count_all = sum(df_experiment_identifier_counts.ix[:, 'count_all'] >= df_experiment_identifier_counts.ix['hic_1', 'count_all'])/float(len(df_experiment_identifier_counts))

###################################### WRITE CSV TO FILE: df_experiment_identifier_counts ######################################
### Sorting, inplace
#tmp_sort_order_list = ["count_significant_pruned_EIID", "count_significant"]
tmp_sort_order_list = ["count_significant_pruned_hemani", "count_significant"]
df_experiment_identifier_counts.sort(tmp_sort_order_list, ascending=False, inplace=True)
### Writing file
df_experiment_identifier_counts.to_csv(file_epistatic_counts_csv) # sep='\t', index=True, header=True


###################################### WRITE CSV TO FILE: df_probe_counts ######################################
### Sorting, inplace
tmp_sort_order_list = ["count_assigned_significant"]
df_probe_counts.sort(tmp_sort_order_list, ascending=False, inplace=True)
### Writing file
df_probe_counts.to_csv(file_probe_counts_csv) # sep='\t', index=True, header=True


###################################### Write stats file ######################################
################## Create ... ##################

################## Write ##################
with open(file_epistatic_stats, 'w') as fh:

	### Probes
	probe_stats_count_assigned_significant = sum(df_probe_counts["count_assigned_significant"] > 0)
	fh.write( "PROBE count_assigned_significant: {} ({:.2f} % of total number of probes in probes_dict)\n".format( probe_stats_count_assigned_significant, probe_stats_count_assigned_significant/float(n_probes)*100 ) )	
		

	### Hi-C
	fh.write( "HIC count_significant: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant']) )
	fh.write( "HIC count_significant_pruned_EIID: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant_pruned_EIID']) )
	fh.write( "HIC count_significant_pruned_hemani: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant_pruned_hemani']) )

	fh.write( "p_value_dict[count_significant]: {}\n".format(p_value_dict['count_significant']) )
	fh.write( "p_value_dict[count_significant_pruned_EIID]: {}\n".format(p_value_dict['count_significant_pruned_EIID']) )
	fh.write( "p_value_dict[count_significant_pruned_hemani]: {}\n".format(p_value_dict['count_significant_pruned_hemani']) )


###################################### Plotting ######################################


################## Histogram ##################
print "will plot histogram now..."
for column_name, series in df_experiment_identifier_counts.iteritems(): # Iterator over (column, series) pairs
	print "PLOTTING HISTOGRAM FOR COLUMN_NAME: {}".format(column_name)
	fig_filename = "{base}/fig_{plot}.{ext}".format(base=path_main_out, plot=column_name, ext="pdf")

	x = series
	plt.figure()
	ax = plt.gca()  # NEEDED for "transform=ax.transAxes"
	## OR --> fig, ax = plt.subplots(1)
	## OR --> fig = plt.figure(); ax = fig.add_subplot(111)
	_ = plt.hist(x, bins=500) # 100 looks ok # ---> The return value is a tuple (n, bins, patches) 

	### Adding VERTICAL LINE for OBSERVED statistic (HiC)
	xv1 = df_experiment_identifier_counts.ix['hic_1', column_name]
	plt.axvline(xv1, color='blue', linestyle='dashed', linewidth=2)  # SEE --> http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.axhline
	
	### Adding p-value text
	textstr = "p-value={p_value}".format(p_value=p_value_dict[column_name])
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5) # # these are matplotlib.patch.Patch properties
	# place a text box in upper left in axes coords
	plt.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props) # transform=ax.transAxes
	### Our Favorite Recipes - Placing text boxes
	# http://matplotlib.org/users/recipes.html#placing-text-boxes
	# http://matplotlib.org/users/transforms_tutorial.html#transforms-tutorial
	# http://matplotlib.org/api/text_api.html#matplotlib.text.Text

	### Adding titles and labels
	plt.title(column_name, fontsize=18)
	plt.xlabel(column_name)
	plt.ylabel('Counts')
	plt.savefig(fig_filename)
	plt.close() # or plt.close(fig)


################## Subplot histogram - EPISTASIS ENRICHMENT ##################
print "GENERATION PLOT: subplot EPISTASIS ENRICHMENT"
subplot_n_row = len(df_experiment_identifier_counts.columns)
print "subplot_n_row={}".format(subplot_n_row)
#for subplot_no, (column_name, series) in enumerate(df_experiment_identifier_counts.iteritems(), start=1): # Iterator over (column, series) pairs
for subplot_no, column_name in enumerate(epistatic_counts_dict_master_keys, start=1): # *OBS: looping a in specific order*
	series = df_experiment_identifier_counts[column_name]
	print "PLOTTING SUBPLOT HISTOGRAM FOR COLUMN_NAME: {}".format(column_name)
	print "subplot_no={}".format(subplot_no)

	x = series
	

	ax = plt.subplot(subplot_n_row, 1, subplot_no) # .subplot(nrow, ncol, plot_idx)
	### matplotlib.pyplot.subplot(*args, **kwargs)
	# Return a subplot axes positioned by the given grid definition.
	# Typical call signature:
	# subplot(nrows, ncols, plot_number)


	### OR ###
	# fig, ax = plt.subplots(10, 10)
	# where ax will contain one hundred axis in a list (of lists).
	#Definition: plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, **fig_kw)
	#Create a figure with a set of subplots already made.

	_ = plt.hist(x, bins=500) # 100 looks ok # ---> The return value is a tuple (n, bins, patches) 
	### CONSIDER: ax.hist()

	### Adding VERTICAL LINE for OBSERVED statistic (HiC)
	xv1 = df_experiment_identifier_counts.ix['hic_1', column_name]
	plt.axvline(xv1, color='blue', linestyle='dashed', linewidth=2)  # SEE --> http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.axhline


	### Adding p-value text
	textstr = "p-value={p_value:.3f}".format(p_value=p_value_dict[column_name])
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5) # # these are matplotlib.patch.Patch properties
	# place a text box in upper left in axes coords
	plt.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', horizontalalignment='right', bbox=props) # transform=ax.transAxes
	# plt.text(x, x, string)
	### Our Favorite Recipes - Placing text boxes
	# http://matplotlib.org/users/recipes.html#placing-text-boxes
	# http://matplotlib.org/users/transforms_tutorial.html#transforms-tutorial
	# http://matplotlib.org/api/text_api.html#matplotlib.text.Text

	### Adding axis labels and title
	plt.xlabel(column_name)
	plt.ylabel('Counts')
	plt.title(column_name, fontsize=18)

plt.tight_layout()
#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
## http://matplotlib.org/users/tight_layout_guide.html
## http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.tight_layout
## http://matplotlib.org/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels

fig_filename = "{base}/fig_{plot}.{ext}".format(base=path_main_out, plot="subplot_Hi-C_epistasis_enrichment", ext="pdf")
### saving fig
plt.savefig(fig_filename)
plt.close() # or plt.close(fig)


################## Subplot histogram - PROBES ##################
print "GENERATION PLOT: subplot PROBES"
subplot_n_row = len(df_probe_counts.columns)
print "subplot_n_row={}".format(subplot_n_row)
for flag_logscale in [True, False]:
	print "flag_logscale: {}".format(flag_logscale)
	#for subplot_no, (column_name, series) in enumerate(df_probe_counts.iteritems(), start=1): # Iterator over (column, series) pairs
	for subplot_no, column_name in enumerate(df_probe_counts_master_keys, start=1): # *OBS: looping a in specific order*
		series = df_probe_counts[column_name]
		print "PLOTTING SUBPLOT HISTOGRAM FOR COLUMN_NAME: {}".format(column_name)
		print "subplot_no={}".format(subplot_no)

		x = series
		ax = plt.subplot(subplot_n_row, 1, subplot_no) # .subplot(nrow, ncol, plot_idx)
		_ = plt.hist(x, bins=500) # 100 looks ok # ---> The return value is a tuple (n, bins, patches) 

		### Adding VERTICAL LINE
		xv1 = x.mean()
		x_min = x.min()
		x_max = x.max()
		plt.axvline(xv1, color='blue', linestyle='dashed', linewidth=2)

		### About log-scale plots ###
		#REF: http://stackoverflow.com/a/17952890
		#http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.yscale
		#The issue is with the bottom of bars being at y=0 and the default
		#is to mask out in-valid points (log(0) -> undefined) when doing
		#the log transformation (there was discussion of changing this,
		#but I don't remember which way it went) so when it tries to draw
		#the rectangles for you bar plot, the bottom edge is masked out ->
		#no rectangles.
		# ax.set_ylim(1, SOMETING)
		
		# nonposx/nonposy: ['mask' | 'clip' ]: non-positive values in x or y can be masked as invalid, or clipped to a very small positive number

		if flag_logscale:
			plt.yscale('log', basey=10) # what is the default --> log10.
			#plt.yscale('log', nonposy='clip')

		### Adding p-value text
		textstr = "mean={mean:.3f}\nmax={max:.1f}\nmin={min:.1f}".format(mean=xv1, max=x_max, min=x_min)
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5) # # these are matplotlib.patch.Patch properties
		# place a text box in upper left in axes coords
		plt.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', horizontalalignment='right', bbox=props) # transform=ax.transAxes

		### Adding axis labels and title
		plt.xlabel(column_name)
		plt.ylabel('Counts')
		plt.title(column_name, fontsize=18)

	plt.tight_layout()
	fig_filename = "{base}/fig_{plot}.{ext}".format(base=path_main_out, plot="subplot_probes_distribution_ylogscale-{}".format(flag_logscale), ext="pdf")
	### saving fig
	plt.savefig(fig_filename)
	plt.close() # or plt.close(fig)


###################################### Finishing ######################################
print "The script is complete!"
