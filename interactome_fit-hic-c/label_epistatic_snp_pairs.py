#!/usr/bin/env python

import os
import sys
import time

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

import glob

import resource # to allow enough filehandles to be opened

import pdb

###################################### CONSTANTS ######################################

N_PHENOTYPES_TESTED = 9269
ALPHA = 0.05

###################################### Functions ######################################

def get_blacklisted_snps():
	""" 
	Function to define and retrieve a set of blacklisted SNPs. SNPs can be blacklisted if they are e.g. duplicates from the .bim file.

	### BLACKLIST LOG ###
	04/23/2015: Added 15 unique SNPs from running "get_duplicates.py" on ~/Dropbox/5_Data/EGCUT_DATA/geno/all_clean/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.duplicates_unique_rsID.txt
	
	"""

	blacklisted_snps_set = set([
						"rs4935071",
						"rs1072594",
						"rs206146",
						"rs7383287",
						"rs2581",
						"rs9460309",
						"rs241448",
						"rs2523608",
						"rs11244",
						"rs2289150",
						"rs2070121",
						"rs1042337",
						"rs2249255",
						"rs241447"
						]) # len() --> 15



	return blacklisted_snps_set

def makehash():
	""" Function to create a perl-like hash. That is, a dict with autovivification """
	return collections.defaultdict(makehash) 

def function_sort_EID(EID):
	""" 
	Function to sort EID (experiment_identifier) for dicts. This is a UTILITY FUNCTION for sorted(a_dict, key=function_sort_EID)
	EID type: string.
	EID format example(s): hic_1, null_1, null_2

	NOTICE that hic_1 and null_1 will get equal ranking because both return 1.
	"""
	experiment = EID.split("_")[0] # "hic" or "null"
	
	### We set a numeric "code" for the experiments, so "hic" has HIGHER PRECEDENCE for sorting.
	if experiment == "hic":
		experiment_ranking = 1
	elif experiment == "null":
		experiment_ranking = 2
	else:
		raise Exception( "Received unexpected experiment: {}. The function only expects experiments 'hic' or 'null'".format(experiment) )
	### 
	experiment_no = int(EID.split("_")[1]) # e.g. "1", "15" or "293" [range is 1...n_perm].
	
	tuple_to_sort_by = (experiment, experiment_no)
	# SEE: Python secondary sorting: http://stackoverflow.com/a/16193637
	return tuple_to_sort_by

def function_sort_EIID(EIID):
	"""
	Function to sort EIID.
	EIID type: string.
	EIID format example(s): hic_1_2013, hic_1_2574, null_2_1074
	
	See function_sort_EID() for additional documentation.
	"""
	interaction_no = int(EIID.split("_")[-1]) # e.g. "2013" [range is 1..n_interactions]
	return interaction_no

###################################### USAGE ######################################

#python label_epistatic_snp_pairs.py --path_main_input XXX


################## Broad ##################
##### hIMR90 #####
### hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10
#python label_epistatic_snp_pairs.py --path_main_input /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10/fastEpi_compiled
	# ---> runtime: 10-40 min

### hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8
#python label_epistatic_snp_pairs.py --path_main_input /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8/fastEpi_compiled

### hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8
#python label_epistatic_snp_pairs.py --path_main_input /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8/fastEpi_compiled

### hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8
#python label_epistatic_snp_pairs.py --path_main_input /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8/fastEpi_compiled

##### hESC #####
### hESC_width_1000_maf_5_q_1e-12_epi1_1e-10
#python label_epistatic_snp_pairs.py --path_main_input /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hESC_width_1000_maf_5_q_1e-12_epi1_1e-10/fastEpi_compiled

### hESC_width_500_maf_5_q_1e-14_epi1_1e-8
#python label_epistatic_snp_pairs.py --path_main_input /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hESC_width_500_maf_5_q_1e-14_epi1_1e-8/fastEpi_compiled

### hESC_width_500_maf_5_q_1e-16_epi1_1e-8
#python label_epistatic_snp_pairs.py --path_main_input /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hESC_width_500_maf_5_q_1e-16_epi1_1e-8/fastEpi_compiled

### hESC_width_2500_maf_5_q_1e-13_epi1_1e-10
#python label_epistatic_snp_pairs.py --path_main_input /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hESC_width_2500_maf_5_q_1e-13_epi1_1e-10/fastEpi_compiled



################## OSX ##################
### hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8
#python label_epistatic_snp_pairs.py --path_main_input /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8

### hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8
#python label_epistatic_snp_pairs.py --path_main_input /Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8



###################################### SYNOPSIS ######################################
# ???

### Output files
# ??


###################################### TODO ######################################

### Analysis to run
# 1) Run label_epistatic_snp_pairs.py on 
	# a) hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8

### Add-on for the script
# 1) Calculate empirical p-value based on different criteria:
	# a) No filter. Do not consider bonferroni correction
	# b) Count only significant interactions based on bonferroni calculations
	# c) Same as b) but "prune" to keep only one pair from each "experiment_interaction_identifer"
# 2) Generate (scripted) histogram for empirical p-value. Indicate position of observed hic statistic

### Additinal
# x) Run on IMR90


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
# THIS IS A TAB SEPERATED FILE

################## file_SNP2interaction_map ##################
### FORMAT: pickled! (or perhaps text)

################## file_bonferroni_correction ##################
### Example
# hic_1   4609175
# null_1  4486950
# null_10 4511100
# null_100        4541756
### FORMAT: two column, tab seperated. [BTW: the file is alphabetically sorted by the first column]

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
	*Main* path where the FastEpistasis *COMPILED* results lives.
	Expected folder organization in path_main_input:
	<path_main_input>/ (e.g. fastEpi_compiled)
			XXX___probe_epistatic_count.txt
			XXX___results.epi.qt.lm.combined.txt	[*required*]
			XXX___stat_file.txt
			SNP2interaction_map/					[*required*]
				SNP2interaction_map.pickle 			[*required*]
				bonferroni_correction.txt 			[*required*] - used for reading the "experiments" (i.e. get the names of the null and hic)
				MORE_FILES...

	The program relies on the *specific names* in path_main_input (fastEpi_compiled) to match the specified pattern.
	""")
args = arg_parser.parse_args()

path_main_input = os.path.abspath(args.path_main_input)

################## Set output filenames ##################
path_main_out = path_main_input + "/assigned" # dir # hic_hypothesis_testing
path_out_epistatic_results = path_main_out + "/epistatic_results" # dir
extension_epistatic_results = "lm.combined.txt" # file extension name | e.g. XXX/fastEpi_compiled/assigned/epistatic_results/null_996.lm.combined.txt

### Epistasis: epistasis tables
file_epistasis_table = path_main_out + "/epistasis_table.txt"
file_epistasis_table_pruned_EIID = path_main_out + "/epistasis_table_pruned_EIID.txt"
file_epistasis_table_pruned_hemani = path_main_out + "/epistasis_table_pruned_hemani.txt"

### Epistasis: experiment counts and stats file
file_epistatic_counts = path_main_out + "/epistatic_counts.txt"
file_epistatic_counts_csv = path_main_out + "/epistatic_counts.csv"
file_epistatic_stats = path_main_out + "/epistatic_stats.txt"

### Probes
file_probe_counts_csv = path_main_out + "/probe_counts.csv"

file_null_false_negatives = path_main_out + "/null_false_negatives.txt" # 
file_epistatic_intrachromosomal = path_main_out + "/epistatic_intrachromosomal.txt" # 

#file_epistatic_counts_significant = path_main_out + "/epistatic_counts_significant.txt"

################## Get input files ##################
### Get FastEpistasis compiled result file - GLOBBING ON FILENAME PATTERN
glob_lm_pattern = path_main_input+"/*lm.combined.txt"
glob_lm_file = glob.glob(glob_lm_pattern)
if not len(glob_lm_file) == 1:
	print "glob_lm_pattern: {}".format(glob_lm_pattern)
	raise Exception( "glob_lm_file does not contain EXACTLY one matching file. Matches are:\n[{}]".format("\n".join(glob_lm_file)) )
file_fastepistasis_lm_combined = glob_lm_file[0] # because of the above check, we know that the list only contains one element

### Probe list
glob_probe_list_pattern = path_main_input+"/*probe_list.txt"
glob_probe_list_file = glob.glob(glob_probe_list_pattern)
assert(len(glob_probe_list_file)==1) # error check
file_probe_list = glob_probe_list_file[0] # because of the above check, we know that the list only contains one element

### Maps and bonferroni correction
file_SNP2interaction_map = path_main_input + "/SNP2interaction_map/SNP2interaction_map.pickle"
file_bonferroni_correction = path_main_input + "/SNP2interaction_map/bonferroni_correction.txt"

###################################### INITIAL Error checks ######################################

################## Check that files exists ##################
if not os.path.exists(file_SNP2interaction_map):
	raise Exception("file_SNP2interaction_map does not exists: %s" % file_SNP2interaction_map)

###################################### OUTPUT ######################################
if os.path.exists(path_main_out):
	print "WARNING: path_main_out={} exists".format(path_main_out)
	
	if flag_mode_interactive: # only ask for user input if the script is running from command line
		ans = ""
		while ans != "yes":
			ans = raw_input("Do you want overwrite the content (type 'yes' to continue)? ")
		print # add extra newline

################## *Make OUTPUT dirs* ##################
for path in [path_main_out, path_out_epistatic_results]:
	if not os.path.exists(path):
		os.makedirs(path)

################## Clearing any existing files that will be *APPENDED* to ##################

for filename in [file_epistatic_intrachromosomal, file_null_false_negatives]:
	if os.path.exists(filename):
		print "Removing existing file: {}".format(filename)
		os.remove(filename)

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
assert(n_probes == N_PHENOTYPES_TESTED) # 9269

###################################### Read/Load input filesfile_SNP2interaction_map ######################################
filesize_SNP2interaction_map = os.path.getsize(file_SNP2interaction_map) >> 20 # BitWise Operations: shift-right | x >> y :Returns x with the bits shifted to the right by y places.
print "Loading SNP2interaction_dict via pickled file (size={} MB) ...".format(filesize_SNP2interaction_map)
time_start_tmp = time.time()
with open(file_SNP2interaction_map, 'r') as fh: # perhaps read the pickle file in binary mode.
	SNP2interaction_dict = pickle.load(fh)
print "Done loading pickled file in {:.2f} min".format( (time.time()-time_start_tmp)/60 )

###################################### Get *BLACKLISTED SNPs* ######################################

blacklisted_snps_set = get_blacklisted_snps() # returns set()
print "Retrieved blacklisted_snps_set: n={}".format(len(blacklisted_snps_set))

###################################### Open files ######################################

#resource.getrlimit(resource.RLIMIT_NOFILE)
#resource.setrlimit(resource.RLIMIT_NOFILE, (500,-1))

soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
print "SYSTEM INFO: number of files this process can open: (resource.RLIMIT_NOFILE, soft limit): {}".format(soft)
print "SYSTEM INFO: number of files this process can open: (resource.RLIMIT_NOFILE, hard limit): {}".format(hard)
resource.setrlimit(resource.RLIMIT_NOFILE, (4000, hard)) # *CHANGING SOFT LIMIT* to 9999 (abitrary high number) # Broad RHEL ISH has a hardlimit of 4096
soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
print "SYSTEM INFO - UPDATED: number of files this process can open: (resource.RLIMIT_NOFILE, soft limit): {}".format(soft)

################## Open filehandles ##################
experiment_identifier_fh_dict = {}
bonferroni_correction_dict = {}

print "Will now open filehandles in writing mode for assigning results.."
with open(file_bonferroni_correction, 'r') as fh:
	for line in fh:
		fields = line.strip().split()
		
		### Get identifier
		experiment_identifier = fields[0]

		### Set bonferroni correction - *IMPORTANT*
		bonferroni_correction_dict[experiment_identifier] = ALPHA/(int(fields[1])*N_PHENOTYPES_TESTED) # <-- FLOAT. alpha/(N_epistasis_tests*N_phenotypes_tests)

		### Set filename
		tmp_filename = "{path}/{basename}.{ext}".format(path=path_out_epistatic_results, basename=experiment_identifier, ext=extension_epistatic_results)
		tmp_fh = open(tmp_filename, 'w')

		experiment_identifier_fh_dict[experiment_identifier] = tmp_fh # e.g. experiment_identifier_fh_dict['null_1'] = <FILE_HANDLE>

print "Opened {} filehandles".format(len(experiment_identifier_fh_dict))

################## Open *ADDITONAL* filehandles ##################

### Filenames
file_multiple_assignments = "{path}/{basename}.{ext}".format(path=path_out_epistatic_results, basename="multiple_assigments", ext=extension_epistatic_results)
file_unassigned = "{path}/{basename}.{ext}".format(path=path_out_epistatic_results, basename="unassigned", ext=extension_epistatic_results)
### Filehandles
fh_multiple_assignments = open(file_multiple_assignments, 'w')
fh_unassigned = open(file_unassigned, 'w')

###################################### Main Loop ######################################

### Counters
count_total = 0
count_unassigned = 0
count_assigned = 0
count_multiple_assignments = 0
count_null_false_positives = 0
count_intrachromosomal_interactions = 0

### Counters/dict
counter_snp_not_in_SNP2interaction = collections.Counter()
counter_snp_blacklisted = collections.Counter()

################## Dict for epistatic SNP-pairs counts ##################
#epistatic_counts_dict = collections.defaultdict(lambda: collections.defaultdict(int)) # e.g. epistatic_counts_dict['MASTER_KEY']['EXPERIMENT_IDENTIFIER'] = INTEGER
	# the expression inside the parenthesis must be a 'callable'. That is why I use a 'lambda function'
epistatic_counts_dict_master_keys = ["count_all", "count_significant", "count_significant_pruned_EIID", "count_significant_pruned_hemani"] # *OBS* this is used for defining the order of the subplot histogram
epistatic_counts_dict = dict() #collections.defaultdict(dict)
epistatic_counts_dict["count_all"] = collections.defaultdict(int)
epistatic_counts_dict["count_significant"] = collections.defaultdict(int)
epistatic_counts_dict["count_significant_pruned_EIID"] = collections.defaultdict(set) # *OBS*: set()
epistatic_counts_dict["count_significant_pruned_hemani"] = collections.defaultdict(int)
	### MASTER KEYS IN DICT:
	# count_all
	# count_significant
	# count_significant_pruned_EIID
	# count_significant_pruned_hemani


################## Dict for epistasis_table_*.txt ##################
# EID = Experimental IDentifier. E.g. hic_1, null_1
# EIID = Experimental Interaction IDentifier. E.g. hic_1_1923, null_1_119
epistasis_table_dict = collections.defaultdict(list) # *OBS*: list
	# D[EID].append(line)
epistasis_table_pruned_EIID_dict = makehash() # "hash"/auto-vivification
	# D[EID][EIID]['best_pvalue'] = pvalue
	# D[EID][EIID]['best_line'] = line
epistasis_table_pruned_hemani_dict = makehash() # "hash"/auto-vivification
	# D[EID][illumina_probe_id][chr_pair]['best_pvalue'] = pvalue
	# D[EID][illumina_probe_id][chr_pair]['best_line'] = line
	# *OBS*: 
		# chr_pair is formed by sorting chrA and chrB. Example: "2:10" or "10:19"
		# this is IMPORTANT to represent a chromosme pair uniquely and unambiguously




################## Dict for distribution of significant probe counts ##################
# probe_counts_dict = dict()
# probe_counts_dict['count_all'] = collections.defaultdict(int)
# probe_counts_dict['count_assigned'] = collections.defaultdict(int)
# probe_counts_dict['count_assigned_significant'] = collections.defaultdict(int)
### DataFrame for distribution of significant probe counts
df_probe_counts_master_keys = ["count_all", "count_assigned", "count_assigned_significant"]
df_probe_counts = pd.DataFrame(columns=df_probe_counts_master_keys, index=probes_dict.keys()) # initialize data frame with NaNs.
	# INDEX: illumina_probe_id, e.g. ILMN_2221
	# COLUMNS: df_probe_counts_master_keys
df_probe_counts = df_probe_counts.fillna(0) # fill the data with 0s rather than NaNs





### Count the number of lines in the file
n_lines_fastepistasis_lm_combined = sum(1 for line in open(file_fastepistasis_lm_combined)) # this is pretty fast (only a bit slower than "wc -l"). 3.5 sec for a 1 GB file

time_start_loop = time.time()

### READ and WRITE line-by-line
with open(file_fastepistasis_lm_combined, 'r') as fh_compiled:
	next(fh_compiled) # SKIPPING THE FIRST LINE!
	for line_no, line in enumerate(fh_compiled, start=1):
		count_total += 1

		time_elapsed_loop = time.time() - time_start_loop # <type 'float'>
		if line_no % 5000 == 0: 
			print "Main loop | #{line_no}/#{n_lines} | {pct_complete:.2f} % done | {sec:.2f} sec [{min:.2f} min]".format(line_no=line_no, n_lines=n_lines_fastepistasis_lm_combined, pct_complete=(line_no/float(n_lines_fastepistasis_lm_combined))*100, sec=time_elapsed_loop, min=time_elapsed_loop/60)
			sys.stdout.flush()


		#chr1, snp_A, chr2, snp_B, beta, chisq, pvalue, phenotype
		# CHR	SNP_A	CHR	SNP_B	BETA	CHISQ	PVALUE	PHENOTYPE
		# 11	rs1508531	11	rs4939298	-0.13990	33.32995	7.77755E-09	ILMN_1716816
		line = line.strip() # strip() without argument will remove BOTH leading and trailing whitespace.
		fields = line.split() 
		chr_A, snp_A, chr_B, snp_B, pvalue, illumina_probe_id = fields[0], fields[1], fields[2], fields[3], fields[6], fields[7]


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


		################## *EXCLUDING SNPs* ##################
		# Pascal detected 04/22/2015 that the .bim file [Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.maf5.bim] contained DUPLICATED SNPs.
		# This gives problems downstream analysis (e.g. R snpStats). Rownames etc.
		# Here we exclude a list of "black_listed_rsID"
		### ***Write something HERE! ***
		# 2) check if snp_A or snp_B is in the list
			# 2a) if yes, exclude it/them
		flag_blacklisted = False
		if snp_A in blacklisted_snps_set: 
			counter_snp_blacklisted[snp_A] += 1
			flag_blacklisted = True
		elif snp_B in blacklisted_snps_set:
			counter_snp_blacklisted[snp_B] += 1
			flag_blacklisted = True

		if flag_blacklisted: # we continue loop if the flag is raised.
			continue

		################## SNP map and sets ##################
		### NOTE
		# We use "try/except" to catch exception occurring if a SNP is not in the SNP2interaction_dict.
		# This could happen if we for some reason have QC'ed the genotypes *AFTER* running FastEpistasis (.bim --> set files + XXX/stats/df_interaction_table_snps.txt --> SNP2interaction_dict).
		# We avoid errors by setting the set_A/B to an empty set. In this way, we EFFECTIVELY IGNORE SNPs present in the epistasis concatenated epistasis results.
		try:
			set_A = SNP2interaction_dict[snp_A] # TODO: put this into a try/except or use dict.get() for a default value
		except KeyError as e:
			set_A = set() # empty set | this is nicer to do, compared to a "continue"
			counter_snp_not_in_SNP2interaction[snp_A] += 1 # incrementing collections.Counter()

		try:
			set_B = SNP2interaction_dict[snp_B] # TODO: put this into a try/except or use dict.get() for a default value
		except KeyError as e:
			set_B = set() # empty set | this is nicer to do, compared to a "continue"
			counter_snp_not_in_SNP2interaction[snp_B] += 1 # incrementing collections.Counter()
		

		set_AB = set_A & set_B # intersection
		### REMEMBER: the set_AB contains strings of the following format:
			# experiment_interaction_identifier = {experiment_type}_{experiment_no}_{interaction_no}
		# e.g. 			
			# null_256_22233
			# null_1000_26312
			# ....


		### count_all
		#probe_counts_dict['count_all'][illumina_probe_id] += 1
		df_probe_counts.ix[illumina_probe_id, "count_all"] += 1


		### Mapping SNP-pairs to *UNASSIGNED* group
		if len(set_AB)==0:
			# assign to unassigned group
			count_unassigned += 1
			fh_unassigned.write(line + "\n")
		else: # assigned state
			### count_assigned 
			#probe_counts_dict['count_assigned'][illumina_probe_id] += 1
			df_probe_counts.ix[illumina_probe_id, "count_assigned"] += 1

			if pvalue <= bonferroni_correction_dict['hic_1']: # *OBS: hardcoded *experiment_identifier*!!!
				# ^^ It makes sense to use the cut-off for hic_1 experiment as significant
				### count_assigned_significant
				#probe_counts_dict['count_assigned_significant'][illumina_probe_id] += 1
				df_probe_counts.ix[illumina_probe_id, "count_assigned_significant"] += 1


		##############################################################################
		### Checking for *ASSIGNED* INTRA chromosomal interactions
		if (chr_A == chr_B) and (len(set_AB) > 0): # len(set_AB) > 0: only considering assigned SNPs:
			count_intrachromosomal_interactions += 1
			with open(file_epistatic_intrachromosomal, 'a') as fh: # OBS: append mode!
				fh.write( "{}\t{}\t{}\n".format(line, len(set_AB), ";".join(set_AB)) )
			continue # *IMPORTANT*: we skip *ASSINGED* *INTRACHROMOSOMAL* interactions
		##############################################################################

		### Mapping SNP-pairs to their experiments
		if len(set_AB)==1:
			### assign to one group
			count_assigned += 1

			experiment_interaction_identifier = iter(set_AB).next() # GET first (AND ONLY) element in set. e.g. null_256_22233
			#experiment_interaction_identifier = list(set_AB)[0]
			tmp_parts_list = experiment_interaction_identifier.split('_')
			experiment_identifier = tmp_parts_list[0] + "_" + tmp_parts_list[1] # e.g. null_256

			if pvalue <= bonferroni_correction_dict[experiment_identifier]:
				flag_significant = True
				### Count dicts
				epistatic_counts_dict['count_significant'][experiment_identifier] += 1
				epistatic_counts_dict['count_significant_pruned_EIID'][experiment_identifier].add(experiment_interaction_identifier)
			else:
				flag_significant = False

			epistatic_counts_dict['count_all'][experiment_identifier] += 1

			### Line output
			line_out = "{}\t{}\t{}".format(line, experiment_interaction_identifier, flag_significant)
			experiment_identifier_fh_dict[experiment_identifier].write( line_out + "\n" )
			
			### Populating epistasis_table_dict* ##
			if flag_significant:
				### Epistasis table dicts - ONLY CONSIDERING SIGNIFICANT epistatic SNP-pairs ###
				# epistasis_table_dict
				epistasis_table_dict[experiment_identifier].append( line_out )
				# epistasis_table_pruned_EIID_dict
				if pvalue < epistasis_table_pruned_EIID_dict[experiment_identifier][experiment_interaction_identifier]['best_pvalue']:
					# D[EID][EIID]['best_pvalue'] = pvalue
					epistasis_table_pruned_EIID_dict[experiment_identifier][experiment_interaction_identifier]['best_pvalue'] = pvalue
					epistasis_table_pruned_EIID_dict[experiment_identifier][experiment_interaction_identifier]['best_line'] = line_out
				# epistasis_table_pruned_hemani_dict
				if pvalue < epistasis_table_pruned_hemani_dict[experiment_identifier][illumina_probe_id][chr_pair]['best_pvalue']:
					# D[EID][illumina_probe_id][chr_pair]['best_pvalue'] = pvalue
					epistasis_table_pruned_hemani_dict[experiment_identifier][illumina_probe_id][chr_pair]['best_pvalue'] = pvalue
					epistasis_table_pruned_hemani_dict[experiment_identifier][illumina_probe_id][chr_pair]['best_line'] = line_out

			
		elif len(set_AB)>1:
			# "warn"
			# assigned to multiple groups
			# VALIDATION: check that hic_1 is *NOT* co-occurring with ANY null experiments
			count_assigned += 1 # OBS: SNP pairs with multiple assignments will only be counted once!
			#print "Found SNP pair with multiple assignments (#{})".format(len(set_AB))
			flag_null_seen = False
			for experiment_interaction_identifier in set_AB:
				tmp_parts_list = experiment_interaction_identifier.split('_')
				experiment_identifier = tmp_parts_list[0] + "_" + tmp_parts_list[1] # e.g. null_256

				if experiment_identifier.startswith("null"): # hmmm, there is a slight performance loss when doing this testing this 'if statement' many time. However, this will not affect the speed of the code seriously.
					flag_null_seen = True

				if (experiment_identifier == "hic_1") and (flag_null_seen): # checking for co-occurring of hic with ANY null experiments
					print "List of experiment_interaction_identifiers in set_AB:\n{}".format( "\n".join(set_AB) )
					print "DAMN, got hic_1 for a SNP-pair with 'multiple_assigments' that co-occurs with a null experiment. This means that the *NULL CONTAINS AT TRUE POSITIVE*."
					count_null_false_positives += 1
					with open(file_null_false_negatives, 'a') as fh: # OBS: append mode!
						fh.write( "{}\t{}\t{}\n".format(line, len(set_AB), ";".join(set_AB)) )

					# The reason for this 'false negatives' in the null, is that some interactions may be overlaping because of the 'width' of the regions.
					# Remember we use the same 'interationA' for the HiC and Null. Thus if a Null interaction for 'interactionB' falls suficiently close to the 'interactionB' from the HiC, they will have an overlap and the SNPs in this overlap will be in both a Null and HiC experiment.
					# Using a large width or many interaction will cause this 'false negative' problem to be bigger.
					# Remember that the R-script that contructs the null tables *ADVANCED* in the sense that it uses *chr:pos information* to ensure that *NO NULL [in interactionB] OVERLAP WITH THE HIC*.
					# *CONCLUSION*: this is not really a big problem if it happens rarely. We just need to know how often this occur to monitor the 'false negative' rate in the Null

					### IMPORTANT - recent Questions:
					# How is it possible to have two HiC experiments listed in the set_AB for two SNPs? (Should the code check be changed?)
					# How come two SNPs are listed as interchromosomal interactions?
						# --> the null must *NOT HAVE ANY INTRA-CHROMOSOMAL INTERACTIONS*.


					#raise Exception("See above print messages...")

				if pvalue <= bonferroni_correction_dict[experiment_identifier]:
					flag_significant = True
					epistatic_counts_dict['count_significant'][experiment_identifier] += 1
					epistatic_counts_dict['count_significant_pruned_EIID'][experiment_identifier].add(experiment_interaction_identifier)
				else:
					flag_significant = False

				epistatic_counts_dict['count_all'][experiment_identifier] += 1

				### Line output ###
				line_out = "{}\t{}\t{}".format(line, experiment_interaction_identifier, flag_significant)
				experiment_identifier_fh_dict[experiment_identifier].write( line_out + "\n" )
				
				### Populating epistasis_table_dict* ##
				if flag_significant:
					### Epistasis table dicts - ONLY CONSIDERING SIGNIFICANT epistatic SNP-pairs ###
					# epistasis_table_dict
					epistasis_table_dict[experiment_identifier].append( line_out )
					# epistasis_table_pruned_EIID_dict
					if pvalue < epistasis_table_pruned_EIID_dict[experiment_identifier][experiment_interaction_identifier]['best_pvalue']:
						# D[EID][EIID]['best_pvalue'] = pvalue
						epistasis_table_pruned_EIID_dict[experiment_identifier][experiment_interaction_identifier]['best_pvalue'] = pvalue
						epistasis_table_pruned_EIID_dict[experiment_identifier][experiment_interaction_identifier]['best_line'] = line_out
					# epistasis_table_pruned_hemani_dict
					if pvalue < epistasis_table_pruned_hemani_dict[experiment_identifier][illumina_probe_id][chr_pair]['best_pvalue']:
						# D[EID][illumina_probe_id][chr_pair]['best_pvalue'] = pvalue
						epistasis_table_pruned_hemani_dict[experiment_identifier][illumina_probe_id][chr_pair]['best_pvalue'] = pvalue
						epistasis_table_pruned_hemani_dict[experiment_identifier][illumina_probe_id][chr_pair]['best_line'] = line_out


			### increment count
			count_multiple_assignments += 1

			### multiple assignment file
			fh_multiple_assignments.write(line + "\t" + ";".join(set_AB) + "\n")


print "Done with main loop!"

###################################### PROCESS and WRITE epistasis_table_* content ######################################

################## POPULATE *epistatic_counts_dict* AND WRITE epistasis_table_* ##################
# Populate dict with data from epistasis_table_*


### epistasis_table_dict ###
## epistasis_table_dict[EID] = LIST of "line_out"
with open(file_epistasis_table, 'w') as fh:
	for EID in sorted(epistasis_table_dict, key=function_sort_EID):
		list_of_line_out = epistasis_table_dict[EID]
		for line_out in list_of_line_out:
			fh.write( line_out + "\n" )

		

### epistasis_table_pruned_EIID_dict ###
## epistasis_table_pruned_EIID_dict[EID][EIID]['best_pvalue'] = pvalue
with open(file_epistasis_table_pruned_EIID, 'w') as fh:
	for EID in sorted(epistasis_table_pruned_EIID_dict, key=function_sort_EID):
		EID_dict = epistasis_table_pruned_EIID_dict[EID]
			# EID_dict.keys() is EIID
			# EXAMPLE --> ['hic_1_2013', 'hic_1_2574']

		### The below can be used for validation/extra check ###
		##n_epistatic_counts = len(EID_dict[EID].keys()) # this is the number of EIID (experiment_interaction_identifier) that is significant for a given EID
		##epistatic_counts_dict["count_significant_pruned_EIID__alternative_counting_method"][EID] = n_epistatic_counts # no need for incrementing ("+=")

		for EIID in sorted(EID_dict, key=function_sort_EIID):
			EIID_dict = EID_dict[EIID]
				# EIID_dict.keys()
				# EXAMPLE [*FIXED*] --> ['best_pvalue', 'best_line']
			line_out = EIID_dict["best_line"]
			fh.write( line_out + "\n" )


### epistasis_table_pruned_hemani_dict ###
## epistasis_table_pruned_hemani_dict[EID][illumina_probe_id][chr_pair]['best_pvalue'] = pvalue
#for EID, EID_dict in epistasis_table_pruned_hemani_dict.items():
with open(file_epistasis_table_pruned_hemani, 'w') as fh:
	for EID in sorted(epistasis_table_pruned_hemani_dict, key=function_sort_EID):
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

################## epistatic_counts.txt ##################
# *OBS*: this file is similar to file_epistatic_counts_csv, except for the following:
# 1) file_epistatic_counts is not sorted
# 2) file_epistatic_counts contains a column about bonferroni correction
# 2) file_epistatic_counts has no header.
# 2) file_epistatic_counts is tab seperated.
with open(file_epistatic_counts, 'w') as fh:
	#for experiment_identifier in sorted(epistatic_counts_dict, key=epistatic_counts_dict.get, reverse=True):
	for experiment_identifier in bonferroni_correction_dict:
		fh.write( "{}\t{}\t{}\t{}\t{}\n".format(experiment_identifier, 
												epistatic_counts_dict['count_all'][experiment_identifier], 
												epistatic_counts_dict['count_significant'][experiment_identifier], 
												len(epistatic_counts_dict['count_significant_pruned_EIID'][experiment_identifier]), 
												bonferroni_correction_dict[experiment_identifier]),
												) # epistatic_counts_dict + epistatic_counts_dict['count_significant']
		#fh.write( "{}\t{}\n".format(experiment_identifier, epistatic_counts_dict['count_all'][experiment_identifier]) ) # 


###################################### CREATING DataFrame: Calculating Empirical P-value ######################################

#df_experiment_identifier_counts = pd.DataFrame()
df_experiment_identifier_counts = pd.DataFrame(columns=epistatic_counts_dict.keys(), index=bonferroni_correction_dict.keys()) # initialize data frame with NaNs.
df_experiment_identifier_counts = df_experiment_identifier_counts.fillna(0) # fill the data with 0s rather than NaNs
### COLUMNS IN DATAFRAME: "keys of epistatic_counts_dict"
	# count_all
	# count_significant
	# count_significant_pruned_EIID


for key_master in epistatic_counts_dict:
	for experiment_identifier in epistatic_counts_dict[key_master]:
	#for experiment_identifier in bonferroni_correction_dict:
		#df_experiment_identifier_counts[key_master][]
		if key_master == 'count_significant_pruned_EIID':
			n_epistatic_counts = len(epistatic_counts_dict[key_master][experiment_identifier]) # length of set()
			### Validation: the two methods of counting should give the same
			assert( n_epistatic_counts == len(epistasis_table_pruned_EIID_dict[experiment_identifier].keys()) )
			#if not n_epistatic_counts == len(epistasis_table_pruned_EIID_dict[experiment_identifier].keys()):
				#pdb.set_trace()
			
			### assignment to data frame
			df_experiment_identifier_counts.ix[experiment_identifier, key_master] = n_epistatic_counts
		else:
			df_experiment_identifier_counts.ix[experiment_identifier, key_master] = epistatic_counts_dict[key_master][experiment_identifier]
		#print key_master, experiment_identifier
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
tmp_sort_order_list = ["count_assigned_significant", "count_all", "count_assigned"]
df_probe_counts.sort(tmp_sort_order_list, ascending=False, inplace=True)
### Writing file
df_probe_counts.to_csv(file_probe_counts_csv) # sep='\t', index=True, header=True


###################################### Write stats file ######################################
################## Create ... ##################

################## Write ##################
with open(file_epistatic_stats, 'w') as fh:
	### Counts
	fh.write( "count_unassigned: {}\n".format(count_unassigned) )
	fh.write( "count_assigned: {} ({:.2f} % of total)\n".format(count_assigned, count_assigned/float(count_total)*100) )
	fh.write( "count_multiple_assignments: {} ({:.2f} % of count_assigned)\n".format(count_multiple_assignments, count_multiple_assignments/float(count_assigned)*100) )
	fh.write( "count_null_false_positives: {}\n".format(count_null_false_positives) )
	fh.write( "count_intrachromosomal_interactions: {} ({:.2f} % of total)\n".format(count_intrachromosomal_interactions, count_intrachromosomal_interactions/float(count_total)*100) )

	### Counters
	fh.write( "counter_snp_blacklisted: [{}]\n".format(counter_snp_blacklisted) )
	fh.write( "sum(counter_snp_blacklisted.values()): {}\n".format(sum(counter_snp_blacklisted.values())) )

	fh.write( "counter_snp_not_in_SNP2interaction: [{}]\n".format(counter_snp_not_in_SNP2interaction) )
	fh.write( "sum(counter_snp_not_in_SNP2interaction.values()): {}\n".format(sum(counter_snp_not_in_SNP2interaction.values())) )


	### Probes
	probe_stats_count_all = sum(df_probe_counts["count_all"] > 0) # --> could also use "(df_probe_counts["count_all"] > 0).sum()". Gives the same. TESTED.
	probe_stats_count_assigned = sum(df_probe_counts["count_assigned"] > 0)
	probe_stats_count_assigned_significant = sum(df_probe_counts["count_assigned_significant"] > 0)
	fh.write( "PROBE count_all: {} ({:.2f} % of total number of probes in probes_dict)\n".format( probe_stats_count_all, probe_stats_count_all/float(n_probes)*100 ) )	
	fh.write( "PROBE count_assigned: {} ({:.2f} % of total number of probes in probes_dict)\n".format( probe_stats_count_assigned, probe_stats_count_assigned/float(n_probes)*100 ) )	
	fh.write( "PROBE count_assigned_significant: {} ({:.2f} % of total number of probes in probes_dict)\n".format( probe_stats_count_assigned_significant, probe_stats_count_assigned_significant/float(n_probes)*100 ) )	
		

	### Hi-C
	fh.write( "HIC count_all: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_all']) )
	fh.write( "HIC count_significant: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant']) )
	fh.write( "HIC count_significant_pruned_EIID: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant_pruned_EIID']) )
	fh.write( "HIC count_significant_pruned_hemani: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant_pruned_hemani']) )

	fh.write( "p_value_dict[count_all]: {}\n".format(p_value_dict['count_all']) )
	fh.write( "p_value_dict[count_significant]: {}\n".format(p_value_dict['count_significant']) )
	fh.write( "p_value_dict[count_significant_pruned_EIID]: {}\n".format(p_value_dict['count_significant_pruned_EIID']) )
	fh.write( "p_value_dict[count_significant_pruned_hemani]: {}\n".format(p_value_dict['count_significant_pruned_hemani']) )
	### OLD
	# fh.write( "p_value_count_significant_pruned: {}\n".format(p_value_count_significant_pruned) )
	# fh.write( "p_value_count_significant: {}\n".format(p_value_count_significant) )
	# fh.write( "p_value_count_all: {}\n".format(p_value_count_all) )

###################################### Closing file handle ######################################

print "Closing filehandles..."
fh_multiple_assignments.close()
fh_unassigned.close()

#for experiment_interaction_identifier in experiment_identifier_fh_dict:
for fh in experiment_identifier_fh_dict.values():
	fh.close()


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

###################################### GARBAGE ######################################

	# *Main* path where the FastEpistasis results lives.
	# Expected folder organization in path_main_input:
	# <path_main_input>/
	# 	errors/ 
	# 	LSF_jobs_ids/								[optional]
	# 	cat_epistasis_YYMMDD_HHMMSS/				[optional]
	# 	config.json									[optional]
	# 	epi/										[optional]
	# 	fastEpi_compiled/							[*required*]
	# 		XXX___probe_epistatic_count.txt
	# 		XXX___results.epi.qt.lm.combined.txt	[*required*]
	# 		XXX___stat_file.txt
	# 		SNP2interaction_map/
	# 			SNP2interaction_map.pickle 			[*required*]
	# 			MORE_FILES...
	# 	jobs/										[optional]
	# 	link_bim									[optional]
	# 	link_probes									[optional]
	# 	link_set									[optional]
	# 	logs/										[optional]
	# 	logs_worker/								[optional]

	# The program relies on the *specific names* in fastEpi_compiled to match the specified pattern.

###################################### ABOUT MATPLOTLIB Axes ######################################
### Ref: http://stackoverflow.com/a/11786526

# There is a helpful page here, which provides an overview of the classes in matplotlib.

# Essentially, the process is:

# Create a figure which can hold Axes instances (and other artists)
# Create a canvas for the figure to draw itself to
# Create an Axes instance, ax, to which plotted lines/patches etc can be added. e.g. ax.plot(range(10)) or ax.contourf(array)
# I think your confusion comes from the understanding of what an Axes is. It is "the rectangular area which holds the basic elements" (for rectilinear plots). By default there is only one Axes in a figure, no matter how many times you run the command plt.plot(range(10)), although you may decide to use plt.subplot to have sub-plots in your figure, in which case you would have many Axes instances in your figure.

# HTH,


###################################### ?? ######################################

# def example_plot(ax, fontsize=12):
#      ax.plot([1, 2])
#      ax.locator_params(nbins=3)
#      ax.set_xlabel('x-label', fontsize=fontsize)
#      ax.set_ylabel('y-label', fontsize=fontsize)
#      ax.set_title('Title', fontsize=fontsize)

# plt.close('all')
# fig, ax = plt.subplots()
# example_plot(ax, fontsize=24)

