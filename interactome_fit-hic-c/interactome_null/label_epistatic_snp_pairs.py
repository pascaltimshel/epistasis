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

###################################### USAGE ######################################

#python label_epistatic_snp_pairs.py --path_main_input XXX


################## Broad ##################
### hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10
#python label_epistatic_snp_pairs.py --path_main_input /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10/fastEpi_compiled
	# ---> runtime: XXX

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
extension_epistatic_results = "lm.combined.txt" # file extension name

file_epistatic_counts = path_main_out + "/epistatic_counts.txt"
file_epistatic_counts_csv = path_main_out + "/epistatic_counts.csv"
file_epistatic_stats = path_main_out + "/epistatic_stats.txt"

file_null_false_negatives = path_main_out + "/null_false_negatives.txt" # 
file_epistatic_intrachromosomal = path_main_out + "/epistatic_intrachromosomal.txt" # 

#file_epistatic_counts_significant = path_main_out + "/epistatic_counts_significant.txt"

################## Get input files ##################
### Get FastEpistasis compiled result file - GLOBBING ON FILENAME PATTERN
glob_lm_file = glob.glob(path_main_input+"/*lm.combined.txt")
if not len(glob_lm_file) == 1:
	raise Exception( "glob_lm_file does not contain EXACTLY one matching file. Matches are:\n[{}]".format("\n".join(glob_lm_file)) )
file_fastepistasis_lm_combined = glob_lm_file[0] # because of the above check, we know that the list only contains one element

file_SNP2interaction_map = path_main_input + "/SNP2interaction_map/SNP2interaction_map.pickle"
file_bonferroni_correction = path_main_input + "/SNP2interaction_map/bonferroni_correction.txt"

###################################### INITIAL Error checks ######################################

################## Check that files exists ##################
if not os.path.exists(file_SNP2interaction_map):
	raise Exception("file_SNP2interaction_map does not exists: %s" % file_SNP2interaction_map)

###################################### OUTPUT ######################################
if os.path.exists(path_main_out):
	print "path_main_out={} exists".format(path_main_out)
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

###################################### Read/Load input filesfile_SNP2interaction_map ######################################
filesize_SNP2interaction_map = os.path.getsize(file_SNP2interaction_map) >> 20 # BitWise Operations: shift-right | x >> y :Returns x with the bits shifted to the right by y places.
print "Loading SNP2interaction_dict via pickled file (size={} MB) ...".format(filesize_SNP2interaction_map)
time_start_tmp = time.time()
with open(file_SNP2interaction_map, 'r') as fh: # perhaps read the pickle file in binary mode.
	SNP2interaction_dict = pickle.load(fh)
print "Done loading pickled file in {:.2f} min".format( (time.time()-time_start_tmp)/60 )

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

### READ and WRITE line-by-line
count_total = 0
count_unassigned = 0
count_assigned = 0
count_null_false_positives = 0
count_intrachromosomal_interactions = 0
#epistatic_counts_dict = collections.defaultdict(lambda: collections.defaultdict(int)) # e.g. epistatic_counts_dict['MASTER_KEY']['EXPERIMENT_IDENTIFIER'] = INTEGER
	# the expression inside the parenthesis must be a 'callable'. That is why I use a 'lambda function'
epistatic_counts_dict = dict() #collections.defaultdict(dict)
epistatic_counts_dict['count_significant_pruned'] = collections.defaultdict(set)
epistatic_counts_dict['count_significant'] = collections.defaultdict(int)
epistatic_counts_dict['count_all'] = collections.defaultdict(int)
	### MASTER KEYS IN DICT:
	# count_significant_pruned
	# count_significant
	# count_all

#epistatic_counts_dict['count_significant'] = collections.defaultdict(int)

### Count the number of lines in the file
n_lines_fastepistasis_lm_combined = sum(1 for line in open(file_fastepistasis_lm_combined)) # this is pretty fast (only a bit slower than "wc -l"). 3.5 sec for a 1 GB file

time_start_loop = time.time()

with open(file_fastepistasis_lm_combined, 'r') as fh_compiled:
	next(fh_compiled) # SKIPPING THE FIRST LINE!
	for line_no, line in enumerate(fh_compiled, start=1):
		count_total += 1

		time_elapsed_loop = time.time() - time_start_loop # <type 'float'>
		if line_no % 5000 == 0: 
			print "Main loop | #{line_no}/#{n_lines} | {pct_complete:.2f} % done | {sec:.2f} sec [{min:.2f} min]".format(line_no=line_no, n_lines=n_lines_fastepistasis_lm_combined, pct_complete=(line_no/float(n_lines_fastepistasis_lm_combined))*100, sec=time_elapsed_loop, min=time_elapsed_loop/60)
		


		#chr1, snp_A, chr2, snp_B, beta, chisq, pvalue, phenotype
		# CHR	SNP_A	CHR	SNP_B	BETA	CHISQ	PVALUE	PHENOTYPE
		# 11	rs1508531	11	rs4939298	-0.13990	33.32995	7.77755E-09	ILMN_1716816
		line = line.strip() # strip() without argument will remove BOTH leading and trailing whitespace.
		fields = line.split() 
		chr_A, snp_A, chr_B, snp_B, pvalue = fields[0], fields[1], fields[2], fields[3], float(fields[6])

		set_A = SNP2interaction_dict[snp_A] # TODO: put this into a try/except or use dict.get() for a default value
		set_B = SNP2interaction_dict[snp_B]

		set_AB = set_A & set_B # intersection
		### REMEMBER: the set_AB contains strings of the following format:
			# experiment_interaction_identifier = {experiment_type}_{experiment_no}_{interaction_no}
		# e.g. 			
			# null_256_22233
			# null_1000_26312
			# ....

		### Mapping SNP-pairs to *UNASSIGNED* group
		if len(set_AB)==0:
			# assign to unassigned group
			count_unassigned += 1
			fh_unassigned.write(line + "\n")

		##############################################################################
		### Checking for *ASSINGED* INTRA chromosomal interactions
		if (chr_A == chr_B) and (len(set_AB) > 0):
			count_intrachromosomal_interactions += 1
			with open(file_epistatic_intrachromosomal, 'a') as fh: # OBS: append mode!
				fh.write( "{}\t{}\t{}\n".format(line, len(set_AB), ";".join(set_AB)) )
			continue # *IMPORTANT*: we skip *ASSINGED* *INTRACHROMOSOMAL* interactions
		##############################################################################

		### Mapping SNP-pairs to their experiments
		if len(set_AB)==1:
			# assign to one group
			count_assigned += 1
			experiment_interaction_identifier = iter(set_AB).next() # e.g. null_256_22233
			tmp_parts_list = experiment_interaction_identifier.split('_')
			experiment_identifier = tmp_parts_list[0] + "_" + tmp_parts_list[1] # e.g. null_256

			if pvalue <= bonferroni_correction_dict[experiment_identifier]:
				flag_significant = True
				epistatic_counts_dict['count_significant'][experiment_identifier] += 1
				epistatic_counts_dict['count_significant_pruned'][experiment_identifier].add(experiment_interaction_identifier)
			else:
				flag_significant = False

			epistatic_counts_dict['count_all'][experiment_identifier] += 1

			#experiment_interaction_identifier = list(set_AB)[0]
			experiment_identifier_fh_dict[experiment_identifier].write( "{}\t{}\t{}\n".format(line, experiment_interaction_identifier, flag_significant) )
			
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
					epistatic_counts_dict['count_significant_pruned'][experiment_identifier].add(experiment_interaction_identifier)
				else:
					flag_significant = False

				epistatic_counts_dict['count_all'][experiment_identifier] += 1

				experiment_identifier_fh_dict[experiment_identifier].write( "{}\t{}\t{}\n".format(line, experiment_interaction_identifier, flag_significant) )
				
					
			### multiple assignment file
			fh_multiple_assignments.write(line + "\t" + ";".join(set_AB) + "\n")

print "Done with main loop!"


###################################### Writing files ######################################

with open(file_epistatic_counts, 'w') as fh:
	#for experiment_identifier in sorted(epistatic_counts_dict, key=epistatic_counts_dict.get, reverse=True):
	for experiment_identifier in bonferroni_correction_dict:
		fh.write( "{}\t{}\t{}\t{}\t{}\n".format(experiment_identifier, 
												epistatic_counts_dict['count_all'][experiment_identifier], 
												epistatic_counts_dict['count_significant'][experiment_identifier], 
												len(epistatic_counts_dict['count_significant_pruned'][experiment_identifier]), 
												bonferroni_correction_dict[experiment_identifier]),
												) # epistatic_counts_dict + epistatic_counts_dict['count_significant']
		#fh.write( "{}\t{}\n".format(experiment_identifier, epistatic_counts_dict['count_all'][experiment_identifier]) ) # 


###################################### CREATING DataFrame: Calculating Empirical P-value ######################################

#df_experiment_identifier_counts = pd.DataFrame()
df_experiment_identifier_counts = pd.DataFrame(columns=epistatic_counts_dict.keys(), index=bonferroni_correction_dict.keys()) # initialize data frame with zeroes.
df_experiment_identifier_counts = df_experiment_identifier_counts.fillna(0) # fill the data with 0s rather than NaNs
### COLUMNS IN DATAFRAME: "keys of epistatic_counts_dict"
	# count_all
	# count_significant
	# count_significant_pruned


for key_master in epistatic_counts_dict:
	for experiment_identifier in epistatic_counts_dict[key_master]:
	#for experiment_identifier in bonferroni_correction_dict:
		#df_experiment_identifier_counts[key_master][]
		if key_master == 'count_significant_pruned':
			df_experiment_identifier_counts.ix[experiment_identifier, key_master] = len(epistatic_counts_dict[key_master][experiment_identifier])
		else:
			df_experiment_identifier_counts.ix[experiment_identifier, key_master] = epistatic_counts_dict[key_master][experiment_identifier]
		#print key_master, experiment_identifier
		#if 'c' in df_experiment_identifier_counts.index: pdb.set_trace()

################## Creating P-value dict ##################
### About p-values: "obtaining a result EQUAL TO or MORE EXTREME than what was actually observed" --> p_val = sum(X >= X_OBS)
p_value_dict = {} # KEYS will be: count_significant_pruned, count_significant, count_all
for key_master in epistatic_counts_dict:
	p_value_dict[key_master] = sum(df_experiment_identifier_counts.ix[:, key_master] >= df_experiment_identifier_counts.ix['hic_1', key_master])/float(len(df_experiment_identifier_counts))
### OLD:
# p_value_count_significant_pruned = sum(df_experiment_identifier_counts.ix[:, 'count_significant_pruned'] >= df_experiment_identifier_counts.ix['hic_1', 'count_significant_pruned'])/float(len(df_experiment_identifier_counts))
# p_value_count_significant = sum(df_experiment_identifier_counts.ix[:, 'count_significant'] >= df_experiment_identifier_counts.ix['hic_1', 'count_significant'])/float(len(df_experiment_identifier_counts))
# p_value_count_all = sum(df_experiment_identifier_counts.ix[:, 'count_all'] >= df_experiment_identifier_counts.ix['hic_1', 'count_all'])/float(len(df_experiment_identifier_counts))

################## Write csv to file ##################
### Sorting, inplace
df_experiment_identifier_counts.sort(['count_significant_pruned', 'count_significant'], ascending=False, inplace=True)
### Writing file
df_experiment_identifier_counts.to_csv(file_epistatic_counts_csv) # sep='\t', index=True, header=True

################## Write stats file ##################

with open(file_epistatic_stats, 'w') as fh:
	fh.write( "count_unassigned: {}\n".format(count_unassigned) )
	fh.write( "count_assigned: {} ({:.2f} %)\n".format(count_assigned, count_assigned/float(count_unassigned+count_assigned)) )
	fh.write( "count_null_false_positives: {}\n".format(count_null_false_positives) )
	fh.write( "count_intrachromosomal_interactions: {} ({:.2f} %)\n".format(count_intrachromosomal_interactions, count_intrachromosomal_interactions/float(count_total)) )

	fh.write( "HIC count_significant_pruned: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant_pruned']) )
	fh.write( "HIC count_significant: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant']) )
	fh.write( "HIC count_all: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_all']) )

	fh.write( "p_value_dict[count_significant_pruned]: {}\n".format(p_value_dict['count_significant_pruned']) )
	fh.write( "p_value_dict[count_significant]: {}\n".format(p_value_dict['count_significant']) )
	fh.write( "p_value_dict[count_all]: {}\n".format(p_value_dict['count_all']) )

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


################## Subplot histogram ##################
subplot_n_row = len(df_experiment_identifier_counts.columns)
print "subplot_n_row={}".format(subplot_n_row)
for subplot_no, (column_name, series) in enumerate(df_experiment_identifier_counts.iteritems(), start=0): # Iterator over (column, series) pairs
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

fig_filename = "{base}/fig_{plot}.{ext}".format(base=path_main_out, plot="subplot", ext="pdf")
### saving fig
plt.savefig(fig_filename)
plt.close() # or plt.close(fig)


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

