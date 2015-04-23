#!/usr/bin/env python

import os
import sys
import time

import argparse

import re
import collections

try:
	import cPickle as pickle # NO subclass - cPickle is written in C.
except:
	print "loaded pickle - could not load cPickle..."
	import pickle

import glob

import resource # to allow enough filehandles to be opened

###################################### Library ######################################
import shutil # for copying probe_list

import epistasis_table_library


###################################### USAGE ######################################

#python epistasis_table_labeller.py --path_main_input XXX


################## Broad ##################
##### hIMR90 #####

##### hESC #####

################## OSX ##################

##### hIMR90 #####
# python epistasis_table_labeller.py --path_main_input /Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8_test_case

##### hESC #####

###################################### SYNOPSIS ######################################

# ### Output files ### 
# 	<path_main_input>/ [e.g. "fastEpi_compiled/"]
# 		<path_main_out>/ [i.e. "assigned/"]
# 			<path_out_epistatic_results>/ ["epistatic_results/"]
# 				<EID>.lm.combined.txt
# 				multiple_assigments.lm.combined.txt
# 				unassigned.lm.combined.txt
# 			file_epistasis_table [i.e. "epistasis_table.txt"] --> *IMPORTANT*
#			probe_list.txt	--> *IMPORTANT*
# 			experiments_list.txt --> *IMPORTANT*
# 			null_false_negatives.txt
# 			epistatic_intrachromosomal.txt

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
			XXX___probe_list.txt					[*required*] [NEW]
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

### Filenames - epistatic_results multiple assignments + more
file_multiple_assignments = "{path}/{basename}.{ext}".format(path=path_out_epistatic_results, basename="multiple_assigments", ext=extension_epistatic_results)
file_unassigned = "{path}/{basename}.{ext}".format(path=path_out_epistatic_results, basename="unassigned", ext=extension_epistatic_results)

### Epistasis: epistasis tables
file_epistasis_table = path_main_out + "/epistasis_table.txt"

### Epistasis: experiment counts and stats file
file_experiments_list = path_main_out + "/experiments_list.txt"
file_epistatic_stats = path_main_out + "/stats_labeller.txt"

file_null_false_negatives = path_main_out + "/null_false_negatives.txt" # 
file_epistatic_intrachromosomal = path_main_out + "/epistatic_intrachromosomal.txt" # 

### Probes
file_probe_list_out = path_main_out + "/probe_list.txt" # 

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
assert(len(glob_probe_list_file)==1) # Check that file exists AND that ONLY ONE file exists
file_probe_list_in = glob_probe_list_file[0] # because of the above check, we know that the list only contains one element


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


###################################### Read/Load input filesfile_SNP2interaction_map ######################################
filesize_SNP2interaction_map = os.path.getsize(file_SNP2interaction_map) >> 20 # BitWise Operations: shift-right | x >> y :Returns x with the bits shifted to the right by y places.
print "Loading SNP2interaction_dict via pickled file (size={} MB) ...".format(filesize_SNP2interaction_map)
time_start_tmp = time.time()
with open(file_SNP2interaction_map, 'r') as fh: # perhaps read the pickle file in binary mode.
	SNP2interaction_dict = pickle.load(fh)
print "Done loading pickled file in {:.2f} min".format( (time.time()-time_start_tmp)/60 )

###################################### Get *BLACKLISTED SNPs* ######################################

blacklisted_snps_set = epistasis_table_library.get_blacklisted_snps() # returns set()
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
		bonferroni_correction_dict[experiment_identifier] = epistasis_table_library.ALPHA/(int(fields[1])*epistasis_table_library.N_PHENOTYPES_TESTED) # <-- FLOAT. alpha/(N_epistasis_tests*N_phenotypes_tests)

		### Set filename
		tmp_filename = "{path}/{basename}.{ext}".format(path=path_out_epistatic_results, basename=experiment_identifier, ext=extension_epistatic_results)
		tmp_fh = open(tmp_filename, 'w')

		experiment_identifier_fh_dict[experiment_identifier] = tmp_fh # e.g. experiment_identifier_fh_dict['null_1'] = <FILE_HANDLE>

print "Opened {} filehandles".format(len(experiment_identifier_fh_dict))

################## Open *ADDITONAL* filehandles ##################

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


################## Dict for epistasis_table_*.txt ##################
# EID = Experimental IDentifier. E.g. hic_1, null_1
# EIID = Experimental Interaction IDentifier. E.g. hic_1_1923, null_1_119
epistasis_table_dict = collections.defaultdict(list) # *OBS*: list
# 	# D[EID].append(line)


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

		### Mapping SNP-pairs to *UNASSIGNED* group
		if len(set_AB)==0:
			# assign to unassigned group
			count_unassigned += 1
			fh_unassigned.write(line + "\n")

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
			else:
				flag_significant = False

			### Line output
			line_out = "{}\t{}\t{}".format(line, experiment_interaction_identifier, experiment_identifier)
			experiment_identifier_fh_dict[experiment_identifier].write( line_out + "\n" )
			
			### Populating epistasis_table_dict* ##
			if flag_significant:
				### Epistasis table dicts - ONLY CONSIDERING SIGNIFICANT epistatic SNP-pairs ###
				# epistasis_table_dict
				epistasis_table_dict[experiment_identifier].append( line_out )
				
			
		elif len(set_AB)>1:
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
				else:
					flag_significant = False

				### Line output ###
				line_out = "{}\t{}\t{}".format(line, experiment_interaction_identifier, experiment_identifier)
				experiment_identifier_fh_dict[experiment_identifier].write( line_out + "\n" )
				
				### Populating epistasis_table_dict* ##
				if flag_significant:
					### Epistasis table dicts - ONLY CONSIDERING SIGNIFICANT epistatic SNP-pairs ###
					# epistasis_table_dict
					epistasis_table_dict[experiment_identifier].append( line_out )

			### increment count
			count_multiple_assignments += 1

			### multiple assignment file
			fh_multiple_assignments.write(line + "\t" + ";".join(set_AB) + "\n")


print "Done with main loop!"

###################################### Copy probe list ######################################
### Copy probe_list to the path_main_out
print "Copying probe_list to path_main_out destination"
shutil.copyfile(file_probe_list_in, file_probe_list_out) # shutil.copyfile(src, dst)

###################################### PROCESS and WRITE epistasis_table_* content ######################################

################## WRITE epistasis_table_dict ##################

### DATA Structure --> epistasis_table_dict[EID] = LIST of "line_out"
header_epistasis_table = [
	"CHR_1",
	"SNP_1",
	"CHR_2",
	"SNP_2",
	"BETA",
	"CHISQ",
	"PVALUE",
	"PHENOTYPE",
	"EIID",
	"EID"
	] # --> NCOLS = 10
### EXAMPLE OUTPUT
# CHR_1	SNP_1	CHR_2	SNP_2	BETA	CHISQ	PVALUE	PHENOTYPE	EIID	EID
# 1	rs12564773	3	rs1514423	+0.10238	41.25541	1.33581E-10	ILMN_1653904	hic_1_1197	hic_1
# 1	rs12564755	3	rs1514423	+0.10238	41.25541	1.33581E-10	ILMN_1653904	hic_1_1197	hic_1
# 3	rs17009271	4	rs17676993	-0.41600	57.89248	2.76844E-14	ILMN_1659762	hic_1_1512	hic_1

with open(file_epistasis_table, 'w') as fh:
	fh.write( "\t".join(header_epistasis_table) + "\n" )
	for EID in sorted(epistasis_table_dict, key=epistasis_table_library.function_sort_EID):
		list_of_line_out = epistasis_table_dict[EID]
		for line_out in list_of_line_out:
			fh.write( line_out + "\n" )


###################################### Writing MORE files ######################################

################## experiments_list.txt ##################
with open(file_experiments_list, 'w') as fh:
	for EID in sorted(bonferroni_correction_dict, key=epistasis_table_library.function_sort_EID):
		fh.write( "{}\t{}\n".format(EID, bonferroni_correction_dict[EID]) )

###################################### Write stats file ######################################

################## Write ##################
with open(file_epistatic_stats, 'w') as fh:
	fh.write("N_PHENOTYPES_TESTED = {}\n".format(epistasis_table_library.N_PHENOTYPES_TESTED) )
	fh.write("ALPHA = {}\n".format(epistasis_table_library.ALPHA) )
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


	# ### Probes
	# probe_stats_count_all = sum(df_probe_counts["count_all"] > 0) # --> could also use "(df_probe_counts["count_all"] > 0).sum()". Gives the same. TESTED.
	# probe_stats_count_assigned = sum(df_probe_counts["count_assigned"] > 0)
	# probe_stats_count_assigned_significant = sum(df_probe_counts["count_assigned_significant"] > 0)
	# fh.write( "PROBE count_all: {} ({:.2f} % of total number of probes in probes_dict)\n".format( probe_stats_count_all, probe_stats_count_all/float(n_probes)*100 ) )	
	# fh.write( "PROBE count_assigned: {} ({:.2f} % of total number of probes in probes_dict)\n".format( probe_stats_count_assigned, probe_stats_count_assigned/float(n_probes)*100 ) )	
	# fh.write( "PROBE count_assigned_significant: {} ({:.2f} % of total number of probes in probes_dict)\n".format( probe_stats_count_assigned_significant, probe_stats_count_assigned_significant/float(n_probes)*100 ) )	
		

	# ### Hi-C
	# fh.write( "HIC count_all: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_all']) )
	# fh.write( "HIC count_significant: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant']) )
	# fh.write( "HIC count_significant_pruned_EIID: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant_pruned_EIID']) )
	# fh.write( "HIC count_significant_pruned_hemani: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant_pruned_hemani']) )

	# fh.write( "p_value_dict[count_all]: {}\n".format(p_value_dict['count_all']) )
	# fh.write( "p_value_dict[count_significant]: {}\n".format(p_value_dict['count_significant']) )
	# fh.write( "p_value_dict[count_significant_pruned_EIID]: {}\n".format(p_value_dict['count_significant_pruned_EIID']) )
	# fh.write( "p_value_dict[count_significant_pruned_hemani]: {}\n".format(p_value_dict['count_significant_pruned_hemani']) )

###################################### Closing file handle ######################################

print "Closing filehandles..."
fh_multiple_assignments.close()
fh_unassigned.close()

#for experiment_interaction_identifier in experiment_identifier_fh_dict:
for fh in experiment_identifier_fh_dict.values():
	fh.close()



###################################### Finishing ######################################
print "The script is complete!"
