#!/usr/bin/env python2.7

import os
import subprocess
import glob
import datetime

import time
import argparse

import pdb

###################################### USAGE ######################################
#python cat_fast_epistasis_subprocess.py --path_main_fastepistasis /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10

###################################### PARAMETERs ###################################### 
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--test_n_probes", type=int, help="INTERGER. Number of probes to test. [Parameter not required]") # defaults to None
arg_parser.add_argument("--path_main_fastepistasis", required=True, help="fastEpi main dir. MAY BE RELATIVE OR ABSOLUTE. E.g. 'fastEpi_1e-10_0108_212539' NOT 'fastEpi_1e-10_0108_212539/fastEpi_out'.")
args = arg_parser.parse_args()

path_main_fastepistasis = os.path.abspath(args.path_main_fastepistasis) # IMPORTANT: converts to ABSOLUTE PATH, e.g. 'fastEpi_1e-10_0108_212539' --> '/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpi_1e-10_0108_212539'
test_n_probes = args.test_n_probes
if test_n_probes is not None:
	print "running in test mode: test_n_probes=%s" % test_n_probes
else:
	print "running in full mode (not test mode with a subset of probes)"


###################################### Input checks ######################################
if not os.path.exists(path_main_fastepistasis):
	raise Exception("path_main_fastepistasis does not exists. path_main_fastepistasis=%s" % path_main_fastepistasis)


###################################### Preparing directories ######################################
#time_stamp = datetime.datetime.now().strftime("%a_%b_%d_%Y_%H%M%S")
time_stamp = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
#work_dir = "fastEpiCalc_epi1_%s_%s" % (epi1, time_stamp)

pathname_cat_epistasis = "cat_epistasis_%s" % time_stamp
path_cat_epistasis = os.path.join(path_main_fastepistasis, pathname_cat_epistasis) 
	# e.g. --> XXX/fastEpistasis_fit-hi-ci/hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10/cat_epistasis_<TIMESTAMP>

### Create directory
os.mkdir(path_cat_epistasis) # we KNOW that this dir does not exists BECAUSE WE USE A TIMESTAMP

#### IMPORTANT - SETTING RESULT PATH #### - used in TWO PLACES
# 1) reading finished probes (if switch enabled)
# 2) locating post_epi files
path_epi_results = path_main_fastepistasis + "/epi"

###################################### Filenames ######################################

### Output
# file_lm_combined = "{path}/results.epi.qt.lm.combined___{time_stamp}___{parameters}.txt".format( path=path_cat_epistasis, time_stamp=time_stamp, parameters=os.path.basename(path_main_fastepistasis) )
# file_no_post_out_file = "{path}/log_no_postFastEpistasis_file___{time_stamp}___{parameters}.txt".format( path=path_cat_epistasis, time_stamp=time_stamp, parameters=os.path.basename(path_main_fastepistasis) )
# file_stat = "{path}/stat_file___{time_stamp}___{parameters}.txt".format( path=path_cat_epistasis, time_stamp=time_stamp, parameters=os.path.basename(path_main_fastepistasis) )

file_lm_combined = "{path}/{parameters}__{time_stamp}___results.epi.qt.lm.combined.txt".format( path=path_cat_epistasis, parameters=os.path.basename(path_main_fastepistasis), time_stamp=time_stamp )
file_no_post_out_file = "{path}/{parameters}__{time_stamp}___log_no_postFastEpistasis_file.txt".format( path=path_cat_epistasis, parameters=os.path.basename(path_main_fastepistasis), time_stamp=time_stamp )
file_stat = "{path}/{parameters}__{time_stamp}___stat_file.txt".format( path=path_cat_epistasis, parameters=os.path.basename(path_main_fastepistasis), time_stamp=time_stamp )
file_probe_epistatic_count = "{path}/{parameters}__{time_stamp}___probe_epistatic_count.txt".format( path=path_cat_epistasis, parameters=os.path.basename(path_main_fastepistasis), time_stamp=time_stamp )
file_probe_list = "{path}/{parameters}__{time_stamp}___probe_list.txt".format( path=path_cat_epistasis, parameters=os.path.basename(path_main_fastepistasis), time_stamp=time_stamp )


###################################### LOADING PROBES ######################################
flag_probe_loading_method = "NOT_SET_YET" # default value

################## OPTION #1: read all probes in original dir ##################
if True: # TURNED ON
	flag_probe_loading_method = "ORIGINAL_DIR"
	print "Reading list of all probes in original dir"
	## READ SYMLINK FROM "path_main_fastepistasis" dir (result dir from fastEpi run)
	probe_dir_symlink = path_main_fastepistasis + "/link_probes"
	### ** READING all FILENAMEs in probe_dir ** ###
	probes = glob.glob(probe_dir_symlink+"/*") # *OBS* this is DIRTY. When we change dir later this will not work. Thus the probes CANNOT be changed (which is actually ok.)
	# e.g. /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpiCalc_epi1_1e-6_Thu_Jan_08_2015_194154/symlink_probes/ILMN_1745256.pheno

################## OPTION #2: read only probes with finished *.epi.qt.lm files ##################
if False: # The only reason to read from this dir, is that "you will only get the files you need"
	flag_probe_loading_method = "FASTEPISTASIS_RESULT_DIR"
	print "Reading only list of probes with finished *.epi.qt.lm files"
	#file_epi_qt_lm_summary = path_main_fastepistasis + "/fastEpi_out"
	file_epi_qt_lm_summary = path_epi_results
	
	### Globbing files: OBS: it is preferred to glob on *.epi.qt.lm files instead of *.epi.qt.lm.summary. 
	### This is because *.epi.qt.lm are DOWNSTREAM (writting by postFastEpistasis) from *.epi.qt.lm.summary (written by smpFastEpistasis)
	files_epi_qt_lm = glob.glob(file_epi_qt_lm_summary+"/*.epi.qt.lm") # reading all *.epi.qt.lm files
		# --> e.g. /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpi_1e-10_0108_212539/fastEpi_out/ILMN_1798104.pheno.epi.qt.lm
	### Stripping 'epi.qt.lm' from filename to produce the same 'probes' list as option#1. This will enable OPTION#1 and OPTION#2 to work in the same for loop ("for count, probe in enumerate(probes, start=1):")
	#probes_fastepistasis_result_dir = [elem.rstrip('epi.qt.lm') for elem in files_epi_qt_lm]
	probes = [elem.rstrip('epi.qt.lm') for elem in files_epi_qt_lm]
		# --> e.g. /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpi_1e-10_0108_212539/fastEpi_out/ILMN_1798104.pheno

print "DONE globbing LIST OF probes from path e.g.: %s" % probes[0] # * ..perhaps delete this.. *



###################################### Loops ######################################

################## EXAMPLE of *.epi.qt.lm file ##################
#  CHR                            SNP_A  CHR                            SNP_B         BETA        CHISQ       PVALUE
# ------------------------------------------------------------------------------------------------------------------
#   19                        rs2276470   21                        rs2834541     -0.11563     26.33892  2.86458E-07
# ...
# ...
# CONCLUSION: SKIP the first TWO lines!

## Open combined file for writing ##
f_combined = open(file_lm_combined, 'w') # write in MAIN DIR
header = "CHR\tSNP_A\tCHR\tSNP_B\tBETA\tCHISQ\tPVALUE\tPHENOTYPE" # NOTICE: the EXTRA FIELD (PHENOTYPE)
f_combined.write(header+"\n")

list_post_out_file_present = []
list_post_out_file_missing = []
list_post_out_file_empty = []
dict_probe_epistatic_count = {}

# SORTING on Illumina probe ID, e.g. "ILMN_2221"
probes_sorted = sorted(probes, key=lambda elem: os.path.basename(elem).split(".")[0])

#### NOW CONCATENATE ALL THE RESULTS #####
for count, probe in enumerate(probes_sorted, start=1):
	epi_out_prefix = os.path.basename(probe)
		# e.g. --> ILMN_1799744.pheno
	#epi_out_prefix = probe # *NEW*: to NOT take os.path.basename(probe). We need THE FULL PATH
		# e.g. --> /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpi_1e-10_0108_212539/fastEpi_out/ILMN_1799744.pheno
	illumina_probe_id = os.path.basename(epi_out_prefix).split(".")[0]
	
	post_out_file = path_epi_results + "/%s.epi.qt.lm" % epi_out_prefix # MAKE SURE TO MATCH this name to the name used in the "postFastEpistasis" section
	print "#%s/#%s: probe=%s - processing output and combining.." % (count, len(probes), epi_out_prefix)

	### DEBUGGING ### 
	# REMEMBER to also enable this break in the RESULT CONCATENATION STEP
	if test_n_probes:
		if count==test_n_probes: break
	#################

	### FLUSHING OUTPUT ###
	if count % 100 == 0:
		f_combined.flush()
		print "Flushing f_combined..."
	#################
	

	if not os.path.exists(post_out_file):
		print "post_out_file not found: %s" % post_out_file
		### append to list
		list_post_out_file_missing.append(illumina_probe_id)
	else:
		### READ and WRITE line-by-line
		post_file_epistatic_count = 0
		with open(post_out_file, 'r') as f:
			try:
				next(f) # SKIPPING THE FIRST TWO LINES!
				next(f)
			except StopIteration: # this expection is raised if the file has LESS THAN two lines (e.g. an empty file)
				print "WARNING: StopIteration execption raise when reading file={}. File has less than two lines.".format(post_out_file)
				list_post_out_file_empty.append(illumina_probe_id)
			for line in f:
				post_file_epistatic_count += 1
				fields = line.strip().split() # IMPORTANT TO USE strip() without argument to remove BOTH leading and trailing whitespace
				f_combined.write( "\t".join(fields) + "\t" + illumina_probe_id + "\n" )

			### append to list
			list_post_out_file_present.append(illumina_probe_id)

			### add count to dict
			dict_probe_epistatic_count[illumina_probe_id] = post_file_epistatic_count

f_combined.close()


################## Write probe counts - SORTED BY VALUE/COUNT ##################
with open(file_probe_epistatic_count, 'w') as f_probe_epistatic_count: # APPENDING TO FILE
	for illumina_probe_id in sorted(dict_probe_epistatic_count, key=dict_probe_epistatic_count.get, reverse=True):
		post_file_epistatic_count = dict_probe_epistatic_count[illumina_probe_id]
		f_probe_epistatic_count.write( "{probe}\t{count}\n".format(probe=illumina_probe_id, count=post_file_epistatic_count) )


################## Print status string ##################
status_string = """
len(list_post_out_file_present) = {}
len(list_post_out_file_missing) = {}
len(list_post_out_file_empty) = {}
""".format( 
	len(list_post_out_file_present), 
	len(list_post_out_file_missing),
	len(list_post_out_file_empty)
	)
print status_string

##### Writing to file_stat
with open(file_stat, 'w') as f_stat:
	f_stat.write(datetime.datetime.now().strftime("%a_%b_%d_%Y_%H%M%S")+'\n')
	f_stat.write( "flag_probe_loading_method=%s\n" % flag_probe_loading_method ) 
	f_stat.write( "Read list of probes from path: %s\n" % os.path.dirname(probes[0]) ) 
	f_stat.write(status_string + '\n')

##### Writing probes list
with open(file_probe_list, 'w') as f_probe_list:
	for probe in probes_sorted: # probes_sorted is sorted by illumina_probe_id, e.g. ILMN_2221
		# probe --> e.g. /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/lan-et-al_K562_width_1000_maf_5_q_OUTLIER_RM_epi1_1e-8/link_probes/ILMN_1343291.pheno
		illumina_probe_id = os.path.basename(probe).split(".")[0]
		# illumina_probe_id --> e.g. ILMN_2221
		f_probe_list.write( "{}\n".format(illumina_probe_id) )

with open(file_no_post_out_file, 'w') as f_no_post_file: # write in MAIN DIR
	for illumina_probe_id in list_post_out_file_missing:
		f_no_post_file.write("{}\tresult file did not exist\n".format(illumina_probe_id))
	for illumina_probe_id in list_post_out_file_empty:
		f_no_post_file.write("{}\tresult file was empty\n".format(illumina_probe_id))



print "DONE!"

###################################### TEST/EXPERIMENTING ######################################

################## SORTING PROBES/*.epi.qt.lm: by ILLUMINA ID ##################
#probes_sorted = [os.path.basename(elem) for elem in probes]
# test it:
#probes = ["/BLABLA/fastEpiCalc_epi1_1e-6_Thu_Jan_08_2015_194154/symlink_probes/ILMN_2221.pheno", "/BLABLA/fastEpi_out/ILMN_21.pheno.epi.qt.lm.summary"]



