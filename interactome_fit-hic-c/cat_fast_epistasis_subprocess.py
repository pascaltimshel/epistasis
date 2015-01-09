#!/usr/bin/env python2.7

import os
import subprocess
import glob
import datetime

import time
import argparse

###################################### USAGE ######################################

# python run_fast_epistasis_subprocess.py --job_name ISH.test --test_n_probes 3

###################################### IMPORTANT REMARKs ######################################
# THIS SCRIPT SHOULD BE THE DIRECTORY OUT THE OUTPUT SHOULD BE GENERATED
#	---> this is because FastEpistasis does *NOT* handle LONG PATH NAMES well
# e.g. "Cannot use filename so long: /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpiCalc_epi1_1e-6_Thu_Jan_08_2015_194154/symlink_probes/ILMN_1745256.pheno"
# The script should work anywhere on the filesystem because all the output file paths are relative to the WORKING DIRECTORY.

# 1) This script create a new directory in the directory from where the script is run (current working directory)
# 2) The script then changes its current working directory [os.chdir(path)] to this new directory
# 3) From here we do a lot of things:
# 	a) a symlink to the probe directory is created.
#	b) fastEpistasis is called
# 4) The *.epi.qt.lm file are concatenatted

# Update: *This script uses set files*

###############################################################################################

# preFastEpistasis --bfile /cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.SNPs781 --pheno fastEpiCal
# c_Fri_Dec_12_2014_122914/symlink_probes/ILMN_1678939.pheno --out fastEpiCalc_Fri_Dec_12_2014_122914/ILMN_1678939.pheno
# preFastEpistasis done
# smpFastEpistasis fastEpiCalc_Fri_Dec_12_2014_122914/ILMN_1678939.pheno.bin --method 0 --epi1 1e-6
# smpFastEpistasis done
# postFastEpistasis --pvalue --sort --index ILMN_1678939.pheno.epi.qt.lm.idx
# postFastEpistasis done


#Cannot open file ILMN_1678939.pheno.epi.qt.lm_001.bin
#os.chdir(path)


def test():
	try:
		FNULL = open(os.devnull, 'w')
		subprocess.Popen(["preFastEpistasis", "--help"], stdout=FNULL, stderr=subprocess.STDOUT)
		FNULL.close()
	except Exception as e:
		raise Exception("Could not find <PROGRAM> as executable on path. Please check that you have used 'use <PROGRAM>' Error msg: %s" % e)

test()


###################################### PARAMETERs ###################################### 
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--test_n_probes", type=int, help="INTERGER. Number of probes to test. [Parameter not required]") # defaults to None
arg_parser.add_argument("--path_fastEpi", required=True, help="fastEpi main dir, e.g. 'fastEpi_1e-10_0108_212539' NOT 'fastEpi_1e-10_0108_212539/fastEpi_out'.")
args = arg_parser.parse_args()

path_fastEpi = os.path.abspath(args.path_fastEpi) # IMPORTANT: converts to ABSOLUTE PATH, e.g. 'fastEpi_1e-10_0108_212539' --> '/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpi_1e-10_0108_212539'
test_n_probes = args.test_n_probes
if test_n_probes is not None:
	print "running in test mode: test_n_probes=%s" % test_n_probes
else:
	print "running in full mode (all probes)"


###################################### File - input ######################################
### Input


###################################### Preparing directories ######################################
#time_stamp = datetime.datetime.now().strftime("%a_%b_%d_%Y_%H%M%S")
time_stamp = datetime.datetime.now().strftime("%m%d_%H%M%S")
#work_dir = "fastEpiCalc_epi1_%s_%s" % (epi1, time_stamp)
work_dir = "%s_cat_%s" % (path_fastEpi, time_stamp)
os.mkdir(work_dir) # we assume that this dir does not exists
#### CHANGE DIR  - main work_dir ####
os.chdir(work_dir)
path_work_dir = os.getcwd()
print "Have now changed directory. Current working directory: %s" % os.getcwd()

###################################### LOADING PROBES ######################################

################## OPTION #1: read all probes in original dir ##################
if False: # TURNED OFF
	print "Reading list of all probes in original dir"
	## READ SYMLINK FROM "path_fastEpi" dir (result dir from fastEpi run)
	probe_dir_symlink = path_fastEpi + "/symlink_probes"
	### ** READING all FILENAMEs in probe_dir ** ###
	probes = glob.glob(probe_dir_symlink+"/*") # *OBS* this is DIRTY. When we change dir later this will not work. Thus the probes CANNOT be changed (which is actually ok.)
	# e.g. /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpiCalc_epi1_1e-6_Thu_Jan_08_2015_194154/symlink_probes/ILMN_1745256.pheno

################## OPTION #2: read only probes with finished *.epi.qt.lm.summary files ##################
if True: # TURNED ON
	print "Reading only list of probes with finished *.epi.qt.lm.summary files"
	file_epi_qt_lm_summary = path_fastEpi + "/fastEpi_out"
	### Globbing files
	epi_qt_lm_summaries = glob.glob(file_epi_qt_lm_summary+"/*.epi.qt.lm.summary") # reading all *.epi.qt.lm.summary files
	# --> e.g. /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpi_1e-10_0108_212539/fastEpi_out/ILMN_1798104.pheno.epi.qt.lm.summary
	### Stripping 'epi.qt.lm.summary' from filename to produce the same 'probes' list as option#1. This will enable OPTION#1 and OPTION#2 to work in the same for loop ("for count, probe in enumerate(probes, start=1):")
	probes = [elem.rstrip('epi.qt.lm.summary') for elem in epi_qt_lm_summaries]
	# --> e.g. /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpi_1e-10_0108_212539/fastEpi_out/ILMN_1798104.pheno

print "Reading LIST OF probes from path e.g.: %s" % probes[0] # * ..perhaps delete this.. *

################## SORTING PROBES/*.epi.qt.lm.summary ##################
#probes_sorted = [os.path.basename(elem) for elem in probes]
# test it:
#probes = ["/BLABLA/fastEpiCalc_epi1_1e-6_Thu_Jan_08_2015_194154/symlink_probes/ILMN_2221.pheno", "/BLABLA/fastEpi_out/ILMN_21.pheno.epi.qt.lm.summary"]
# SORTING on Illumina probe ID, e.g. "ILMN_2221"
probes_sorted = sorted(probes, key=lambda elem: os.path.basename(elem).split(".")[0])


###################################### Writing job_summary file ######################################

with open("job_summary.txt", 'w') as f:
	f.write("path_fastEpi: " + path_fastEpi + '\n')

###################################### Output ######################################

### Output
file_lm_combined = path_work_dir+"/results.epi.qt.lm.combined" # e.g. ???
file_no_post_out_file = path_work_dir+"/log_no_postFastEpistasis_file.txt"

###################################### Loops ######################################


################## EXAMPLE of *.epi.qt.lm file ##################
#  CHR                            SNP_A  CHR                            SNP_B         BETA        CHISQ       PVALUE
# ------------------------------------------------------------------------------------------------------------------
#   19                        rs2276470   21                        rs2834541     -0.11563     26.33892  2.86458E-07
# ...
# ...
# CONCLUSION: SKIP the first TWO lines!

## Open combined file for writing ##
f_combined = open(file_lm_combined, 'a') # write in MAIN DIR
header = "CHR\tSNP_A\tCHR\tSNP_B\tBETA\tCHISQ\tPVALUE\tPHENOTYPE" # NOTICE: the EXTRA FIELD (PHENOTYPE)
f_combined.write(header+"\n")

#### NOW CONCATENATE ALL THE RESULTS #####
for count, probe in enumerate(probes, start=1):
	#epi_out_prefix = os.path.basename(probe)
	epi_out_prefix = probe # *NEW*: to NOT take os.path.basename(probe). We need THE FULL PATH
		# e.g. --> /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpi_1e-10_0108_212539/fastEpi_out/ILMN_1799744.pheno
	illumina_probe_id = os.path.basename(epi_out_prefix).split(".")[0]
	
	post_out_file = "%s.epi.qt.lm" % epi_out_prefix # MAKE SURE TO MATCH this name to the name used in the "postFastEpistasis" section
	print "#%s/#%s: probe=%s - processing output and combining.." % (count, len(probes), epi_out_prefix)

	### DEBUGGING ### 
	# REMEMBER to also enable this break in the RESULT CONCATENATION STEP
	if test_n_probes:
		if count==test_n_probes: break
	#################

	if not os.path.exists(post_out_file):
		with open(file_no_post_out_file, 'a') as f_no_post_file: # write in MAIN DIR
			print "post_out_file not found: %s" % post_out_file
			f_no_post_file.write("%s\tpost_out_file did not exist\n" % post_out_file)
	else:
		### READ and WRITE line-by-line
		with open(post_out_file, 'r') as f:
			next(f) # SKIPPING THE FIRST TWO LINES!
			next(f)
			for line in f: 
				fields = line.strip().split() # IMPORTANT TO USE strip() without argument to remove BOTH leading and trailing whitespace
				f_combined.write( "\t".join(fields) + "\t" + illumina_probe_id + "\n" )

f_combined.close()

print "DONE!"

