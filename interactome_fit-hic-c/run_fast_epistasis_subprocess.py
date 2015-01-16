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
arg_parser.add_argument("--job_name", required=True, help="Job name. This name will be written in the job_summary.txt file")
arg_parser.add_argument("--test_n_probes", type=int, help="INTERGER. Number of probes to test. [Parameter not required]") # defaults to None
args = arg_parser.parse_args()

job_name = args.job_name
test_n_probes = args.test_n_probes
if test_n_probes is not None:
	print "running in test mode: test_n_probes=%s" % test_n_probes
else:
	print "running in full mode (all probes)"

print "job_name: %s" % job_name

###################################### File - input ######################################
### Input
#bfile = "/cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean" # DO NOT ADD EXTENSION to file
bfile = "/cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.maf5" # DO NOT ADD EXTENSION to file
probe_dir = "/cvar/jhlab/timshel/egcut/ETypes_probes_norm_peer/phenofile_log2_k50.top50_mean_top50_var_refseq"
#file_set = "/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/10000_snppool_hIMR90_q_1e-10/snp_sets/set_AB.txt"
file_set = "/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/10000_snppool_hIMR90_q_1e-07/snp_sets/set_AB.txt"


make_subprocess_calls = True # this controls whether or not to make subprocess calls.

###################################### Significance threshold: epi1 ######################################
#epi1 = "1" # ALL OUTPUT!
#epi1 = "1e-6" # string to set the significance thresshold
epi1 = "1e-10" # string to set the significance thresshold
epi2 = epi1


###################################### Preparing directories ######################################
#time_stamp = datetime.datetime.now().strftime("%a_%b_%d_%Y_%H%M%S")
time_stamp = datetime.datetime.now().strftime("%m%d_%H%M%S")
#work_dir = "fastEpiCalc_epi1_%s_%s" % (epi1, time_stamp)
work_dir = "fastEpi_%s_%s" % (epi1, time_stamp)
os.mkdir(work_dir) # we assume that this dir does not exists
#### CHANGE DIR  - main work_dir ####
os.chdir(work_dir)
path_work_dir = os.getcwd()
print "Have now changed directory. Current working directory: %s" % os.getcwd()

## create symlink to probes
probe_dir_symlink = path_work_dir + "/symlink_probes"
os.symlink(probe_dir, probe_dir_symlink)
### ** READING all FILENAMEs in probe_dir ** ###
probes = glob.glob(probe_dir_symlink+"/*") # *OBS* this is DIRTY. When we change dir later this will not work. Thus the probes CANNOT be changed (which is actually ok.)

###################################### Writing job_summary file ######################################

with open("job_summary.txt", 'w') as f:
	f.write("job_name: " + job_name + '\n')
	f.write("bfile: " + bfile + '\n')
	f.write("probe_dir: " + probe_dir + '\n')
	f.write("file_set: " + file_set + '\n')
	f.write("epi1: " + epi1 + '\n')

###################################### Changing dir to EPISTASIS file ######################################
#### make dirs 
fastEpi_dir = "fastEpi_out"
subprocess_out_dir = "subprocess_out"
#fastEpi_bin = "fastEpi_bin"
os.mkdir(fastEpi_dir)
os.mkdir(subprocess_out_dir)
#### CHANGE DIR - fastEpi_out ####
os.chdir(fastEpi_dir)
print "Have now changed directory. Current working directory: %s" % os.getcwd()

###################################### Output ######################################

### Output
file_lm_combined = path_work_dir+"/results.epi.qt.lm.combined"
file_no_post_out_file = path_work_dir+"/log_no_postFastEpistasis_file.txt"

###################################### Loops ######################################

for count, probe in enumerate(probes, start=1):
	time_probe_start = time.time()

	probe_full_path = probe
	epi_out_prefix = os.path.basename(probe)
	print "#%s/#%s: probe=%s" % (count, len(probes), epi_out_prefix)
	
	### DEBUGGING ### 
	# REMEMBER to also enable this break in the RESULT CONCATENATION STEP
	if test_n_probes:
		if count==test_n_probes: break
	#################

	################## preFastEpistasis ##################
	time_preFastEpistasis_start = time.time()
	cmd = "preFastEpistasis --bfile %s --pheno %s --set %s --out %s" % (bfile, probe_full_path, file_set, epi_out_prefix) # # Output name will be probe
	print cmd
	subprocess_file_out = "../%s/%s.preFastEpistasis.subprocess.out" % (subprocess_out_dir, epi_out_prefix)
	with open(subprocess_file_out, 'w') as fsub:
		if make_subprocess_calls:
			p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
			p.wait()
	time_preFastEpistasis_elapsed = time.time() - time_preFastEpistasis_start
	print "preFastEpistasis done: %s" % time_preFastEpistasis_elapsed
	
	################## smpFastEpistasis ##################
	time_smpFastEpistasis_start = time.time()
	# Remember: smpFastEpistasis does not take an "--out" argument
	bin_file = "%s.bin" % (epi_out_prefix, ) # e.g. ILMN_1678939.pheno.bin
	cmd = "smpFastEpistasis %s --method 0 --epi1 %s --epi2 %s" % (bin_file, epi1, epi2)
	print cmd
	subprocess_file_out = "../%s/%s.smpFastEpistasis.subprocess.out" % (subprocess_out_dir, epi_out_prefix)
	with open(subprocess_file_out, 'w') as fsub:
		if make_subprocess_calls:
			p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
			p.wait()
	time_smpFastEpistasis_elapsed = time.time() - time_smpFastEpistasis_start
	print "smpFastEpistasis done: %s" % time_smpFastEpistasis_elapsed
	# ----> this will create a index file like "<probe>.epi.qt.lm.idx"

	################## postFastEpistasis ##################
	index_file = "%s.epi.qt.lm.idx" % epi_out_prefix # e.g. ILMN_1678939.pheno.epi.qt.lm.idx
	post_out_file = "%s.epi.qt.lm" % epi_out_prefix # e.g. ILMN_1678939.pheno.epi.qt.lm - NOTE THAT USING "--out <FILENAME>" sets the EXACT FILENAME - no extensions will be added
	cmd = "postFastEpistasis --pvalue --sort --index %s --out %s" % (index_file, post_out_file)
	print cmd
	subprocess_file_out = "../%s/%s.postFastEpistasis.subprocess.out" % (subprocess_out_dir, epi_out_prefix)
	with open(subprocess_file_out, 'w') as fsub:
		if make_subprocess_calls:
			p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
			p.wait()
	print "postFastEpistasis done"

	################## Cleaning ##################
	### Removing .bin file
	bin_file_size_mb = os.path.getsize(bin_file)/(1024*1024.0)
	print "Removing .bin file: %s (size=%s)" % (bin_file, bin_file_size_mb)
	os.remove(bin_file)

	################## Finishing ##################
	time_probe_elapsed = time.time() - time_probe_start
	print "RUNTIME FOR PROBE: %s" % time_probe_elapsed


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
	epi_out_prefix = os.path.basename(probe)
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
				f_combined.write( "\t".join(fields) + "\t" + epi_out_prefix + "\n" )

f_combined.close()

print "DONE!"


################## BEFORE reading line by line ##################
	# if not os.path.exists(post_out_file):
	# 	with open(file_no_post_out_file, 'a') as f_no_post_file: # write in MAIN DIR
	# 		print "post_out_file not found: %s" % post_out_file
	# 		f_no_post_file.write("%s\tpost_out_file did not exist\n" % post_out_file)
	# else:
	# 	with open(post_out_file, 'r') as f: 
	# 		lines = f.readlines()
	# 		for line in lines[2:]: # SKIPPING THE FIRST TWO LINES!
	# 			fields = line.strip().split() # IMPORTANT TO USE strip() without argument to remove BOTH leading and trailing whitespace
	# 			f_combined.write( "\t".join(fields) + "\t" + epi_out_prefix + "\n" )
