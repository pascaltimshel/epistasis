#!/usr/bin/env python2.7

import os
import subprocess
import glob
import datetime


###################################### IMPORTANT REMARKs ######################################
# THIS SCRIPT IS ORIGINALLY ONLY MEANT TO BE RUN IN THE DIRECTORY: /cvar/jhlab/timshel/test_fastEpistasis/hemani_test
# However, it should work anywhere on the filesystem because all the output file paths are relative to the WORKING DIRECTORY.

# 1) This script create a new directory in the directory from where the script is run (current working directory)
# 2) The script then changes its current working directory [os.chdir(path)] to this new directory
# 3) From here we do a lot of things:
# 	a) a symlink to the probe directory is created.
#	b) fastEpistasis is called
# 4) The *.epi.qt.lm file are concatenatted

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

###################################### File - input ######################################
### Input
bfile = "/cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.SNPs781"
probe_dir = "/cvar/jhlab/timshel/egcut/ETypes_probes_norm_peer/phenofile_log2_k50.top50_mean_top50_var_refseq_hemani_probes_unique_102"


make_subprocess_calls = True # this controls whether or not to make subprocess calls.
epi1 = "1e-5" # string to set the significance thresshold

###################################### Preparing directories ######################################
time_stamp = datetime.datetime.now().strftime("%a_%b_%d_%Y_%H%M%S")
work_dir = "fastEpiCalc_%s_%s" % (epi1, time_stamp)
os.mkdir(work_dir) # we assume that this dir does not exists
#### CHANGE DIR  - main work_dir ####
os.chdir(work_dir)
path_work_dir = os.getcwd()
print "Have now changed directory. Current working directory: %s" % os.getcwd()

## create symlink to probes
probe_dir_symlink = path_work_dir + "/symlink_probes"
os.symlink(probe_dir, probe_dir_symlink)
probes = glob.glob(probe_dir_symlink+"/*") # *OBS* this is DIRTY. When we change dir later this will not work. Thus the probes CANNOT be changed (which is actually ok.)

## make dir for Epistasis output
fastEpi_dir = "fastEpi_out"
os.mkdir(fastEpi_dir)
#### CHANGE DIR - fastEpi_out ####
os.chdir(fastEpi_dir)
print "Have now changed directory. Current working directory: %s" % os.getcwd()

###################################### Output ######################################

### Output
file_lm_combined = path_work_dir+"/results.epi.qt.lm.combined"
file_no_post_out_file = path_work_dir+"/log_no_postFastEpistasis_file.txt"

###################################### Loops ######################################

for count, probe in enumerate(probes, start=1):
	probe_full_path = probe
	epi_out_prefix = os.path.basename(probe)
	print "#%s/#%s: probe=%s" % (count, len(probes), epi_out_prefix)
	
	#if count==10: break

	################## preFastEpistasis ##################
	cmd = "preFastEpistasis --bfile %s --pheno %s --out %s" % (bfile, probe_full_path, epi_out_prefix) # # Output name will be probe
	print cmd
	subprocess_file_out = "%s.preFastEpistasis.subprocess.out" % epi_out_prefix
	with open(subprocess_file_out, 'w') as fsub:
		if make_subprocess_calls:
			p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
			p.wait()
	print "preFastEpistasis done"
	
	################## smpFastEpistasis ##################
	# Remember: smpFastEpistasis does not take an "--out" argument
	bin_file = "%s.bin" % (epi_out_prefix, ) # e.g. ILMN_1678939.pheno.bin
	cmd = "smpFastEpistasis %s --method 0 --epi1 %s" % (bin_file, epi1)
	print cmd
	subprocess_file_out = "%s.smpFastEpistasis.subprocess.out" % epi_out_prefix
	with open(subprocess_file_out, 'w') as fsub:
		if make_subprocess_calls:
			p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
			p.wait()
	print "smpFastEpistasis done"
	# ----> this will create a index file like "<probe>.epi.qt.lm.idx"

	################## postFastEpistasis ##################
	index_file = "%s.epi.qt.lm.idx" % epi_out_prefix # e.g. ILMN_1678939.pheno.epi.qt.lm.idx
	post_out_file = "%s.epi.qt.lm" % epi_out_prefix # e.g. ILMN_1678939.pheno.epi.qt.lm - NOTE THAT USING "--out <FILENAME>" sets the EXACT FILENAME - no extensions will be added
	cmd = "postFastEpistasis --pvalue --sort --index %s --out %s" % (index_file, post_out_file)
	print cmd
	subprocess_file_out = "%s.postFastEpistasis.subprocess.out" % epi_out_prefix
	with open(subprocess_file_out, 'w') as fsub:
		if make_subprocess_calls:
			p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
			p.wait()
	print "postFastEpistasis done"


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

	#if count==10: break

	if not os.path.exists(post_out_file):
		with open(file_no_post_out_file, 'a') as f_no_post_file: # write in MAIN DIR
			print "post_out_file not found: %s" % post_out_file
			f_no_post_file.write("%s\tpost_out_file did not exist\n" % post_out_file)
	else:
		with open(post_out_file, 'r') as f: 
			lines = f.readlines()
			for line in lines[2:]: # SKIPPING THE FIRST TWO LINES!
				fields = line.strip().split() # IMPORTANT TO USE strip() without argument to remove BOTH leading and trailing whitespace

				f_combined.write( "\t".join(fields) + "\t" + epi_out_prefix + "\n" )

f_combined.close()

