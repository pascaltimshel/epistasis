#!/usr/bin/env python2.7

import os
import sys
import subprocess
import glob
import datetime

import time
import argparse

import json

###################################### USAGE ######################################

# python run_fast_epistasis_subprocess.py --file_job <FILE> --path_epi_results <PATH>

###################################### IMPORTANT REMARKs ######################################
### *ALL STDOUT IS REDIRECTED TO A FILE*

# THIS SCRIPT SHOULD BE THE DIRECTORY OUT THE OUTPUT SHOULD BE GENERATED
#	---> this is because FastEpistasis does *NOT* handle LONG PATH NAMES well
# e.g. "Cannot use filename so long: /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/fastEpiCalc_epi1_1e-6_Thu_Jan_08_2015_194154/symlink_probes/ILMN_1745256.pheno"

# Update: *This script uses set files*



def read_jobs(file_job):
	with open(file_job, 'r') as f:
		jobs = f.read().splitlines()
		return jobs

def read_config_file():
	with open(file_json_config, 'r') as f_json:
		config_dict = json.load(f_json)
	return config_dict

###################################### PARAMETERs ###################################### 
arg_parser = argparse.ArgumentParser()
### Read arguments
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--file_job", help="file with jobs. one line per job")
arg_parser.add_argument("--path_epi_results", help="path to the OUTPUT directory FastEpistasis calculations. E.g XXX/MAIN_DIR/epi")
# path_epi_results/
	# epi/
	# jobs/
	# logs/
	# POTENTIAL: subprocess_out/
	# link_bim
	# link_probes
	# link_set
 	# config.json
args = arg_parser.parse_args()

################## Arguments ##################
file_job = args.file_job ## this file MUST exists. Otherwise script will crash
path_epi_results = os.path.abspath(args.path_epi_results)

### Config file
file_json_config = os.path.dirname(path_epi_results) + "/config.json" # OBS: we use the BACKSLASH ('/') because we are sure that any additional backslash from path_epi_results is removed by .abspath()



### OBS: notice the big difference: 
# >>> os.path.abspath('/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/') --> CONCLUSION: .abspath() removes trailing '/'
# '/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci'
# >>> os.path.dirname('/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci') 
# '/cvar/jhlab/timshel/egcut' ---> LOOK HERE
# >>> os.path.dirname('/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/') 
# '/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci' ---> LOOK HERE


### Generating file names for output
file_out = os.path.join(os.path.dirname(file_job), os.path.splitext(os.path.basename(file_job))[0] + "_worker.txt") # extract basename without extenstion. this works even if the file does not have an extension
file_bad_out = os.path.join(os.path.dirname(file_job), "fail_jobs_" + os.path.splitext(os.path.basename(file_job))[0] + "_worker.txt") 
#sys.stdout = open(file_out, 'w') # Redirect output
#sys.stdrr = open(file_out, 'w') # Redirect errors output (so it does not end up in the bsub out)


###################################### Starting ######################################
time_start = time.time()
print "#file_out: {}".format(file_out)
print "#file_bad_out: {}".format(file_bad_out)

################## Reading config file ##################
config_dict = read_config_file()

file_bim = config_dict["link_bim"]
# path_files = config_dict["link_probes"]
# file_set = config_dict["link_set"]
#file_bim = os.path.dirname(path_epi_results) + "/link_bim" # SYMLINK TO THE BIM FILE DOES NOT WORK! FastEpistasis will think the files are /XXX/link_bim.bim, /XXX/link_bim.fam, /XXX/link_bim.
path_files = os.path.dirname(path_epi_results) + "/link_probes"
file_set = os.path.dirname(path_epi_results) + "/link_set"
epi1 = config_dict["epi1"]
epi2 = config_dict["epi2"]
nthreads = config_dict["nthreads"]

for key, value in config_dict.items():
	print "%s:%s" % (key, value)

################## Reading jobs ##################
print "#Reading jobs to process..."
jobs = read_jobs(file_job)
print "#Number of jobs: %s" % len(jobs)

make_subprocess_calls = True
###################################### Changing dir to EPISTASIS file ######################################
#### make dirs 
#subprocess_out_dir = "subprocess_out"
#os.mkdir(subprocess_out_dir)

#### CHANGE DIR - fastEpi_out ####
os.chdir(path_epi_results)
print "Have now changed directory. Current working directory: %s" % os.getcwd()

FastEpistasis_bin = "/cvar/jhlab/timshel/bin/FastEpistasis-Linux-x86_64-2.07/bin"
###################################### Loops ######################################

for count, probe in enumerate(jobs, start=1):
	time_probe_start = time.time()

	### OBS: dirty-code. MODIFYING probe inside for-loop
	probe = probe+".pheno"

	probe_full_path = path_files + "/{probe}".format(probe=probe) # OBS: THIS IS EXTREMELY IMPORTANT!
	epi_out_prefix = os.path.basename(probe) # E.g. ILMN_1678939.pheno
	print "#%s/#%s: job=%s" % (count, len(jobs), epi_out_prefix)
	

	################## preFastEpistasis ##################
	time_preFastEpistasis_start = time.time()
	cmd = "%s/preFastEpistasis --bfile %s --pheno %s --set %s --out %s" % (FastEpistasis_bin, file_bim, probe_full_path, file_set, epi_out_prefix) # # Output name will be probe
	#cmd = "preFastEpistasis --bfile %s --pheno %s --set %s --out %s" % (file_bim, probe_full_path, file_set, epi_out_prefix) # # Output name will be probe
	print cmd

	# subprocess_file_out = os.devnull
	# subprocess_file_out = "../%s/%s.preFastEpistasis.subprocess.out" % (subprocess_out_dir, epi_out_prefix)
	# with open(subprocess_file_out, 'w') as fsub:
	# 	p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
	# 	p.wait()

	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	p.wait()

	time_preFastEpistasis_elapsed = time.time() - time_preFastEpistasis_start
	if p.returncode == 0:
		print "preFastEpistasis done: %s" % time_preFastEpistasis_elapsed
		(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		print stdoutdata
	else:
		print "preFastEpistasis FAILED"
		(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		print stdoutdata
	




	################## smpFastEpistasis ##################
	time_smpFastEpistasis_start = time.time()
	# Remember: smpFastEpistasis does not take an "--out" argument
	bin_file = "%s.bin" % (epi_out_prefix, ) # e.g. ILMN_1678939.pheno.bin
	cmd = "%s/smpFastEpistasis %s --method 0 --epi1 %s --epi2 %s --nthreads %s" % (FastEpistasis_bin, bin_file, epi1, epi2, nthreads)
	#cmd = "smpFastEpistasis %s --method 0 --epi1 %s --epi2 %s --nthreads %s" % (bin_file, epi1, epi2, nthreads)
	print cmd

	# subprocess_file_out = os.devnull
	# subprocess_file_out = "../%s/%s.smpFastEpistasis.subprocess.out" % (subprocess_out_dir, epi_out_prefix)
	# with open(subprocess_file_out, 'w') as fsub:
	# 	p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
	# 	p.wait()

	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	p.wait()

	time_smpFastEpistasis_elapsed = time.time() - time_smpFastEpistasis_start
	if p.returncode == 0:
		print "smpFastEpistasis done: %s" % time_smpFastEpistasis_elapsed
		(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		print stdoutdata
	else:
		print "smpFastEpistasis FALIED"
		(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		print stdoutdata
	# ----> this will create a index file like "<probe>.epi.qt.lm.idx"






	################## postFastEpistasis ##################
	index_file = "%s.epi.qt.lm.idx" % epi_out_prefix # e.g. ILMN_1678939.pheno.epi.qt.lm.idx
	post_out_file = "%s.epi.qt.lm" % epi_out_prefix # e.g. ILMN_1678939.pheno.epi.qt.lm - NOTE THAT USING "--out <FILENAME>" sets the EXACT FILENAME - no extensions will be added
	cmd = "%s/postFastEpistasis --pvalue --sort --index %s --out %s" % (FastEpistasis_bin, index_file, post_out_file)
	#cmd = "postFastEpistasis --pvalue --sort --index %s --out %s" % (index_file, post_out_file)
	print cmd

	# subprocess_file_out = os.devnull
	# subprocess_file_out = "../%s/%s.postFastEpistasis.subprocess.out" % (subprocess_out_dir, epi_out_prefix)
	# with open(subprocess_file_out, 'w') as fsub:
	# 	p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
	# 	p.wait()

	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	p.wait()

	if p.returncode == 0:
		print "postFastEpistasis done"
		(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		print stdoutdata
	else:
		print "postFastEpistasis FAILED"
		(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		print stdoutdata

	################## Cleaning ##################
	### Removing .bin file
	if os.path.exists(bin_file):
		bin_file_size_mb = os.path.getsize(bin_file)/(1024*1024.0)
		print "Removing .bin file: %s (size=%s MB)" % (bin_file, bin_file_size_mb)
		os.remove(bin_file)
	else:
		print "WARNING: no .bin file found!"

	################## Finishing ##################
	time_probe_elapsed = time.time() - time_probe_start
	print "RUNTIME FOR PROBE: %s" % time_probe_elapsed


print "DONE!"