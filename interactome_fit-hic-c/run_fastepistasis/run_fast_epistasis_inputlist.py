#!/usr/bin/env python2.7

import os
import sys
import subprocess
import glob
import datetime

import time
import argparse

import json


sys.path.insert(1, '/cvar/jhlab/snpsnap/snpsnap') # do not use sys.path.insert(0, 'somepath'). path[0], is the directory containing the script that was used to invoke the Python interpreter.
import pplaunch
import pphelper
import pplogger

###################################### USAGE ######################################

# python run_fast_epistasis_subprocess.py --file_job <FILE> --path_epi_results <PATH>

###################################### IMPORTANT REMARKs ######################################
### THIS SCRIPT *CHANGES DIRECTORY* TO THE PATH WHERE THE FastEpistasis OUTPUT IS WRITTEN
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


def LogArguments():
	# PRINT RUNNING DESCRIPTION 
	now = datetime.datetime.now()
	logger.critical( '# ' + ' '.join(sys.argv) )
	logger.critical( '# ' + now.strftime("%a %b %d %Y %H:%M") )
	logger.critical( '# CWD: ' + os.getcwd() )
	logger.critical( '# COMMAND LINE PARAMETERS SET TO:' )
	for arg in dir(args):
		if arg[:1]!='_':
			logger.critical( '# \t' + "{:<30}".format(arg) + "{:<30}".format(getattr(args, arg)) ) ## LOGGING


###################################### PARAMETERs ###################################### 
arg_parser = argparse.ArgumentParser()
### Read arguments
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--file_job", required=True, help="file with jobs. one line per job")
arg_parser.add_argument("--path_epi_results", required=True, help="path to the OUTPUT directory FastEpistasis calculations. E.g XXX/MAIN_DIR/epi")
arg_parser.add_argument("--file_log", help="path+filename. If a argument is passed a logger instance (from pplogger.py) will be created. The argument value should be the complete filename and path of the logging file.")
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
file_log = args.file_log

### Config file
file_json_config = os.path.dirname(path_epi_results) + "/config.json" # OBS: we use the BACKSLASH ('/') because we are sure that any additional backslash from path_epi_results is removed by .abspath()


FastEpistasis_bin = "/cvar/jhlab/timshel/bin/FastEpistasis-Linux-x86_64-2.07/bin"

###################################### SETUP logging ######################################
timestamp = datetime.datetime.now().strftime("%Y_%m_%d-%H.%M.%S") # e.g. '2014_10_13-14.08.47'
current_script_name = os.path.basename(__file__).replace('.py','')

if file_log:
	## Remember the following: os.path.dirname(filename) + os.path.basename(filename) == filename
	log_dir = os.path.dirname(file_log)
	if not os.path.exists(log_dir): 
		raise Exception("log_dir does not exists. log_dir=%s" % log_dir)
	#log_name = timestamp + '-' + current_script_name # e.g. 2014_10_13-14.08.47-wrapper_bsub_prune_sets.py. Uses timestamp!
	log_name = os.path.basename(file_log)
	logger = pplogger.Logger(name=log_name, log_dir=log_dir, log_format=1, enabled=True, log2stream=True, log2file=True).get()
	def handleException(excType, excValue, traceback, logger=logger):
		logger.error("Logging an uncaught exception", exc_info=(excType, excValue, traceback))
	#### TURN THIS ON OR OFF: must correspond to enabled='True'/'False'
	sys.excepthook = handleException
	logger.info( "INSTANTIATION NOTE: placeholder" )
else:
	logger = pplogger.Logger(enabled=False).get()
###########################################################################################

###################################### Run initial functions ######################################
LogArguments()



### OBS: notice the big difference: 
# >>> os.path.abspath('/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/') --> CONCLUSION: .abspath() removes trailing '/'
# '/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci'
# >>> os.path.dirname('/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci') 
# '/cvar/jhlab/timshel/egcut' ---> LOOK HERE
# >>> os.path.dirname('/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/') 
# '/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci' ---> LOOK HERE

### Generating file names for output
#file_out = os.path.join(os.path.dirname(file_job), os.path.splitext(os.path.basename(file_job))[0] + "_worker.txt") # extract basename without extenstion. this works even if the file does not have an extension
#file_bad_out = os.path.join(os.path.dirname(file_job), "fail_jobs_" + os.path.splitext(os.path.basename(file_job))[0] + "_worker.txt") 
#sys.stdout = open(file_out, 'w') # Redirect output
#sys.stdrr = open(file_out, 'w') # Redirect errors output (so it does not end up in the bsub out)


###################################### Starting ######################################
#logger.info( "#file_out: {}".format(file_out) )
#logger.info( "#file_bad_out: {}".format(file_bad_out) )

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
	logger.info( "%s:%s" % (key, value) )

################## Reading jobs ##################
logger.info( "Reading jobs to process..." )
jobs = read_jobs(file_job)
logger.info( "Number of jobs: %s" % len(jobs) )

make_subprocess_calls = True
###################################### Changing dir to EPISTASIS file ######################################
#### make dirs 
#subprocess_out_dir = "subprocess_out"
#os.mkdir(subprocess_out_dir)

#### CHANGE DIR - fastEpi_out ####
os.chdir(path_epi_results)
logger.info( "Have now changed directory. Current working directory: %s" % os.getcwd() )

###################################### Loops ######################################
time_run_fast_epistasis_start = time.time() # start time
time_run_fast_epistasis_list = []

for count, probe in enumerate(jobs, start=1):
	time_probe_start = time.time()

	### OBS: dirty-code. MODIFYING probe inside for-loop
	probe = probe+".pheno"

	probe_full_path = path_files + "/{probe}".format(probe=probe) # OBS: THIS IS EXTREMELY IMPORTANT!
	epi_out_prefix = os.path.basename(probe) # E.g. ILMN_1678939.pheno
	logger.info( "#%s/#%s: job=%s" % (count, len(jobs), epi_out_prefix) )
	

	################## preFastEpistasis ##################
	time_preFastEpistasis_start = time.time()
	cmd = "%s/preFastEpistasis --bfile %s --pheno %s --set %s --out %s" % (FastEpistasis_bin, file_bim, probe_full_path, file_set, epi_out_prefix) # # Output name will be probe
	#cmd = "preFastEpistasis --bfile %s --pheno %s --set %s --out %s" % (file_bim, probe_full_path, file_set, epi_out_prefix) # # Output name will be probe
	logger.info( cmd )

	# subprocess_file_out = os.devnull
	# subprocess_file_out = "../%s/%s.preFastEpistasis.subprocess.out" % (subprocess_out_dir, epi_out_prefix)
	# with open(subprocess_file_out, 'w') as fsub:
	# 	p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
	# 	p.wait()

	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	p.wait()

	time_preFastEpistasis_elapsed = time.time() - time_preFastEpistasis_start
	if p.returncode == 0:
		logger.info( "preFastEpistasis done: %s" % time_preFastEpistasis_elapsed )
		#(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		#logger.info( stdoutdata )
	else:
		logger.critical( "preFastEpistasis FAILED" )
		(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		logger.critical( stdoutdata )
	




	################## smpFastEpistasis ##################
	time_smpFastEpistasis_start = time.time()
	# Remember: smpFastEpistasis does not take an "--out" argument
	bin_file = "%s.bin" % (epi_out_prefix, ) # e.g. ILMN_1678939.pheno.bin
	cmd = "%s/smpFastEpistasis %s --method 0 --epi1 %s --epi2 %s --nthreads %s" % (FastEpistasis_bin, bin_file, epi1, epi2, nthreads)
	#cmd = "smpFastEpistasis %s --method 0 --epi1 %s --epi2 %s --nthreads %s" % (bin_file, epi1, epi2, nthreads)
	logger.info( cmd )

	# subprocess_file_out = os.devnull
	# subprocess_file_out = "../%s/%s.smpFastEpistasis.subprocess.out" % (subprocess_out_dir, epi_out_prefix)
	# with open(subprocess_file_out, 'w') as fsub:
	# 	p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
	# 	p.wait()

	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	p.wait()

	time_smpFastEpistasis_elapsed = time.time() - time_smpFastEpistasis_start
	if p.returncode == 0:
		logger.info( "smpFastEpistasis done: %s" % time_smpFastEpistasis_elapsed )
		#(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		#logger.info( stdoutdata )
	else:
		logger.critical( "smpFastEpistasis FAILED" )
		(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		logger.critical( stdoutdata )
	# ----> this will create a index file like "<probe>.epi.qt.lm.idx"






	################## postFastEpistasis ##################
	index_file = "%s.epi.qt.lm.idx" % epi_out_prefix # e.g. ILMN_1678939.pheno.epi.qt.lm.idx
	post_out_file = "%s.epi.qt.lm" % epi_out_prefix # e.g. ILMN_1678939.pheno.epi.qt.lm - NOTE THAT USING "--out <FILENAME>" sets the EXACT FILENAME - no extensions will be added
	cmd = "%s/postFastEpistasis --pvalue --sort --index %s --out %s" % (FastEpistasis_bin, index_file, post_out_file)
	#cmd = "postFastEpistasis --pvalue --sort --index %s --out %s" % (index_file, post_out_file)
	logger.info( cmd )

	# subprocess_file_out = os.devnull
	# subprocess_file_out = "../%s/%s.postFastEpistasis.subprocess.out" % (subprocess_out_dir, epi_out_prefix)
	# with open(subprocess_file_out, 'w') as fsub:
	# 	p = subprocess.Popen(cmd, stdout=fsub, stderr=subprocess.STDOUT, shell=True)
	# 	p.wait()

	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	p.wait()

	if p.returncode == 0:
		logger.info( "postFastEpistasis done" )
		#(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		#logger.info( stdoutdata )
	else:
		logger.critical( "postFastEpistasis FAILED" )
		(stdoutdata, stderrdata) = p.communicate() # note stderrdata is empty
		logger.critical( stdoutdata )

	################## Cleaning ##################
	### Removing .bin file
	if os.path.exists(bin_file):
		bin_file_size_mb = os.path.getsize(bin_file)/(1024*1024.0)
		logger.info( "Removing .bin file: %s (size=%s MB)" % (bin_file, bin_file_size_mb) )
		os.remove(bin_file)
	else:
		logger.critical( "WARNING: no .bin file found!" )

	core_dumps = glob.glob("core.*")
	if core_dumps: # if not empty
		for core_dump in core_dumps:
			core_dump_file_size_mb = os.path.getsize(core_dump)/(1024*1024.0)
			logger.critical( "WARNING: found core dump: %s (size=%s MB)" % (core_dump, core_dump_file_size_mb) )
			os.remove(core_dump)

	################## Finishing ##################
	time_probe_elapsed = time.time() - time_probe_start
	logger.info( "RUNTIME FOR PROBE: %s" % time_probe_elapsed )
	# calculate mean time
	time_run_fast_epistasis_list.append(time_probe_elapsed) # type(time.time()) --> type 'float'
	time_run_fast_epistasis_list_avg = sum(time_run_fast_epistasis_list)/float(len(time_run_fast_epistasis_list))
	logger.info( "AVG RUNTIME FOR PROBE: %s s (%s min)" % (time_run_fast_epistasis_list_avg, time_run_fast_epistasis_list_avg/60) )

	### More timing
	time_run_fast_epistasis_elapsed = time.time() - time_run_fast_epistasis_start
	logger.info( "TIME ELAPSED FOR PROGRAM: %s s (%s min)" % (time_run_fast_epistasis_elapsed, time_run_fast_epistasis_elapsed/60) )

logger.info( "DONE RUNNING PROGRAM" )