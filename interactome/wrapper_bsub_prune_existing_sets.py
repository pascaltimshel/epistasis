#!/usr/bin/env python2.7

import sys
import glob
import os

import datetime
import time
import argparse

sys.path.insert(1, '/cvar/jhlab/snpsnap/snpsnap') # do not use sys.path.insert(0, 'somepath'). path[0], is the directory containing the script that was used to invoke the Python interpreter.
import pplaunch
import pphelper
import pplogger

import re
import subprocess
import logging

###################################### Description ######################################
## IMPORATNT: This script has three important variables to check:
## 1) interaction_width		for interactome parameter
## 2) path_snp_sets 		for input sets. 
## 							**OBS: this path must ONLY contain snp_sets files. OTHERWISE THE PROGRAM MAY FAIL!
##							consider using the following to check the folder integrity: ls -1af | grep -v 'interaction.*txt'						
##							TODO: make a 'static list' of files to process instead
## 3) path_snp_results	 	for caching property
## The path_snp_sets is combined with interaction_width to search for input files.
## CACHING PROPERTY
## The script searches the path_snp_results for previous results
## Thus you need to check that the expected OUTPUT files are still as described in this program. SEE find_completed_results()
## ***TODO: find out how to check if the PLINK job completed SUCCESSFULLY (try to find pattern in log file)
## ADDITIONAL INFO:
## File are sensitive for extensions, i.e. must match interaction_{}_A.txt
## CONTROL PARAMETERS:
## 1) n_jobs_per_bsub:
## 		number of jobs listed in the text file to be parsed to worker_multiprocess_prune_sets.py
##		this number should be a trade-off between paralizing efficiently (n_jobs_per_bsub low; many bsub_jobs) and not generating too many bsub_jobs.
## OUTPUT:
## path_jobs		job files are written to this path.
##					job files ae structured as <DATE>_<TIMESTAMP>_JOBS_<n_jobs_in_batch>_<job_no>.txt, e.g.
##						2014_10_13-12-56-38_JOBS_100_1.txt
##						2014_10_13-12-56-38_JOBS_100_2.txt
##						...
## log_dir			this path is dependent on the width parameter
## bsub_output 		uses SAME path as log_dir


###################################### NEW FUNCITONS ######################################


def validate_interaction_file_pairs():
	"""Function to validate that all interaction contains both files. This is good practice to ensure that all files exists"""
	unpaired_interactions = []
	for no, interaction in enumerate(snp_interactions, start=1):
		logger.info( "validating interaction number %s" % no )
		A_file = path_snp_sets + "/interaction_{}_A.txt".format(interaction)
		B_file = path_snp_sets + "/interaction_{}_B.txt".format(interaction)
		if not (os.path.exists(A_file) and os.path.exists(B_file)): # BOTH file must exists. DeMorgan's Law: !(AandB) <=> (!Aor!B)
			unpaired_interactions.append(interaction)
	if unpaired_interactions:
		raise Exception("Detected unpaired interaction files: %s" % "\n".join(unpaired_interactions))
	else:
		logger.info( "All interaction files validated and paired." )

def find_completed_results(snp_interactions):
	#### ATTEMPT TO READ LOG FILE - pattern did work... ###
	# PLINK1.7
	# Analysis finished: Sun Oct 12 17:49:43 2014
	# PLINK1.9
	# End time: Sun Oct 12 17:53:38 2014
	# success_string_plink_1_7 = 'Analysis finished' 
	# success_string_plink_1_9 = 'End time'
	# for interaction in snp_interactions:
	# 	A_plink_log_file = path_snp_results + "/interaction_{}_A.log".format(interaction) # /cvar/jhlab/timshel/egcut/interactome/5000_pruned/snp_sets/interaction_1_A.log
	# 	B_plink_log_file = path_snp_results + "/interaction_{}_B.log".format(interaction)
	# 	A_submit_job = True
	# 	B_submit_job = True
	# 	if os.path.exists(A_plink_log_file):
	# 		plink_log_file_content = open(A_plink_log_file, 'r').read()
	# 		if (success_string_plink_1_7 in plink_log_file_content) or (success_string_plink_1_9 in plink_log_file_content):
	# 			A_submit_job = False # job is processed, do not submit again.
	# 	if os.path.exists(B_plink_log_file):
	# 		plink_log_file_content = open(B_plink_log_file, 'r').read()
	# 		if (success_string_plink_1_7 in plink_log_file_content) or (success_string_plink_1_9 in plink_log_file_content):
	# 			B_submit_job = False # job is processed, do not submit again.
	
	### CHECK IF FILE EXISTS - *** assumes that PLINK jobs was NOT KILLED while writing file ***
	### OBS: both A and B are submitted if just one of them does not exists
	snp_interactions_submit = []
	snp_interactions_no_previous_files = []
	for interaction in sorted(snp_interactions): # note that snp_interactions might not be sorted. That is why we sort it
		logger.info( "checking for previous results for interaction %s" % interaction )
		A_plink_pruned_file = path_snp_results + "/interaction_{}_A.prune.in".format(interaction) # /cvar/jhlab/timshel/egcut/interactome/5000_pruned/snp_sets/interaction_1_A.prune.in
		B_plink_pruned_file = path_snp_results + "/interaction_{}_B.prune.in".format(interaction)
		run_this_interaction = False

		if not os.path.exists(A_plink_pruned_file):
			snp_interactions_no_previous_files.append(A_plink_pruned_file)
			logger.info( "did not find paths for interaction %s: %s" % (interaction, A_plink_pruned_file) )
			run_this_interaction = True

		if not os.path.exists(B_plink_pruned_file):
			snp_interactions_no_previous_files.append(B_plink_pruned_file)
			logger.info( "did not find paths for interaction %s: %s" % (interaction, B_plink_pruned_file) )
			run_this_interaction = True

		if run_this_interaction:
			snp_interactions_submit.append(interaction)
	return (snp_interactions_submit, snp_interactions_no_previous_files)


def write_cmd_file(param, file_job):
	### CONSIDERATIONS: this could also be written as a plink script (format is different)
	""" Function to write out commands out to file."""
	# /cvar/jhlab/timshel/bin/plink1.9_linux_x86_64/plink --bfile /cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.chr_infered --extract SOMETHING --indep-pairwise 50 5 0.5
	plink_executable = "/cvar/jhlab/timshel/bin/plink1.9_linux_x86_64/plink"
	bfile = "/cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.chr_infered"
	with open(file_job, 'w') as f:
		for interaction in param:
			## extract files
			A_file_extract = path_snp_sets + "/interaction_{}_A.txt".format(interaction)
			B_file_extract = path_snp_sets + "/interaction_{}_B.txt".format(interaction)
			
			## out-prefix. NOTICE that there is NO extension.
			A_file_out = path_snp_results + "/interaction_{}_A".format(interaction)
			B_file_out = path_snp_results + "/interaction_{}_B".format(interaction)

			## commands
			A_cmd = "{} --bfile {} --extract {} --out {} --indep-pairwise 50 5 0.5 --noweb".format(plink_executable, bfile, A_file_extract, A_file_out)
			B_cmd = "{} --bfile {} --extract {} --out {} --indep-pairwise 50 5 0.5 --noweb".format(plink_executable, bfile, B_file_extract, B_file_out)

			## Write lines
			f.write(A_cmd+"\n")
			f.write(B_cmd+"\n")


def chunks_generator(l, n):
	""" Yield successive n-sized chunks from l. """
	for i in xrange(0, len(l), n):
		yield l[i:i+n]

###################################### TEMPLETE FUNCITONS ######################################

def test():
	try:
		FNULL = open(os.devnull, 'w')
		subprocess.Popen(["plink", "--silent", "--noweb"], stdout=FNULL, stderr=subprocess.STDOUT)
		FNULL.close()
	except Exception as e:
		raise Exception("Could not find plink as executable on path. Please check that you have used 'use Plink' (this is version 1.07 [USE IT!]; 'use PLINK' gives version v1.08p). Error msg: %s" % e.message)


def submit(params):
	""" 
	params: 	A list of generators/lists. 
				Each 'outer'/toplevel list contains INTERGERs for the interaction pairs that should be submitted.

	"""
	script2call = "/cvar/jhlab/timshel/git/epistasis/interactome/worker_multiprocess_v1.py"

	processes = []
	for job_no, param in enumerate(params, start=1):
		logger.info( "RUNNING: type=%s" % job_no )
		
		### write commands to file
		file_job_basename = "{timestamp}_JOBS_{n_jobs_in_batch}_{job_no}.txt".format(timestamp=timestamp, n_jobs_in_batch=n_jobs_in_batch, job_no=job_no)
		file_job = path_jobs + "/" + file_job_basename
		file_job_bsub_out = os.path.splitext(file_job_basename)[0] + "_bsub.out" # ONLY SPECIFY BASENAME (not full path)

		write_cmd_file(param, file_job)

		cmd = "python {} --file_job {}".format(script2call, file_job)
		logger.info( "adding call to processes:\n%s" % cmd )

		jobname = "prune_existing_" + str(job_no)

		processes.append( pplaunch.LaunchBsub(cmd=cmd, queue_name=queue_name, mem=mem, cpu=cpu, jobname=jobname, projectname='epistasis', path_stdout=log_dir, file_output=file_job_bsub_out, no_output=False, email=email, email_status_notification=email_status_notification, email_report=email_report, logger=logger) ) #
		#### Changes to pplaunch:
		## BEFORE file_output=None, gives path_stdout/bsub_outfile_ID630_prune_existing_630.out| NOW: file_output=file_job_bsub_out, gives 2014_10_13-18.33.42_JOBS_912_744_bsub.out
		##
	return processes


def check_jobs(processes, logger):
	logger.info("PRINTING IDs")
	list_of_pids = []
	for p in processes:
		logger.info(p.id)
		list_of_pids.append(p.id)

	logger.info( " ".join(list_of_pids) )

	if args.multiprocess:
		logger.info( "Running report_status_multiprocess " )
		pplaunch.LaunchBsub.report_status_multiprocess(list_of_pids, logger) # MULTIPROCESS
	else:
		logger.info( "Running report_status" )
		pplaunch.LaunchBsub.report_status(list_of_pids, logger) # NO MULTIPROCESS


def ParseArguments():
	arg_parser = argparse.ArgumentParser(description="Python submission Wrapper")
	arg_parser.add_argument("--logger_lvl", help="Set level for logging", choices=['debug', 'info', 'warning', 'error'], default='info') # TODO: make the program read from STDIN via '-'
	arg_parser.add_argument("--multiprocess", help="Swtich; [default is false] if set use report_status_multiprocess. Requires interactive multiprocess session", action='store_true')
	#TODO: implement formatting option
	arg_parser.add_argument("--format", type=int, choices=[0, 1, 2, 3], help="Formatting option parsed to pplaunch", default=1)
	arg_parser.add_argument("--pause", type=int, help="Sleep time after run", default=2)
	
	args = arg_parser.parse_args()
	return args

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


###################################### Global params ######################################
#queue_name = "week" # [bhour, bweek] priority
queue_name = "hour" # [bhour, bweek] priority
# priority: This queue has a per-user limit of 10 running jobs, and a run time limit of three days.
mem="5" # gb      
cpu="12"
email='pascal.timshel@gmail.com' # [use an email address 'pascal.timshel@gmail.com' or 'False'/'None']
email_status_notification=False # [True or False]
email_report=False # # [True or False]

current_script_name = os.path.basename(__file__).replace('.py','')

###################################### ARGUMENTS ######################################
args = ParseArguments()


###################################### Paramters ######################################
n_jobs_per_bsub = 200 # ***OBS: this is the NUMBER of interactions to run per job. The number of PLINK jobs is twice this number, since there are A and B files

interaction_width = 5000 # [5000, 500000, 10000000]
path_snp_sets = "/cvar/jhlab/timshel/egcut/interactome/{}/snp_sets".format(interaction_width)

path_output_base = "/cvar/jhlab/timshel/egcut/interactome/{}_pruned".format(interaction_width)
path_snp_results = path_output_base + "/snp_sets"
path_jobs = path_output_base + "/jobs" # OBS: the worker will also write to this directory.

timestamp = datetime.datetime.now().strftime("%Y_%m_%d-%H.%M.%S") # e.g. '2014_10_13-14.08.47'

### Make paths if they do not exists
for path in [path_snp_results, path_jobs]: #KEEP ME UPDATED!
	if not os.path.exists(path):
		os.makedirs(path)

###################################### SETUP logging ######################################
current_script_name = os.path.basename(__file__).replace('.py','')
log_dir = path_output_base + "/logs".format(interaction_width) #OBS VARIABLE!!!! e.g "/cvar/jhlab/timshel/egcut/interactome/5000_pruned/logs
if not os.path.exists(log_dir):
	os.mkdir(log_dir)
log_name = timestamp + '-' + current_script_name # e.g. 2014_10_13-14.08.47-wrapper_bsub_prune_sets.py. Uses timestamp!
logger = pplogger.Logger(name=log_name, log_dir=log_dir, log_format=1, enabled=True).get()
def handleException(excType, excValue, traceback, logger=logger):
	logger.error("Logging an uncaught exception", exc_info=(excType, excValue, traceback))
#### TURN THIS ON OR OFF: must correspond to enabled='True'/'False'
sys.excepthook = handleException
logger.info( "INSTANTIATION NOTE: placeholder" )
###########################################################################################

###################################### Run initial functions ######################################
LogArguments()
test() # test that things are ok
###################################### Look up input files ######################################

snp_files = glob.glob(path_snp_sets+"/*") # files will look like: /cvar/jhlab/timshel/egcut/interactome/500000/snp_sets/interaction_1_A.txt
# ^^the glob.glob scrables the order. files are NOT sorted...
snp_files_basename = [os.path.basename(x) for x in snp_files] # e.g interaction_1_A.txt, interaction_1_B.txt, ...
snp_interactions = {int(re.search(r"(\d+)", x).group(0)) for x in snp_files_basename} # SET comprehension! Entities in set is INTEGERs!

###################################### Check files ######################################
# CHECK that FILES ARE PAIRED!
#validate_interaction_file_pairs()  # EDIT

### log info
logger.info( "len(snp_interactions): %s" % len(snp_interactions) )
logger.info( "max(snp_interactions): %s" % max(snp_interactions) )
logger.info( "min(snp_interactions): %s" % min(snp_interactions) )

###################################### Remove files already processed ######################################
logger.info( "looking for previous calculated files..." )
#snp_interactions_submit = find_completed_results(snp_interactions) # returns LIST
(snp_interactions_submit, snp_interactions_no_previous_files) = find_completed_results(snp_interactions) # returns LIST

snp_interactions_chunks = list(chunks_generator(snp_interactions_submit, n_jobs_per_bsub)) # MUST use list() and not []

################## Number of jobs in this submit batch ##################
n_jobs_in_batch = len(snp_interactions_chunks)

###################################### RUN FUNCTIONS ######################################
processes = submit(snp_interactions_chunks)
####  ** INSERT THE BELOW INTO wrapper_bsub_gen_and_prune_interaction_SNP_sets.py ***

#logger.info( "interaction pairs (jobs) to submit:\n%s" % ( "\n".join(map(str, sorted(snp_interactions_submit))) ) )
logger.info( "interaction pairs (jobs) to submit:\n%s" % ( "\n".join(snp_interactions_no_previous_files) ) )
logger.info( "number of interaction pairs (jobs) to submit: %s" % len(snp_interactions_submit) )
logger.info( "IMPORTANT: number of jobs in this batch: %s (n_jobs_in_batch)" % n_jobs_in_batch )
logger.info( "IMPORTANT: n_jobs_per_bsub=%s" % n_jobs_per_bsub )


if True:
	ans = ""
	print "*** SAFETY CHECK ***"
	print "You are about to submit %s jobs to the Broad Cluster" % n_jobs_in_batch #(%s jobs to run; %s jobs per bsub)
	print "Plese confirm that this is want you want to do by typing 'yes'"
	while ans != 'yes':
	 	ans = raw_input("Confirm: ")
	print "Ok let's start..."

###################################### RUN PROCESSES ######################################
for p in processes:
	p.run()
	time.sleep(args.pause)


### Finishing
start_time_check_jobs = time.time()
check_jobs(processes, logger) # TODO: parse multiprocess argument?
elapsed_time = time.time() - start_time_check_jobs
logger.info( "Total Runtime for check_jobs: %s s (%s min)" % (elapsed_time, elapsed_time/60) )
logger.critical( "%s: finished" % current_script_name)





###### TODO #######
## make a promt to ask before you submit the jobs
## remember to use multiple CPUs

	# file_interactions = "/cvar/jhlab/timshel/egcut/interactome/lift_findItersection.intersection.paste.clean.nosex.updatedIDs"
	# with open(file_interactions, 'r') as f:
	# 	for line_no, line in enumerate(f, start=1):

	# 		line = line.strip()
	# 		fields = line.split()
	# 		(chr_A, pos_A, chr_B, pos_B, interaction_ID) = (fields[0], fields[1], fields[4], fields[5], fields[8])

	# 		try:
	# 			chr_A = int(chr_A.lstrip('chr'))
	# 			chr_B = int(chr_B.lstrip('chr'))
	# 			pos_A = int(pos_A)
	# 			pos_B = int(pos_B)
	# 		except Exception as e: # KeyboardInterrupt and SystemExit
	# 			print "line_no %s | warning: could not convert chr OR pos for A OR B." % line_no
	# 			print "line: %s" % line
	# 			print "exception instance: %s" % e
	# 			print "will log this and continue..."
	# 			f_error.write("exception:\t{}\nline_no={} | line:\t{}\n".format(e,line_no,line))
	# 			continue



