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

import json

import pdb

###################################### ABOUT the CONFIG FILE ######################################
# To ensure that the fastepistasis runs are run using the same parameters, the config file is NOT overwritten if it exists.
# If you REALLY need to make adjustments to the parameters (MAF, q_value, ..) should MUST MANUALLY MODIFY the config file.

###################################### Description ######################################
## IMPORATNT: This script has XXX important variables to check:
## 2) path_files 			for input sets. USED FOR GLOBBING
## 							**OBS: this path must ONLY contain "INPUT" files. OTHERWISE THE PROGRAM MAY FAIL!
##							consider using the following to check the folder integrity: ls -1af | grep -v '<PATTERN>.*txt'						
##							TODO: make a 'static list' of files to process instead
## 3) path_epi_results	 	for caching property
## CACHING PROPERTY
## The script searches the path_epi_results for previous results
## Thus you need to check that the expected OUTPUT files are still as described in this program. SEE find_completed_results()
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



def find_completed_results(variable_elements, file_path_format_holder):
	""" 
	Caching property
	
	INPUTS
		variable_components: a iterable (e.g. set, list)
		file_path_format_holder: a string used to give the format of the file. FULL PATH TO FILE. E.g. "BASEPATH/VARIABLE/resultXX_{}_A.out"
	"""
	results_missing = []
	results_present = []
	n_variable_elements = len(variable_elements)
	logger.info( "checking for previous results for variable_element..." )
	for i, variable_element in enumerate(sorted(variable_elements), start=1): # note that variable_elements might not be sorted. That is why we sort it
		#if i == 1000: break
		#if i == 2: break

		logger.info( "#%s/#%s: %s" % (i, n_variable_elements, variable_element) )
		
		file_result = file_path_format_holder.format(variable_element) # FORMATTING STRING
		
		##pdb.set_trace()
		
		if not os.path.exists(file_result):
			results_missing.append(variable_element)
			logger.info( "NO previous file for variable_element %s: %s" % (variable_element, file_result) )
		else:
			results_present.append(variable_element)
	return (results_missing, results_present)


def chunks_generator(l, n):
	""" Yield successive n-sized chunks from l. """
	for i in xrange(0, len(l), n):
		yield l[i:i+n]


def write_cmd_file(chunck, file_job):
	""" Function to write out commands out to file."""
	with open(file_job, 'w') as f:
		for variable_element in chunck:
			f.write(variable_element+"\n") # e.g. ILMN_1691436 (no .pheno)

###################################### TEMPLETE FUNCITONS ######################################

def test():
	""" The .fastepistasis-2.05 needs to be loaded because this script submits jobs that requires it to be in the path """
	
	try:
		FNULL = open(os.devnull, 'w')
		subprocess.Popen(["smpFastEpistasis"], stdout=FNULL, stderr=subprocess.STDOUT)
		FNULL.close()
	except Exception as e:
		raise Exception("Could not find <PROGRAM OF INTEREST> as executable on path. Please check that you have LOADED THE PROGRAM. Error msg: %s" % e)
	# reuse .fastepistasis-2.05


def submit(chuncks):
	""" 
	The jobs are added to the LaunchBsub class. The jobs are run LATER using the "p.run()" method

	chuncks: 	A LIST of GENERATORS/LISTS (Note: this a NESTED structure)
				params = [[batch1_x1, batch1_x2, ..., batch1_xN], [batch2_x1, batch2_x2, ..., batch2_xN], ..., NUMBER_OF_BATCHES]
				Each 'inner'/nested list contains *VARIABLE_ELEMENTS* for the jobs that should be submitted.

	"""

	processes = []
	for job_no, chunck in enumerate(chuncks, start=1):
		logger.info( "RUNNING: type=%s" % job_no )
		
		### 
		# file_job_basename = "{timestamp}_JOBS_{n_jobs_in_batch}_{job_no}.txt".format(timestamp=timestamp, n_jobs_in_batch=n_jobs_in_batch, job_no=job_no)
		# file_job = path_jobs + "/" + file_job_basename
		# file_job_bsub_out = os.path.splitext(file_job_basename)[0] + "_bsub.out" # ONLY SPECIFY BASENAME (not full path)
		
		file_job_basename = "{timestamp}_JOBS_{n_jobs_in_batch}_{job_no}".format(timestamp=timestamp, n_jobs_in_batch=n_jobs_in_batch, job_no=job_no)
		file_job = path_jobs + "/{base}.{ext}".format(base=file_job_basename, ext="txt")
		### write commands/variable_elements to file
		write_cmd_file(chunck, file_job)

		### Bsub output file: *NOT full path* - only basename
		file_job_bsub_out = "{base}.{ext}".format(base=file_job_basename, ext="bsub.out") # *OBS* ONLY SPECIFY BASENAME (not full path). This variable is passed to pplaunch.LaunchBsub()
		
		### APPLICATION SPECIFIC ###
		## run_fast_epistasis_inputlist.py --> file_log
		## *OBS* notice the WORKER LOG PATH "path_logs_worker"
		file_job_worker_file_log = path_logs_worker + "/{base}.{ext}".format(base=file_job_basename, ext="worker.out")


		cmd = "python {} --file_job {} --path_epi_results {} --file_log {}".format(script2call, file_job, path_epi_results, file_job_worker_file_log) # OBS: path_epi_results is NOT variable
		logger.info( "adding call to processes:\n%s" % cmd )

		#jobname = "epi_" + str(job_no)
		jobname = "epi-job_no_%s-%s" % (job_no, path_output_name) # Remember: job_no is an INT

		processes.append( pplaunch.LaunchBsub(cmd=cmd, queue_name=queue_name, mem=mem, proc=proc, shared_mem=shared_mem, jobname=jobname, projectname='epistasis', path_stdout=log_dir, file_output=file_job_bsub_out, no_output=False, email=email, email_status_notification=email_status_notification, email_report=email_report, logger=logger) ) #
		## CONSIDER COMSTUM COMMAND
		#processes.append( pplaunch.LaunchBsub(cmd=cmd, queue_name=queue_name, mem=mem, proc=proc, cmd_custom=SOMETING))
		
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
	arg_parser.add_argument("--pause", type=int, help="Sleep time after run", default=0)
	
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

###################################### TEST environment variables ######################################
test() # test that things are ok

###################################### Global params ######################################
queue_name = "week" # --> REL max_jobs_per_user=200 | CentOS max_jobs_per_user=300
#queue_name = "hour" # --> REL max_jobs_per_user=400 | CentOS max_jobs_per_user=500
#queue_name = "priority" # --> 10
#queue_name = "MEDPOP" # --> max_jobs_per_user=NA
# priority: This queue has a per-user limit of 10 running jobs, and a run time limit of three days.
mem="5" # gb      
proc=20 # number of CPUs
shared_mem=True # if True: sets -R 'span[hosts=1]' | use this for multithreaded SMP jobs
#shared_mem=False

email='pascal.timshel@gmail.com' # [use an email address 'pascal.timshel@gmail.com' or 'False'/'None']
email_status_notification=False # [True or False]
email_report=False # # [True or False]

current_script_name = os.path.basename(__file__).replace('.py','')

###################################### ARGUMENTS ######################################
args = ParseArguments()

###################################### Paramters ######################################
# ***OBS: this is the NUMBER of interactions to run per job.
n_jobs_per_bsub = 200 #500 # 60 # --> RHEL WEEK
#n_jobs_per_bsub = 25 # *NEW*
#n_jobs_per_bsub = 4 #25 # 100 # --> RHEL HOUR
#n_jobs_per_bsub = 100 # --> RHEL HOUR v. 2
#n_jobs_per_bsub = 2000 # --> MEDPOP
#n_jobs_per_bsub = 20 # --> CentOS HOUR

### IMPORTANT: define the script to call
script2call = "/cvar/jhlab/timshel/git/epistasis/interactome_fit-hic-c/run_fast_epistasis_inputlist.py"

######################################  ######################################
### Overall parameters
maf = "5"
#interaction_width = "50000" # --> RHEL hour
#interaction_width = "500" # --> RHEL MEDPOP
#interaction_width = "500" # --> RHEL hour v. 2
#interaction_width = "10000" # --> CentOS hour
#interaction_width = "10000"
#interaction_width = "2500" # --> RHEL week - LASTEST hIMR90
#interaction_width = "2500" # --> RHEL week - LASTEST hESC
interaction_width = "1000" # --> RHEL week

hic_dataset = "hIMR90"
#hic_dataset = "hESC"

#q_threshold = "1e-09" # --> RHEL hour
#q_threshold = "1e-08" # --> CentOS hour
#q_threshold = "1e-08" # --> RHEL hour v. 2
q_threshold = "1e-06" # --> RHEL week
#q_threshold = "1e-10"
#q_threshold = "1e-07" # --> RHEL week - LASTEST hIMR90
#q_threshold = "1e-13" # --> RHEL week - LASTEST hESC
#q_threshold = "1e-14" # --> RHEL week - LASTEST hESC

### Input BIM file: *UPS: keep this in sync!* ###
#bfile = "/cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean" # DO NOT ADD EXTENSION to file
file_bim = "/cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.maf5" # DO NOT ADD EXTENSION to file

path_files = "/cvar/jhlab/timshel/egcut/ETypes_probes_norm_peer/phenofile_log2_k50.top50_mean_top50_var_refseq"
#file_set = "/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_5_sets/10000_snppool_hIMR90_q_1e-07/snp_sets/set_AB.txt"

### Formatting file_set
file_set = "/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/maf_{maf}_sets/{interaction_width}_snppool_{hic_dataset}_q_{q_threshold}/snp_sets/set_AB.txt".format(maf=maf, interaction_width=interaction_width, hic_dataset=hic_dataset, q_threshold=q_threshold)
if not os.path.exists(file_set): # checks for a valid symbolic link. That is, it cannot be broken
	raise Exception( "The file_set symlink is broken - path does not exists: {}".format(file_set) )
else:
	print "OK - found the file_set file. Will create symlink soon..."

### Significance threshold: epi1 ###
epi1 = "1e-10" # string to set the significance thresshold
epi2 = "1e-10" # string to set the significance thresshold

#epi1 = "1e-8" # MEDPOP + RHEL hour v. 2
#epi2 = "1e-8" # MEDPOP + RHEL hour v. 2



################## OUTPUT dirs/files ##################
path_output_name = "{hic_dataset}_width_{interaction_width}_maf_{maf}_q_{q_threshold}_epi1_{epi1}".format(hic_dataset=hic_dataset, interaction_width=interaction_width, maf=maf, q_threshold=q_threshold, epi1=epi1)
path_output_base = "/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/{name}".format(name=path_output_name) # inserting name


path_epi_results = path_output_base + "/epi"
### *OBS: IMPORTANT*
file_path_format_holder = path_epi_results + "/{}.pheno.epi.qt.lm" # OBS: THIS STRING NEEDS TO BE COMPATIBLE WITH THE .format() string formatting. That is, the "%s" formatting will *NOT* work
	# e.g. --> PATH_EPI_RESULTS/ILMN_1751898.pheno.epi.qt.lm


path_jobs = path_output_base + "/jobs" # OBS: the worker will also write to this directory.
path_lsf_job_ids = path_output_base + "/LSF_jobs_ids"
path_logs_worker = path_output_base + "/logs_worker" # OBS: this MUST be in sync with the name in "run_fast_epistasis.py"


file_json_config = path_output_base + "/config.json" # json file!


################## MISC ##################
timestamp = datetime.datetime.now().strftime("%Y_%m_%d-%H.%M.%S") # e.g. '2014_10_13-14.08.47'
file_lsf_job_id = path_lsf_job_ids + "/" + timestamp + '-LSF_JOB_IDs.txt' # e.g. 2014_10_13-14.08.47-LSF_JOB_IDs.txt 
###################################### MAKE DIRS/SYMLINKs ######################################
### Safety check - 
if not os.path.exists(path_output_base):
	print "WARNING: you are about to create a new result directory. path_output_base does not exists: %s" % path_output_base
	ans = ""
	while ans != 'yes':
		ans = raw_input("Is this what you want to do ('yes')?: ")
	print "Ok. Will continue..."
	time.sleep(1)


### Make paths if they do not exists
for path in [path_output_base, path_epi_results, path_jobs, path_lsf_job_ids, path_logs_worker]: #KEEP ME UPDATED!
	if not os.path.exists(path):
		os.makedirs(path)

### Make symlinks
symlink_dict = {"link_bim":file_bim, "link_probes":path_files, "link_set":file_set}
for symlink_name, src_path in symlink_dict.items():
	symlink_path = path_output_base + "/" + symlink_name # e.g. /cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/PATH_OUTPUT_BASE/link_bim

	if not os.path.lexists(symlink_path): # .lexists() - CHECK FOR EXISTANCE OF SYMLINK: Return True if path refers to an existing path. Returns True for broken symbolic links
		# we NEED to make sure that there is no symlink created already. If there is an existing symlink, an error will be raised!
		# note that os.path.exists() could also be used. It just returns False if the link i broken. This is actually ok, since in case of a broken link, to code to fail anyway
		
		# os.symlink(src,dest)
		os.symlink(src_path, symlink_path)


###################################### SETUP logging ######################################
current_script_name = os.path.basename(__file__).replace('.py','')
log_dir = path_output_base + "/logs" #OBS VARIABLE!!!! 
if not os.path.exists(log_dir):
	os.mkdir(log_dir)
log_name = timestamp + '-' + current_script_name # e.g. 2014_10_13-14.08.47-wrapper_bsub_prune_sets.py. Uses timestamp!
logger = pplogger.Logger(name=log_name, log_dir=log_dir, log_format=1, enabled=True, log2stream=True, log2file=True).get()
def handleException(excType, excValue, traceback, logger=logger):
	logger.error("Logging an uncaught exception", exc_info=(excType, excValue, traceback))
#### TURN THIS ON OR OFF: must correspond to enabled='True'/'False'
sys.excepthook = handleException
logger.info( "INSTANTIATION NOTE: placeholder" )
###########################################################################################

###################################### Run initial functions ######################################
LogArguments()


###################################### WRITE CONFIG FILE ######################################

config_dict = {
		"link_bim":file_bim, 
		"link_probes":path_files, 
		"link_set":file_set,
		"epi1":epi1,
		"epi2":epi2,
		"n_jobs_per_bsub": n_jobs_per_bsub,
		"maf": maf,
		"interaction_width": interaction_width,
		"hic_dataset": hic_dataset,
		"q_threshold": q_threshold,
		"nthreads": proc,
		"queue_name": queue_name
}

# check for existance of a PREVIOUS config file
if not os.path.exists(file_json_config):
	with open(file_json_config, 'w') as f_json:
		json.dump(config_dict, f_json, sort_keys=True, indent=2)
else:
	logger.info( "CONFIG FILE EXISTS. Will read it and check for any changed configs." )
	
	### Checking for changed config
	with open(file_json_config, 'r') as f_json:
		config_dict_previous = json.load(f_json)
	changed_config = [] # this list will contain "keys"
	for key in config_dict:
		if not config_dict[key] == config_dict_previous.get(key): # .get() is needed if key is not present in old config file
			changed_config.append(key)
	#pdb.set_trace()
	if changed_config: # list is non-empty
		logger.warning( "CONFIG SETTINGS HAVE CHANGED." )
		logger.warning( "Printing the complete new config file:" )
		for key in config_dict:
			logger.warning( "{}:{}".format(key, config_dict[key]) )
		logger.warning( "Here are the changes in the config file:" )
		for key in changed_config:
			logger.info( "Config={} | New value=[{}] | Old value=[{}]".format(key, config_dict.get(key), config_dict_previous.get(key)) )
			#logger.info( "If you want your changes (new values) to take effect you must MANUALLY EDIT the config file." )

		### User input - CONFIRMATION ###
		ans = ""
		while ans != 'yes':
		 	#ans = raw_input("Are you sure you want to continue with the old values (type 'yes')?")
		 	ans = raw_input("Are you sure you want to OVERWRITE the old values (type 'yes')? ")
	 	### Renaming old config file
	 	file_json_config_old_file = "%s.before_%s" % (file_json_config, timestamp)
	 	logger.info( "Renaming old config file. FROM=%s TO=%s" % ( os.path.basename(file_json_config), os.path.basename(file_json_config_old_file) ) )
	 	os.rename(file_json_config, file_json_config_old_file) # os.rename(src, dst)
	 	### Writing new file
	 	with open(file_json_config, 'w') as f_json:
	 		json.dump(config_dict, f_json, sort_keys=True, indent=2)
	else:
		logger.info( "Config settings are identical to previous config file." )




	#config_set_previous_values = config_dict_previous.value()

#json.load(j)

###################################### Look up input files ######################################

#EXAMPLE: /cvar/jhlab/timshel/egcut/ETypes_probes_norm_peer/egcut.peer_residuals_log2_k50.top50_mean_top50_var_refseq.txt/ILMN_1691436.pheno

list_files = glob.glob(path_files+"/*") # files will look like: XXX
# ^^the glob.glob scrables the order. files are NOT sorted...
list_files_basename = [os.path.basename(x) for x in list_files] # e.g ILMN_1691436.pheno, ILMN_1691436.pheno, ...
variable_elements = {x.split(".")[0] for x in list_files_basename} # SET comprehension!
	# EXAMPLE: ILMN_1691436
#snp_interactions = {int(re.search(r"(\d+)", x).group(0)) for x in list_files_basename} # SET comprehension! Entities in set is INTEGERs!

###################################### Check files ######################################
# CHECK that FILES ARE PAIRED!
#validate_interaction_file_pairs()  # EDIT

### log info
logger.info( "len(variable_elements): %s" % len(variable_elements) )
#logger.info( "max(variable_elements): %s" % max(variable_elements) )
#logger.info( "min(variable_elements): %s" % min(variable_elements) )

###################################### Remove files already processed ######################################
logger.info( "looking for previous calculated files..." )
(results_missing, results_present) = find_completed_results(variable_elements, file_path_format_holder) # returns two LISTs
	# IMPORTANT: ---> the function "find_completed_results()" returns "variable_elements" ("atomic" names)
#pdb.set_trace()
n_result_missing = len(results_missing)
n_result_present = len(results_present)

results_missing_chunks = list(chunks_generator(results_missing, n_jobs_per_bsub)) # MUST use list() and not []
#pdb.set_trace()

################## Number of jobs in this submit batch ##################
n_jobs_in_batch = len(results_missing_chunks)

###################################### RUN FUNCTIONS ######################################
processes = submit(results_missing_chunks)

logger.info( "List of jobs to submit (results_missing):\n%s" % ( "\n".join(results_missing) ) ) # This will likely produce a lot of output. Keep it "ABOVE" the other print statements so it does not fill your terminal with output
logger.info( "Number of total jobs to submit - sum over all batches (results_missing): %s" % n_result_missing )
logger.info( "Number of PREVIOUS COMPLETED jobs: %s out of %s total jobs" % ( n_result_present,  n_result_present+n_result_missing) )
logger.info( "IMPORTANT: n_jobs_in_batch=%s" % n_jobs_in_batch )
logger.info( "IMPORTANT: n_jobs_per_bsub=%s" % n_jobs_per_bsub )


if True:
	ans = ""
	print "*** SAFETY CHECK ***"
	print "You are about to submit %s jobs (%s in total) to the Broad Cluster (queue_name=%s)" % ( n_jobs_in_batch, n_result_missing, queue_name ) #(%s jobs to run; %s jobs per bsub)
	print "Plese confirm that this is want you want to do by typing 'yes'"
	while ans != 'yes':
	 	ans = raw_input("Confirm: ")
	print "Ok let's start..."

###################################### RUN PROCESSES ######################################
#lsf_job_id = []
with open(file_lsf_job_id, 'a') as f_lsf:
	for p in processes:
		p.run() # after the .run() is succesfull, the p.id attribute is set.
		f_lsf.write(str(p.id)+" ") # space seperated
		time.sleep(args.pause)
#bjobs | egrep "week\s+node1747" | perl -ane 'print "$F[0] "'

### Finishing
start_time_check_jobs = time.time()
check_jobs(processes, logger) # TODO: parse multiprocess argument?
elapsed_time = time.time() - start_time_check_jobs
logger.info( "Total Runtime for check_jobs: %s s (%s min)" % (elapsed_time, elapsed_time/60) )
logger.critical( "%s: finished" % current_script_name)



