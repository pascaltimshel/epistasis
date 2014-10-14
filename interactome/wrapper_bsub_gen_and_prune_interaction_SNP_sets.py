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


###################################### file_interactions ######################################
### file_interactions
### IMPORTANT information: 
#1) the first column listing the position will be used (e.g. 67280788 and not 67280789)
#2) the left part of the file (column 1-4) will be denoted interaction_XX_A. 
#3) the right part of the file (column 5-8) wil be denoted interaction_XX_B.
#4) THE last column contains "updated" interaction IDs.
# chr12   67280788        67280789        interaction1    chr12   66561751        66561752        interaction1    interaction_1
# chr3    150566213       150566214       interaction2    chr3    150061015       150061016       interaction2    interaction_2
# chr10   28973034        28973035        interaction3    chr10   28917902        28917903        interaction3    interaction_3
# chr10   28978176        28978177        interaction4    chr10   29702662        29702663        interaction4    interaction_4
# chr17   52725358        52725359        interaction5    chr17   53563372        53563373        interaction5    interaction_5


def read_file_interactions():
	file_interactions = "/cvar/jhlab/timshel/egcut/interactome/lift_findItersection.intersection.paste.clean.nosex.updatedIDs"
	interaction_ID_dict = {}
	with open(file_interactions, 'r') as f:
		for line_no, line in enumerate(f, start=1):

			line = line.strip()
			fields = line.split()
			(chr_A, pos_A, chr_B, pos_B, interaction_ID) = (fields[0], fields[1], fields[4], fields[5], fields[8])

			try:
				interaction = int(interaction_ID.lstrip('interaction_')) # OBS: important
				chr_A = int(chr_A.lstrip('chr'))
				chr_B = int(chr_B.lstrip('chr'))
				pos_A = int(pos_A)
				pos_B = int(pos_B)
			except Exception as e: # KeyboardInterrupt and SystemExit
				logger.critical( "line_no %s | warning: could not convert chr OR pos for A OR B." % line_no )
				logger.critical( "line: %s" % line )
				logger.critical( "exception instance: %s" % e )
				logger.critical( "will re-raise execption..." )
				raise
			
			interaction_ID_dict[interaction] = {'chr_A':chr_A, 'pos_A':pos_A, 'chr_B':chr_B, 'pos_B':pos_B}
	return interaction_ID_dict


# ###################################### NEW FUNCITONS ######################################
######## MAY BE DELETED!! #####
# def find_completed_results(interaction_ID_dict):
# 	### CHECK IF FILE EXISTS - *** assumes that PLINK jobs was NOT KILLED while writing file ***
# 	### OBS: both A and B are submitted if just one of them does not exists
# 	snp_interactions_submit = []
# 	for interaction in sorted(interaction_ID_dict.keys()): # note that interaction_ID_dict might not be sorted
# 		logger.info( "checking for previous results for interaction %s" % interaction )
# 		A_plink_pruned_file = path_snp_results + "/interaction_{}_A.prune.in".format(interaction) # /cvar/jhlab/timshel/egcut/interactome/5000_pruned/snp_sets/interaction_1_A.prune.in
# 		B_plink_pruned_file = path_snp_results + "/interaction_{}_B.prune.in".format(interaction)
# 		if not ( os.path.exists(A_plink_pruned_file) and os.path.exists(B_plink_pruned_file) ):
# 			#interaction_ID_dict_submit[interaction] = interaction_ID_dict[interaction] # COPY dict
# 			snp_interactions_submit.append(interaction)
# 	return snp_interactions_submit

###################################### NEW FUNCITONS ######################################
def find_completed_results(interaction_ID_dict):
	### CHECK IF FILE EXISTS - *** assumes that PLINK jobs was NOT KILLED while writing file ***
	### OBS: both A and B are submitted if just one of them does not exists
	logger.info( "looking for files in %s" % path_snp_results )
	ttmp = time.time()
	#snp_files = glob.glob(path_snp_results+"/*") # the glob.glob scrables the order. files are NOT sorted...
	snp_files = os.listdir(path_snp_results) # the list is in arbitrary order
	logger.info( "done globbing. time: %s" % (time.time()-ttmp,) )
	#snp_files_basename = [os.path.basename(x) for x in snp_files] # e.g interaction_1_A.txt, interaction_1_B.txt, ...
	#snp_interactions = {int(re.search(r"(\d+)", x).group(0)) for x in snp_files_basename} # SET comprehension! Entities in set is INTEGERs!

	snp_interactions_submit = []
	for interaction in sorted(interaction_ID_dict.keys()): # note that interaction_ID_dict might not be sorted
		logger.info( "checking for previous results for interaction %s" % interaction )
		A_plink_pruned_file = path_snp_results + "/interaction_{}_A.prune.in".format(interaction) # /cvar/jhlab/timshel/egcut/interactome/5000_pruned/snp_sets/interaction_1_A.prune.in
		B_plink_pruned_file = path_snp_results + "/interaction_{}_B.prune.in".format(interaction)
		if not ( (A_plink_pruned_file in snp_files) and (B_plink_pruned_file in snp_files) ):
			#interaction_ID_dict_submit[interaction] = interaction_ID_dict[interaction] # COPY dict
			snp_interactions_submit.append(interaction)
	return snp_interactions_submit


def write_cmd_file(param, file_job):
	### CONSIDERATIONS: this could also be written as a plink script (format is different)
	""" Function to write out commands out to file."""
	# /cvar/jhlab/timshel/bin/plink1.9_linux_x86_64/plink --bfile /cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.chr_infered --extract SOMETHING --indep-pairwise 50 5 0.5
	plink_executable = "/cvar/jhlab/timshel/bin/plink1.9_linux_x86_64/plink"
	bfile = "/cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.chr_infered"
	with open(file_job, 'w') as f:
		for interaction in param:
			## get chr and pos from dict
			chr_A = interaction_ID_dict[interaction]['chr_A'] # int
			pos_A = interaction_ID_dict[interaction]['pos_A'] # int, in bp
			pos_A_lower = pos_A - interaction_width
			pos_A_upper = pos_A + interaction_width


			chr_B = interaction_ID_dict[interaction]['chr_B']
			pos_B = interaction_ID_dict[interaction]['pos_B']
			pos_B_lower = pos_B - interaction_width
			pos_B_upper = pos_B + interaction_width
			
			## out-prefix. NOTICE that there is NO extension.
			A_file_out = path_snp_results + "/interaction_{}_A".format(interaction)
			B_file_out = path_snp_results + "/interaction_{}_B".format(interaction)

			### r^2 param: SNP pairs with "r^2 > threshold" are greedily pruned.
			## commands
			A_cmd = "{exe} --bfile {bfile} --chr {chr} --from-bp {from_bp} --to-bp {to_bp} --out {out} --indep-pairwise 50 5 0.8 --noweb".format(exe=plink_executable, bfile=bfile, chr=chr_A, from_bp=pos_A_lower, to_bp=pos_A_upper, out=A_file_out)
			B_cmd = "{exe} --bfile {bfile} --chr {chr} --from-bp {from_bp} --to-bp {to_bp} --out {out} --indep-pairwise 50 5 0.8 --noweb".format(exe=plink_executable, bfile=bfile, chr=chr_B, from_bp=pos_B_lower, to_bp=pos_B_upper, out=B_file_out)

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

		jobname = "genANDprune_" + str(job_no)

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
n_jobs_per_bsub = 150 # ***OBS: this is the NUMBER of interactions to run per job. The number of PLINK jobs is twice this number, since there are A and B files

interaction_width = 5000
path_output_base = "/cvar/jhlab/timshel/egcut/interactome/width_{}_pruned".format(interaction_width)
path_snp_results = path_output_base + "/snp_sets"
path_jobs = path_output_base + "/jobs" # OBS: the worker will also write to this directory.

timestamp = datetime.datetime.now().strftime("%Y_%m_%d-%H.%M.%S") # e.g. '2014_10_13-14.08.47'

### Make paths if they do not exists
for path in [path_snp_results, path_jobs]: #KEEP ME UPDATED!
	if not os.path.exists(path):
		logger.warning( "OBS - path exists: %s" % path )
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

###################################### *** Read in interactions from file *** ######################################
interaction_ID_dict = read_file_interactions() # the number of keys in the dict corresponds to the number of interactions

###################################### Remove files already processed ######################################
logger.info( "looking for previous calculated files..." )
snp_interactions_submit = find_completed_results(interaction_ID_dict) # returns LIST of interactions (integers)
logger.info( "number of interaction pairs (jobs) to submit: %s" % len(snp_interactions_submit) )

snp_interactions_chunks = list(chunks_generator(snp_interactions_submit, n_jobs_per_bsub)) # MUST use list() and not []

################## Number of jobs in this submit batch ##################
n_jobs_in_batch = len(snp_interactions_chunks)

if True:
	ans = ""
	print "*** SAFETY CHECK ***"
	print "You are about to submit %s jobs to the Broad Cluster" % n_jobs_in_batch #(%s jobs to run; %s jobs per bsub)
	print "Plese confirm that this is want you want to do by typing 'yes'"
	while ans != 'yes':
	 	ans = raw_input("Confirm: ")
	print "Ok let's start..."

###################################### RUN FUNCTIONS ######################################
processes = submit(snp_interactions_chunks)

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



