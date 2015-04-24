#!/usr/bin/env python

import os
import sys
import time

###################################### Library ######################################

import argparse
import subprocess

###################################### USAGE ######################################

#python epistasis_table_run.py --epi_job_identifier XXX

### Remarks ###
# This script only works on the Broad file system

###################################### SCRIPT ######################################

# reopen stdout file descriptor with write mode
# and 0 as the buffer size (unbuffered)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

###################################### PARAMETERs ###################################### 
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--epi_job_identifier", required=True, help="""
	e.g.
	hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8
	hESC_width_1000_maf_5_q_1e-12_epi1_1e-10
	""")
args = arg_parser.parse_args()

epi_job_identifier = args.epi_job_identifier

###################################### Run subprocesses ######################################

# subprocess.check_call(args, *, stdin=None, stdout=None, stderr=None, shell=False)
# Run command with arguments. Wait for command to complete. If the return code was zero then return, otherwise raise CalledProcessError. The CalledProcessError object will have the return code in the returncode attribute.


################## epistasis_table_labeller ##################
path_main_input_epistasis_table_labeller = "/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/{job}/fastEpi_compiled".format(job=epi_job_identifier)
cmd = "python epistasis_table_labeller.py --path_main_input {}".format(path_main_input_epistasis_table_labeller)
print "running cmd={}".format(cmd)
p = subprocess.check_call(cmd, shell=True) # Notice that the std/err goes to this scripts stream (stdout)
print "DONE!"

################## epistasis_table_processing ##################


path_main_input_epistasis_table_processing = "/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/{job}/fastEpi_compiled/assigned".format(job=epi_job_identifier)
cmd = "python epistasis_table_processing.py --path_main_input {}".format(path_main_input_epistasis_table_processing)
print "running cmd={}".format(cmd)
p = subprocess.check_call(cmd, shell=True) # Notice that the std/err goes to this scripts stream (stdout)
print "DONE!"


###################################### FINISH ######################################
print "The script is done"