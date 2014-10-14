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

import subprocess
import logging

import itertools
import multiprocessing


def read_jobs(file_job):
	with open(file_job, 'r') as f:
		jobs = f.read().splitlines()
		return jobs


def run(job):
	### TODO: consider making a try except clause here! See http://stackoverflow.com/questions/16450788/python-running-subprocess-in-parallel

	### THROW the output to the bitbucket
	with open('/dev/null', 'w') as fnull:
		return_report = subprocess.call(job, stdout=fnull, stderr=subprocess.STDOUT, shell=True) # returns returncode. hmmm, starts a shell...

	##### READING FROM PIPES. NOT TESTED, but should work! ####
	# p = subprocess.Popen(job, stdout=subprocess.PIPE, subprocess.STDOUT, shell=True) # hmmm, starts a shell...
	# (stdoutdata, stderrdata) = p.communicate() # wait for process to finish and read file-object. note that both stdout and stderr is piped to stdout. Thus stderr is None
	# return_report = None
	# for line in stdoutdata.splitlines():
	# 	if ('Analysis finished' in line) or ('End time' in line):
	# 		return_report = "{};{}".format(p.returncode, line) # return only result line
	
	# if return_report is None: # 'success' line was not found.
	# 	return_report = "{};{}".format(p.returncode, stdoutdata) # returns all output as a string			
	return return_report

	
### Read arguments
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--file_job", help="file with jobs. one line per job")
args = arg_parser.parse_args()

### Output
file_job = args.file_job ## this file MUST exists. Otherwise script will crash
file_out = os.path.join(os.path.dirname(file_job), os.path.splitext(os.path.basename(file_job))[0] + "_worker.txt") # extract basename without extenstion. this works even if the file does not have an extension
file_bad_out = os.path.join(os.path.dirname(file_job), "fail_jobs_" + os.path.splitext(os.path.basename(file_job))[0] + "_worker.txt") 
sys.stdout = open(file_out, 'w') # Redirect output
#sys.stdrr = open(file_out, 'w') # Redirect errors output (so it does not end up in the bsub out)

time_start = time.time()
print "#file_out: {}".format(file_out)
print "#file_bad_out: {}".format(file_bad_out)
print "#Reading jobs to process..."
jobs = read_jobs(file_job)
print "#Number of jobs: %s" % len(jobs)
pool = multiprocessing.Pool() # use all available CPUs, processes=multiprocessing.cpu_count()
print "#Processes in pool: %s" % multiprocessing.cpu_count()
return_reports = pool.map(run, jobs)
for job, return_report in itertools.izip_longest(jobs, return_reports, fillvalue=None):
	print "job={}|return_report={}".format(job, return_report)
	if return_report != 0:
		with open(file_bad_out, 'a') as f_bad:
			f_bad.write( "job={}|return_report={}\n".format(job, return_report) )


time_elapsed = (time.time() - time_start)/60
print "#DONE: time=%s min" % time_elapsed


# >>> os.path.splitext('/path/to/somefile.ext')
# ('/path/to/somefile', '.ext')
# >>> os.path.splitext('/path/to/somefile')
# ('/path/to/somefile', '')

# def safe_run():
#     """Call run(), catch exceptions."""
#     try: 
#     	run()
#     except Exception as e:
#         print("error in run(): %s" % e)

# PLINK1.9
# End time: Sun Oct 12 17:53:38 2014

#PLINK1.7
#Analysis finished: Sun Oct 12 17:49:43 2014
