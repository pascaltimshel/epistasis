#!/usr/bin/env python

import os
import sys
import time

###################################### Library ######################################

import argparse
import subprocess

###################################### USAGE ######################################

#python epistasis_table_processing_run.py

###################################### SCRIPT ######################################

###################################### PARAMETERs ###################################### 

epi_job_identifier_list = [
"hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10",
"hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8",
"hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8",
"hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8",
"hIMR90_width_1000_maf_5_q_1e-06_epi1_1e-10",

###### hESC [4] #######
"hESC_width_1000_maf_5_q_1e-12_epi1_1e-10",
"hESC_width_500_maf_5_q_1e-14_epi1_1e-8",
"hESC_width_500_maf_5_q_1e-16_epi1_1e-8",
"hESC_width_2500_maf_5_q_1e-13_epi1_1e-10",

###### lan-et-al_K562 [2] #######
"lan-et-al_K562_width_1000_maf_5_q_OUTLIER_RM_epi1_1e-8",
"lan-et-al_K562_width_5000_maf_5_q_OUTLIER_RM_epi1_1e-8",

###### contactCount_1 #######
"hESC-contactCount_1_width_1000_maf_5_q_1_epi1_1e-8",
"hIMR90-contactCount_1_width_1000_maf_5_q_1_epi1_1e-8", # trailing comma "," is ok
]


#epi_job_identifier_list= ["hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8"]
#epi_job_identifier_list= ["hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8"]

###################################### Run subprocesses ######################################

# subprocess.check_call(args, *, stdin=None, stdout=None, stderr=None, shell=False)
# Run command with arguments. Wait for command to complete. If the return code was zero then return, otherwise raise CalledProcessError. The CalledProcessError object will have the return code in the returncode attribute.

################## epistasis_table_processing ##################

for epi_job_identifier in epi_job_identifier_list:
	path_main_input_epistasis_table_processing = "/Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/{job}".format(job=epi_job_identifier)
	cmd = "python epistasis_table_processing.py --path_main_input {}".format(path_main_input_epistasis_table_processing)
	print "running cmd={}".format(cmd)
	p = subprocess.check_call(cmd, shell=True) # Notice that the std/err goes to this scripts stream (stdout)
	print "DONE: {}".format(epi_job_identifier)


###################################### FINISH ######################################
print "The script is done"