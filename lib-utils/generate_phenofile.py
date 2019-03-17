#!/usr/bin/env python2.7

import pandas as pd
import numpy as np

import sys
import glob
import os

import datetime
import time
import argparse

################## SYNOPSIS ##################
# This script will split a pheno_matrix file into VALID phenofiles for FastEpistasis/Plink.
# Each generated phenofile will contain three columns: FID, IDD and the "Phenotypename" (e.g. ILMN_XXXXX)

################## DESIGN ##################
# This script uses pandas.
# The script reads the WHOLE phenotype matrix into MEMORY

################## USAGE ##################
## Plain - write out all phenotypes from matrix
# python generate_phenofile.py --pheno_matrix egcut.peer_residuals_log2_k50.top50_mean_top50_var_refseq.txt --path_out phenofile_log2_k50.top50_mean_top50_var_refseq
## Using file subset
# python generate_phenofile.py --pheno_matrix egcut.peer_residuals_log2_k50.top50_mean_top50_var_refseq.txt --path_out phenofile_log2_k50.top50_mean_top50_var_refseq_hemani_probes_unique_102 --file_probe_subset hemani_SNPpair_probes_in_egcut.unique246.txt
## Using different pheno_matrix (HEMANI)
# python generate_phenofile.py --pheno_matrix egcut.peer_residuals_log2_k50.hemani_probes235.txt --path_out phenofile_log2_k50.hemani_probes235

################## INPUT SNIPPETs ##################
### FILE: pheno_matrix
# FID IID ILMN_2055271 ILMN_2321282 ILMN_1653355 ILMN_1717783 ILMN_1755321 ILMN_1698554 ILMN_1814092 ILMN_2061446 
# 93623 93623 -0.0393590927124023 -0.0848455429077148 -0.132477760314941 -0.0088958740234375 -0.0230574607849121 
# 27045 27045 -0.0958895683288574 0.0140762329101562 -0.0893988609313965 -0.0726613998413086 0.00134754180908203 
# 9339 9339 -0.0722055435180664 0.0204095840454102 0.0674304962158203 0.0830326080322266 -0.105680465698242 
# 83740 83740 0.00312042236328125 -0.00959682464599609 -0.238336563110352 -0.0579462051391602 -0.00715970993041992 


arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--pheno_matrix", required=True, help="SPACE seperated matrix file with header. First two columns must be FID and IID")
arg_parser.add_argument("--path_out", required=True, help="Specify a path to output files. *Path will be created if it does not exist*. May be relative or absolute from where the script is run")
arg_parser.add_argument("--file_probe_subset", help="File") # defaults to "None" is argument is not specified
args = arg_parser.parse_args()

pheno_matrix = args.pheno_matrix
path_out = os.path.abspath(args.path_out) # removes any trailing /
file_probe_subset = args.file_probe_subset

### Creating path_out if it does not exists ###
if not os.path.exists(path_out):
	print "path_out did not exist - creating new dir: %s" % path_out
	os.makedirs(path_out)

#pheno_matrix = "egcut.peer_residuals_log2_k50.top50_mean_top50_var_refseq.txt"

print "Reading pheno_matrix"
df = pd.read_csv(pheno_matrix, sep=' ', index_col=['FID', 'IID'])

if file_probe_subset is not None: # if the argument is 
	print "file_probe_subset argument was given: %s" % file_probe_subset
	with open(file_probe_subset, 'r') as f:
		probes = f.read().splitlines() # remove newlines since keep=False by default.
		print "read %s probes to subset" % len(probes)
		probes_set = set(probes)
		print "read unique %s probes to subset" % len(probes_set)

		columns_set = set(df.columns)
		match_intersection = columns_set & probes_set
		print "found #%s/#%s matching user input probes in pheno_matrix out" % (len(match_intersection), len(probes_set))
		print "subsetting data frame..."
		# IMPORTANT: updating data frame
		df = df.ix[:,match_intersection] # list(match_intersection) surely also works
		print "data frame now has %s columns" % len(df.columns)


#for column_name, series in df.ix[:,0:10].iteritems(): # ** FOR TESTING **
for i, (column_name, series) in enumerate(df.iteritems(), start=1): # Iterator over (column, series) pairs
	print "#%s/#%s: %s" % (i,len(df.columns), column_name)
	#print column_name
	filename = "%s/%s.pheno" % (path_out, column_name)
	series.to_csv(filename, sep=' ', header=True)
	# http://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.to_csv.html
	#index_label note needed
	#header IS needed for Series
	# ALTERNATIVE SOLUTION: loop over columnnames (df.columns) and use .to_csv(columns="<CURRENT_COLUMN>")












