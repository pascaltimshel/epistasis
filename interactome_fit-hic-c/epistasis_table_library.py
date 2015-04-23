#!/usr/bin/env python

import re
import collections


###################################### CONSTANTS ######################################

N_PHENOTYPES_TESTED = 9269
ALPHA = 0.05

###################################### Functions ######################################

def get_blacklisted_snps():
	""" 
	Function to define and retrieve a set of blacklisted SNPs. SNPs can be blacklisted if they are e.g. duplicates from the .bim file.

	### BLACKLIST LOG ###
	04/23/2015: Added 15 unique SNPs from running "get_duplicates.py" on ~/Dropbox/5_Data/EGCUT_DATA/geno/all_clean/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.duplicates_unique_rsID.txt
	
	"""

	blacklisted_snps_set = set([
						"rs4935071",
						"rs1072594",
						"rs206146",
						"rs7383287",
						"rs2581",
						"rs9460309",
						"rs241448",
						"rs2523608",
						"rs11244",
						"rs2289150",
						"rs2070121",
						"rs1042337",
						"rs2249255",
						"rs241447"
						]) # len() --> 15



	return blacklisted_snps_set

def makehash():
	""" Function to create a perl-like hash. That is, a dict with autovivification """
	return collections.defaultdict(makehash)

def function_sort_EID(EID):
	""" 
	Function to sort EID (experiment_identifier) for dicts. This is a UTILITY FUNCTION for sorted(a_dict, key=function_sort_EID)
	EID type: string.
	EID format example(s): hic_1, null_1, null_2

	NOTICE that hic_1 and null_1 will get equal ranking because both return 1.
	"""
	experiment = EID.split("_")[0] # "hic" or "null"
	
	### We set a numeric "code" for the experiments, so "hic" has HIGHER PRECEDENCE for sorting.
	if experiment == "hic":
		experiment_ranking = 1
	elif experiment == "null":
		experiment_ranking = 2
	else:
		raise Exception( "Received unexpected experiment: {}. The function only expects experiments 'hic' or 'null'".format(experiment) )
	### 
	experiment_no = int(EID.split("_")[1]) # e.g. "1", "15" or "293" [range is 1...n_perm].
	
	tuple_to_sort_by = (experiment, experiment_no)
	# SEE: Python secondary sorting: http://stackoverflow.com/a/16193637
	return tuple_to_sort_by

def function_sort_EIID(EIID):
	"""
	Function to sort EIID.
	EIID type: string.
	EIID format example(s): hic_1_2013, hic_1_2574, null_2_1074
	
	See function_sort_EID() for additional documentation.
	"""
	interaction_no = int(EIID.split("_")[-1]) # e.g. "2013" [range is 1..n_interactions]
	return interaction_no
