#!/usr/bin/env python

import os
import sys

import collections

import memory_profiler; import psutil # OBS: the program RUN EXTREMELY SLOW when memory_profiling WITHOUT psutil. To install pip in homedir on a server: "pip install -user psutil"

import time

try:
	import cPickle as pickle # NO subclass - cPickle is written in C.
except:
	print "loaded pickle - could not load cPickle..."
	import pickle

@memory_profiler.profile
def load_snp2interaction_dict():
	with open(file_in, 'r') as fh:
		SNP2interaction_dict = pickle.load(fh) # pickle.load(file)
	return SNP2interaction_dict

file_in = "/cvar/jhlab/timshel/egcut/fastEpistasis_fit-hi-ci/hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10/fastEpi_compiled/SNP2interaction_map/SNP2interaction_map.pickle"

SNP2interaction_dict = load_snp2interaction_dict()
print "Loaded load_snp2interaction_dict"

print "Will make test print..."
for key, elem in SNP2interaction_dict.iteritems():
    print "{}:{}".format(key, elem)
    break

#dict.itervalues().next()
#dict[dict.keys()[0]]

print "sleeping for 5 seconds"
time.sleep(5)

print "The script is complete!"

