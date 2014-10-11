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


import pdb

def chunks_generator(l, n):
	""" Yield successive n-sized chunks from l. """
	for i in xrange(0, len(l), n):
		yield l[i:i+n]


interaction_width = 500000
path_snp_sets = "/cvar/jhlab/timshel/egcut/interactome/{}/snp_sets/".format(interaction_width)

snp_files = glob.glob(path_snp_sets+"/*") # files will look like: /cvar/jhlab/timshel/egcut/interactome/500000/snp_sets/interaction1_A.txt
snp_files_basename = [os.path.basename(x) for x in snp_files]
snp_interactions = [int(re.search(r"(\d+)", x).group(0)) for x in snp_files_basename]
snp_interactions.sort() # the glob.glob scrables the order. Sorting is good
print len(snp_interactions)
print max(snp_interactions)
print min(snp_interactions)


snp_interactions_chunks = list(chunks_generator(snp_interactions, 100)) # MUST use list()
for x in snp_interactions_chunks:
	print x
#print snp_interactions_chunks

#snp_interactions_chunks = chunks_list(snp_interactions, 100)
#print snp_interactions_chunks
