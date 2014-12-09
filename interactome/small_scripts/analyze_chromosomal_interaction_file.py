#!/usr/bin/env python

import os
import sys
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg') #Agg backend and not an X-using backend that required an X11 connection. Call use BEFORE importing pyplot!
# REF: http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pyplot as plt

import datetime
import time
import argparse

#####
import collections
import operator


### file_interactions
### IMPORTANT information: 
#1) the first column listing the position will be used (e.g. 67280788 and not 67280789)
#2) the left part of the file (column 1-4) will be denoted interaction_XX_A. 
#3) the right part of the file (column 5-8) wil be denoted interaction_XX_B.
#4) THE last column contains "updated" interaction IDs. There are 92932 of these updated (there where 96137 of the original)
# chr12   67280788        67280789        interaction1    chr12   66561751        66561752        interaction1    interaction_1
# chr3    150566213       150566214       interaction2    chr3    150061015       150061016       interaction2    interaction_2
# chr10   28973034        28973035        interaction3    chr10   28917902        28917903        interaction3    interaction_3
# chr10   28978176        28978177        interaction4    chr10   29702662        29702663        interaction4    interaction_4
# chr17   52725358        52725359        interaction5    chr17   53563372        53563373        interaction5    interaction_5


file_interactions = "/cvar/jhlab/timshel/egcut/interactome/lift_findItersection.intersection.paste.clean.nosex.updatedIDs"

file_out = "analyze_chromosomal_interaction_file.out"
#sys.stdout = open(file_out, 'w')

pair_count_dict = collections.Counter()
#pair_count_dict = collections.defaultdict(int)

pair_trans_count_dict = collections.Counter()

###################################### Processing ######################################
with open(file_interactions, 'r') as f:
	for line_no, line in enumerate(f, start=1):
		line = line.strip()
		fields = line.split()
		(chr_A, pos_A, chr_B, pos_B, interaction_ID) = (fields[0], fields[1], fields[4], fields[5], fields[8])

		chr_A = chr_A.lstrip('chr')
		chr_B = chr_B.lstrip('chr')
		pos_A = pos_A
		pos_B = pos_B

		# Number of pairs of interactions: binom(k=2, n=22) --> 231 # 22 autosomes
		interaction_pair = sorted([chr_A, chr_B]) # IMPORTANT TO SORT!
		interaction_pair_string = ":".join(interaction_pair)
		pair_count_dict.update([interaction_pair_string]) # must be encapsulated in list
		## Counting trans pairs
		if chr_A != chr_B:
			pair_trans_count_dict.update([interaction_pair_string]) # must be encapsulated in list


print "Found %s number of chromosome interaction pairs. Expected around 231" % len(pair_count_dict)

for k,v in sorted(pair_count_dict.items(), key=operator.itemgetter(1)):
	print "%s=%s" % (k,v)

pair_count_series = pd.Series(pair_count_dict)
print pair_count_series.describe()

print "sum of interactions: %s" % pair_count_series.sum()


################## Trans ##################
print "TRANS INTERACTIONS!"
for k in sorted(pair_trans_count_dict, key=pair_trans_count_dict.get):
	print "%s=%s" % (k,pair_trans_count_dict[k])


pair_trans_count_series = pd.Series(pair_trans_count_dict)
print pair_trans_count_series.describe()

print "sum of interactions: %s" % pair_trans_count_series.sum()


##################  ##################
# print pair_count_dict["13:9"]
# print pair_trans_count_dict["13:9"]


################## Plotting ##################
# for column_name, series in df_set_stats.iteritems(): # Iterator over (column, series) pairs
# 	print column_name
# 	file_figure = "{}/fig_{}.pdf".format(path_figs, column_name)

# 	plt.figure()
# 	_ = plt.hist(series.dropna()) # needed for NaN values in distance
# 	plt.title(column_name, fontsize=18)
# 	plt.xlabel(column_name)
# 	plt.ylabel('counts')
# 	#plt.savefig(file_figure)
# 	plt.savefig(file_figure)
# 	plt.close() # or plt.close(fig)


