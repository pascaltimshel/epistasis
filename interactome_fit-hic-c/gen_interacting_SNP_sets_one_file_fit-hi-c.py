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

import collections

###################################### USAGE ######################################

#python gen_interacting_SNP_sets_one_file_fit-hi-c.py --width 500000

###################################### WARNING ######################################

# I am not sure that the "caching" functing of the script works in this script.
# I recommend to ALWAYS delete ANY "old" result directory and re-run this analysis.

###################################### FILE snippets ######################################
### file_bim
# 1       rs11497407      0       558390  0       2
# 1       rs12565286      0       711153  0       2
# 1       rs11804171      0       713682  0       2
# 1       rs2977670       0       713754  0       1


### file_interactions
### IMPORTANT information: 
#0) The file has a *HEADER*
#1) The file is sorted chr(chrA) <= chr(chrB)
#2) The chromosomes are NUMERIC/INTEGER. No sex chromosomes.
#3) THE last column contains "updated" interaction IDs. 
# chrA	posA	chrB	posB	contactCount	p.value	q.value	interactionID	interactionID_final
# 10	100515000	15	71455000	5	2.097017e-16	3.081624e-09	interaction_393	interaction_1
# 5	152585000	10	100665000	5	1.602572e-14	7.197691e-08	interaction_521	interaction_2
# 2	8165000	10	100725000	5	1.821807e-15	1.39996e-08	interaction_555	interaction_3
# 10	100955000	16	69065000	5	1.526798e-15	1.238718e-08	interaction_727	interaction_4


###################################### PARAMETERs ###################################### 
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--width", required=True, help="INTERGER. Width in bp. Width is the distance to boundary from the focal/center point. x - width <= snp_pos <= x + width") # default=10*1000
args = arg_parser.parse_args()

# in bp. try [20kb, 500kb, 2Mb and 10Mb] or [10kb, 250kb, 1Mb and 5Mb], where the numbers are TWICE the distance to the region boundary from the focal point.
interaction_width = int(args.width)
#interaction_width = 500*1000 

###################################### Files - input ######################################
### OSX
#file_bim = "/Users/pascaltimshel/Dropbox/EGCUT_DATA/geno/Prote_370k_251011.no_mixup.chr_infered.bim"
#file_interactions = "/Users/pascaltimshel/p_HiC/Lan_et_al_chromosomal_interactions/lift_findItersection.intersection.paste.updatedIDs"

##############################################################################################
### *IMPORTANT*: setting the specificiers for the data set. THIS IS USED IN path_base_out *###
#hic_data_set = "hIMR90"
hic_data_set = "hESC"

list_of_q_threshold = ["1e-12", "1e-13", "1e-14", "1e-15", "1e-16", "1e-17", "1e-18"]

### Setting threshold
#q_threshold = "1e-06" 
#q_threshold = "1e-07"
#q_threshold = "1e-08"
#q_threshold = "1e-09"
#q_threshold = "1e-10"

for q_threshold in list_of_q_threshold:
	print "running {}".format(q_threshold)

	file_interactions = "/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/interation_table.fit-hi-c.nosex.interchromosomal.{dataset}.q_{q}.txt".format(dataset=hic_data_set, q=q_threshold)

	#file_interactions = "/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-06.txt"
	#file_interactions = "/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-07.txt"
	#file_interactions = "/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-08.txt"
	#file_interactions = "/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-09.txt"
	#file_interactions = "/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-10.txt"

	#file_interactions = "/cvar/jhlab/timshel/egcut/interactome/lift_findItersection.intersection.paste.updatedIDs"
	#file_interactions = "/cvar/jhlab/timshel/egcut/interactome/lift_findItersection.intersection.paste.clean.nosex.updatedIDs.interchromosomal"

	###############################################################################################

	### Broad
	#file_bim = "/cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.bim" # REMEMBER: use the CLEAN FILE later!
	#file_bim = "/cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.bim" 
	file_bim = "/cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.maf5.bim" 

	## wc 12/09/2014
	#  2009553 = Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.bim
	#  2552493 = Prote_370k_251011.no_mixup.with_ETypes.chr_infered.bim


	###################################### SANITY CHECK ######################################
	if not q_threshold in file_interactions:
		raise Exception( "ERROR: could not find '%s' (value of q_threshold) as substring in file_interactions=%s" % ( q_threshold, file_interactions ) )

	else:
		print "OK: q_threshold is in file_interaction"


	#path_base_out = "/Users/pascaltimshel/p_HiC/Lan_et_al_interaction_SNP_sets/{}".format(interaction_width) # OSX
	path_base_out = "/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/{width}_snppool_{dataset}_q_{q}".format(width=interaction_width, dataset=hic_data_set, q=q_threshold) # BROAD

	path_snp_sets = path_base_out+"/snp_sets"
	path_figs = path_base_out+"/figs"
	path_stats = path_base_out+"/stats"
	path_errors = path_base_out+"/errors"


	paths = [path_base_out, path_snp_sets, path_figs, path_stats, path_errors]
	for path in paths:
		if not os.path.exists(path):
			os.makedirs(path)
		else:
			print "warining: path exists %s" % path

	###################################### Files ######################################
	file_stats = path_stats+"/df_set_stats.txt"
	file_summary = path_stats+"/df_set_summary.txt"
	#file_interaction_table = path_stats+"/df_interaction_table.txt"
	file_interaction_table_snps = path_stats+"/df_interaction_table_snps.txt"

	################## Error file ##################
	timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H-%M-%S") # e.g. '2014_10_11_15:04:22'
	file_error = path_errors + '/error_{}'.format(timestamp)
	f_error = open(file_error, 'w')

	###################################### Checking for previous files ######################################
	if os.path.exists(file_stats): # this check could also be for file_interaction_table/file_interaction_table_snps
		df_set_stats = pd.read_csv(file_stats, sep="\t", index_col=0)  # index_col=0 IS needed, but header=0 is NOT needed
		
		#df_interaction_table = pd.read_csv(file_interaction_table, sep="\t", index_col=0)  # index_col=0 IS needed, but header=0 is NOT needed
		df_interaction_table_snps = pd.read_csv(file_interaction_table_snps, sep="\t", index_col=0)  # index_col=0 IS needed, but header=0 is NOT needed
		print "Detected that existing df_set_stats file exists. Will start from here..."
		print "LAST INTERACTION CALCULATED: %s" % df_set_stats.index[-1] 
	else:
		df_set_stats = pd.DataFrame() # has no index
		#df_interaction_table = pd.DataFrame() # has no index
		df_interaction_table_snps = pd.DataFrame() # has no index

		print "No previous df_set_stats file exists. Will start from scratch!"


	###################################### Function for stat files ######################################

	def write_summaries():
		""" Function to (overwrite) existing summary files """
		################## Creating summary ##################
		df_set_summary = pd.Series(
			{
			'setA_size_mean': df_set_stats['setA_size'].mean(),
			'setB_size_mean': df_set_stats['setB_size'].mean(),
			'set_n_interchromosomal_sum': df_set_stats['set_interchromosomal'].sum(),
			'set_n_tests_non_redundant_sum':df_set_stats['set_n_tests_non_redundant'].sum(),
			'set_n_tests_plink_sum':df_set_stats['set_n_tests_plink'].sum(),
			'set_intersect_mean':df_set_stats['set_intersect'].mean(),
			'set_percentage_shared_mean':df_set_stats['set_percentage_shared'].mean(),
			'set_lenght_df_set_stats':len(df_set_stats), 
			'setA_n_empty_sets': sum(df_set_stats['setA_size']==0),
			'setB_n_empty_sets': sum(df_set_stats['setB_size']==0),
			'set_n_empty_sets': sum(df_set_stats['setA_size']==0)+sum(df_set_stats['setB_size']==0)
			})

		################## Write stats to csv ##################
		### Writing file including column wiht list of comma seperated SNPs for setA and setB
		df_interaction_table_snps.to_csv(file_interaction_table_snps, sep="\t", index=False) # DataFrame. This writes the index by default

		### Write the data frame with PER INTERACTION data (large file)
		df_set_stats.to_csv(file_stats, sep="\t", index=False) # DataFrame. This writes the index by default

		### Write current SUMMARY of stats. (small file)
		df_set_summary.to_csv(file_summary, sep="\t", index=True) # Series




	###################################### READ bim file ######################################
	header_bim = ["chr", "snp", "cm", "pos", "dummy1", "dummy2"]
	df_bim = pd.read_csv(file_bim, sep="\t", header=None, names=header_bim)

	time_minutes_count = 0
	time_old = 0
	### Some timing
	list_time_lookup = []
	list_time_set_operations = []

	###################################### Processing ######################################
	### setting file names
	file_setA = path_snp_sets + "/set_A.txt"
	file_setB = path_snp_sets + "/set_B.txt"
	with open(file_interactions, 'r') as f:
		for line_no, line in enumerate(f, start=1):

			### *Skipping header* ###
			if line_no==1: 
				print "skipping header!"
				continue

			## USE THIS FOR TESTING ##
			#if line_no == 100: break

			time_elapsed =  time.time() - time_old
			if time_elapsed > 60: 
				time_minutes_count += 1
				print "Time elapsed: %s minute(s)" % time_minutes_count
				time_old = time.time()

			line = line.strip()
			fields = line.split()
			# KEEP THIS SNIPPED UPDATED
			# chrA	posA	chrB	posB	contactCount	p.value	q.value	interactionID	interactionID_final
			# 10	100515000	15	71455000	5	2.097017e-16	3.081624e-09	interaction_393	interaction_1
			(chr_A, pos_A, chr_B, pos_B, interaction_ID) = (fields[0], fields[1], fields[2], fields[3], fields[8])
			
			
			### IMPORTANT: "BUFFERING". skipping processing interaction if it is already in df_set_stats.index
			### This buffering (i.e. writing summaries every 100 lines/interactions) ensures that at max 99 set files are overwritten (and avoids processing of interactions twice).
			## Note that the buffering relies on interaction_ID NOT to change format. E.g. interaction1 and interaction_1 are not the same!
			if interaction_ID in df_set_stats.index:
				print "%s: already seen in df_set_stats. skipping it..." % interaction_ID
				continue

			### Convert to chr and pos to INTEGERs
			try:
				chr_A = int(chr_A)
				chr_B = int(chr_B)
				pos_A = int(pos_A)
				pos_B = int(pos_B)
			except Exception as e: # KeyboardInterrupt and SystemExit
				print "line_no %s | warning: could not convert chr OR pos for A OR B." % line_no
				print "line: %s" % line
				print "exception instance: %s" % e
				print "will log this and continue..."
				f_error.write("exception:\t{}\nline_no={} | line:\t{}\n".format(e,line_no,line))
				continue
			# Check that chr_A and chr_B are as we except
			#if not (isinstance(chr_A, int) and isinstance(chr_B, int)): # or try is_integer()


			if line_no % 50 == 0: 
				print "line_no: {} |".format(line_no), chr_A, pos_A, chr_B, pos_B, interaction_ID
			
			###################################### LOOK-UP SNPs ######################################
			tmp_time_lookup_start = time.time()
			### A
			df_A_extract = df_bim[df_bim["chr"]==chr_A]
			df_A_extract = df_A_extract[(df_A_extract["pos"] >= pos_A-interaction_width) & (df_A_extract["pos"] <= pos_A+interaction_width)]
			### B
			df_B_extract = df_bim[df_bim["chr"]==chr_B]
			df_B_extract = df_B_extract[(df_B_extract["pos"] >= pos_B-interaction_width) & (df_B_extract["pos"] <= pos_B+interaction_width)]       
			
			list_time_lookup.append(time.time() - tmp_time_lookup_start) # this is a float
			if line_no % 50 == 0:
				print "avg look-up time: %s s" % ( sum(list_time_lookup)/float(len(list_time_lookup)) )
				list_time_lookup = []


			########################################################################################

			###################################### SET COMPUTATIONS ######################################
			tmp_time_set_operations_start = time.time()

			## set calculations
			tmp_intersection = len(set(df_A_extract['snp']).intersection(set(df_B_extract['snp'])))
			tmp_union = len(set(df_A_extract['snp']).union(set(df_B_extract['snp'])))
			## number of tests
			n_tests_plink = len(df_A_extract)*len(df_B_extract)-tmp_intersection
			n_tests_non_redundant = len(df_A_extract)*len(df_B_extract)-tmp_intersection**2+(tmp_intersection-1)*tmp_intersection/2 # int division ok! # N_A*N_B-N_AB^2+(N_AB-1)*N_AB/2
			
			### tmp_df_stats
			tmp_df_stats = pd.DataFrame({
									'interaction_ID':interaction_ID,
									'chr_A':chr_A,
									'pos_A':pos_A,
									'chr_B':chr_B,
									'pos_B':pos_A,
									'set_interchromosomal':1 if chr_A != chr_B else 0,
									'set_distance': abs(pos_A - pos_B) if chr_A == chr_B else None, # or np.nan
									'set_intersect':tmp_intersection,
									'set_union':tmp_union,
									'set_self_interaction':1 if (chr_A == chr_B) and (pos_A == pos_B) else 0,
									'set_percentage_shared':tmp_intersection/float(tmp_union)*100 if tmp_union!=0 else None, # float divison
									'set_n_tests_plink':n_tests_plink,
									'set_n_tests_non_redundant': n_tests_non_redundant,
									'setA_size':len(df_A_extract), 
									'setA_range':df_A_extract["pos"].max()-df_A_extract["pos"].min(), 
									'setB_size':len(df_B_extract), 
									'setB_range':df_B_extract["pos"].max()-df_B_extract["pos"].min()}, 
									index=[interaction_ID])
			df_set_stats = df_set_stats.append(tmp_df_stats) # index and columns are appended to the end of the dataframe
			
			list_time_set_operations.append(float(time.time() - tmp_time_set_operations_start)) # this is a float
			if line_no % 50 == 0:
				print "avg set-operation time: %s s" % ( sum(list_time_set_operations)/float(len(list_time_set_operations)) )
				list_time_set_operations = []
			# If no columns are passed, the columns will be the sorted list of dict keys.
			##############################################################################


			### tmp_df_interaction_table_snps
			colorder = ["interaction_ID", "chr_A", "pos_A", "chr_B", "pos_B", "setA_size", "setB_size", "snps_A", "snps_B"]
			# *IMPORTANT* the column ORDER MUST BE SYNCHRONIZED!!!
			tmp_df_interaction_table_snps = pd.DataFrame({
				'interaction_ID':interaction_ID,
				'chr_A':chr_A,
				'pos_A':pos_A,
				'chr_B':chr_B,
				'pos_B':pos_A,
				'setA_size':len(df_A_extract),
				'setB_size':len(df_B_extract),
				'snps_A': ";".join(df_A_extract["snp"]),
				'snps_B': ";".join(df_B_extract["snp"])},
				index=[interaction_ID], columns=colorder)
			df_interaction_table_snps = df_interaction_table_snps.append(tmp_df_interaction_table_snps) # index and columns are appended to the end of the dataframe


			### NEW: *APPEND* to SET file. NOTE: this causes problems with the CACHING properties.
			### UPDATE: it should NOT cause a problem BECAUSE WE ARE ONLY INTERESTED IN UNIQUE SNPS --> the "set_AB.txt" will only contain unique SNPs.
			with open(file_setA, 'a') as fA:
				for snp in df_A_extract["snp"]:
					fA.write(snp+"\t"+interaction_ID+"\n")
			with open(file_setB, 'a') as fB:
				for snp in df_B_extract["snp"]:
					fB.write(snp+"\t"+interaction_ID+"\n")

			### write summaries every hundred line. this if statement will NOT be triggered at the first time, because line_no starts at 1
			if line_no % 100 == 0: 
				print "writing summary..."
				time_tmp = time.time()
				write_summaries() # updating stat files
				print "done (%s s)" % (time.time()-time_tmp, )

	###################################### Writing final file ######################################
	print "writing summary FOR the LAST time..."
	time_tmp = time.time()
	write_summaries() # IMPORTANT: must make sure that we write file when all is done
	print "done (%s s)" % (time.time()-time_tmp, )


	###################################### Plotting ######################################

	for column_name, series in df_set_stats.iteritems(): # Iterator over (column, series) pairs
		print "making plot: %s" % column_name
		file_figure = "{}/fig_{}.pdf".format(path_figs, column_name)
		series_no_NA = series.dropna() # needed for NaN values in distance
		print "number of observations for plotting: %s" % len(series_no_NA)
		try:
			plt.figure()
			_ = plt.hist(series_no_NA) 
			plt.title(column_name, fontsize=18)
			plt.xlabel(column_name)
			plt.ylabel('counts')
			#plt.savefig(file_figure)
			plt.savefig(file_figure)
			plt.close() # or plt.close(fig)
		except Exception, e:
			print "Exception occurred: %s" % e # 
			#print "Exception message: %s" % e.message # this works. gives ---> Exception message: index out of bounds
			#print "Exception args: %s" % e.args # this works. gives ---> Exception args: index out of bounds
			#print vars(e) # gives *WEIRD* ---> {}
			#print(dir(e)) # dir will give a list of names comprising the methods and attributes of an object
				# gives --> ['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__getitem__', '__getslice__', '__hash__', '__init__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__', '__unicode__', 'args', 'message']
			# There is one more public attribute, message, but if you access that you will get a DeprecationWarning. Deprecated attributes are not always documented since you shouldn't use them in new code.




	###################################### Closing files ######################################
	f_error.close()


	###################################### *NEW*: Combining set_A and set_B ######################################
	### setting file names
	file_setA_duplicates = path_snp_sets + "/set_A_duplicates.txt"
	file_setB_duplicates = path_snp_sets + "/set_B_duplicates.txt"
	file_setAB_stats = path_snp_sets + "/set_AB_stats.tab"
	file_setAB = path_snp_sets + "/set_AB.txt"

	with open(file_setA, 'r') as fA, open(file_setB, 'r') as fB:
		################## SetA ##################
		setA_snps = [line.split("\t")[0] for line in fA.readlines()] # no stripping because we take first field
		#setA_snps_duplicates = [key for key, count in collections.Counter(setA_snps).items() if count > 1]
		print "making collections.Counter() for setA and writing to file"
		n_setA_duplicated_snps = 0
		with open(file_setA_duplicates, 'w') as f:
			for key, count in collections.Counter(setA_snps).most_common(): # could also use .items() for unsorted
				if count > 1: # duplicted SNP
					n_setA_duplicated_snps += 1
					f.write("%s\t%s\n" % (key,count) )

		################## SetB ##################
		setB_snps = [line.split("\t")[0] for line in fB.readlines()] # no stripping because we take first field
		print "making collections.Counter() for setB and writing to file"
		n_setB_duplicated_snps = 0
		with open(file_setB_duplicates, 'w') as f:
			for key, count in collections.Counter(setB_snps).most_common(): # could also use .items() for unsorted
				if count > 1: # duplicted SNP
					n_setB_duplicated_snps += 1
					f.write("%s\t%s\n" % (key,count) )

		################## Computing and Writing stats ##################
		print "Computing and Writing stats... "
		setA_snps_unique = set(setA_snps)
		n_setA_unique_snps = len(setA_snps_unique)
		print "number of unique setA: %s" % n_setA_unique_snps
		setB_snps_unique = set(setB_snps)
		n_setB_unique_snps = len(setB_snps_unique)
		print "number of unique setB: %s" % n_setB_unique_snps
		n_set_AB_intersection = len(setA_snps_unique & setB_snps_unique) # same as setA_snps_unique.intersection(setB_snps_unique)
		n_set_AB_union = len(setA_snps_unique | setB_snps_unique) # same as setA_snps_unique.intersection(setB_snps_unique)
		
		with open(file_setAB_stats, 'w') as f:
			f.write("n_setA_unique_snps\t%s\n" % n_setA_unique_snps)
			f.write("n_setB_unique_snps\t%s\n" % n_setB_unique_snps)
			f.write("n_setA_duplicated_snps\t%s\n" % n_setA_duplicated_snps)
			f.write("n_setB_duplicated_snps\t%s\n" % n_setB_duplicated_snps)
			f.write("n_set_AB_intersection\t%s\n" % n_set_AB_intersection)
			f.write("n_set_AB_union\t%s\n" % n_set_AB_union)
			f.write("RATIO=n_set_AB_intersection/n_set_AB_union*100\t%s%%\n" % ( round(n_set_AB_intersection/float(n_set_AB_union)*100, 3) ,) )

		################## Writing final list ##################
		#Recognized keywords are SET_A, SET_B and END
		# SET_A
		# rs17763185
		# rs1367975
		# ...
		# END
		# SET_B
		# rs17763185
		# rs12544008
		# ...
		# END
		print "Writing out file_setAB: %s" % file_setAB
		with open(file_setAB, "w") as f:
			### SetA
			f.write("SET_A\n")
			for snp in setA_snps_unique:
				f.write(snp+"\n")
			f.write("END\n")
			
			### SetB
			f.write("SET_B\n")
			for snp in setB_snps_unique:
				f.write(snp+"\n")
			f.write("END\n")

	print "The script is complete!"

