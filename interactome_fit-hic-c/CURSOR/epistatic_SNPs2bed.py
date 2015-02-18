#!/usr/bin/env python

import os
import sys
import time

import matplotlib
matplotlib.use('Agg') #Agg backend and not an X-using backend that required an X11 connection. Call use BEFORE importing pyplot!
# REF: http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pyplot as plt

import argparse

import re
import collections
import pandas as pd
import numpy as np # for plotting an array

import glob


import pdb


###################################### USAGE ######################################

#python XXXX.py --path_main_input XXX


################## Broad ##################



################## OSX ##################

###################################### Implementation description ######################################
### Steps
# 1) Read .bim as Pandas DataFrame
# 2) Initialize epistatic_dict = defaultdict(int): default value will be zero.
	# KEYS: SNPs/rsID
	# VALUE: Maximal P-value (may be changed later)
# 3) Read file_fastepistasis_lm_combined line-by-line
# 4) MAIN LOOP: Add epistatic SNPs to dict on the fly. Save the p-value if it is the largest value seen so far.
# 5) When the file have been read, loop over all the keys (SNPs) in the epistatic_dict:
	# a) Look-up the POSITION of the SNP in the .bim file. [Make sure that the CHR from the compiled file matches the CHR in the .bim file]
	# b) Print the line for the .bed file [chr, pos, pos+1, dummy=rsID, score=P-value]
# 6) OPTINAL: Print addition stats...

###################################### TODO ######################################
### Implement options/methods for saving SNPs:
# option 1) make alternative scores
	# a) -log10(P-val)
	# b) count of SNP/rsID. That is, the number of nominal significant epistasis pairs the given SNP forms
	# c) count of SNP for each probe. That is, the number of unique probes the SNP is associated with.
# option 2) ??

### Implement ways of exporting meta data (additional columns) from the retained SNP pairs

###################################### TESTs to make ######################################


###################################### FILE snippets - *INPUT* ######################################
################## file_fastepistasis_lm_combined ##################
### Example ### 
# CHR	SNP_A	CHR	SNP_B	BETA	CHISQ	PVALUE	PHENOTYPE
# 11	rs1508531	11	rs4939298	-0.13990	33.32995	7.77755E-09	ILMN_1716816
# 1	rs4650623	6	rs3798983	-0.22814	33.00383	9.19773E-09	ILMN_1716816
# 1	rs4650623	6	rs3798982	-0.22814	33.00383	9.19773E-09	ILMN_1716816
# 4	rs12501969	21	rs762173	-0.04258	32.96917	9.36321E-09	ILMN_1666935
### Expected format ### 
# THIS IS A TAB SEPERATED FILE

################## file_bim ##################
### Example ###
# 1       rs11497407      0       558390  0       2
# 1       rs12565286      0       711153  0       2
# 1       rs11804171      0       713682  0       2
# 1       rs2977670       0       713754  0       1

###################################### FILE snippets - *OUTPUT* ######################################
################## file_epistatic_snps_bed ##################
### Example ###
# PASTE EXAMPLE HERE
### Expected format ### 
#.bed file
# See more at [HOWEVER THE CURSOR FORMAT IS *LESS* RESTRICTIVE THAN WHAT THE UCSC STATES]: http://genome.ucsc.edu/FAQ/FAQformat.html#format13
# COL#1: chrom - Name of the chromosome (or contig, scaffold, etc.).
# COL#2: chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
# COL#3: chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. If all scores were '0' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
# COL#4: *DUMMY* | name - Name given to a region (preferably unique). Use '.' if no name is assigned.
# COL#5: SCORE | Should not be negative. Will be autoscaled for 98th quantile.
# COL#6-N_meta_cols: any additional columns for meta data.

###################################### PARAMETERs ###################################### 
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--file_fastepistasis_lm_combined", required=True, help="See source code for file description")
arg_parser.add_argument("--file_output", help="Output filename (preferably completely specified: path+filename). DEFAULT filename is enabled", default="epistatic_SNPs_default_name.bed")
args = arg_parser.parse_args()

################## Set input filenames ##################
file_fastepistasis_lm_combined = os.path.abspath(args.path_main_input)

### OSX
#file_bim = "/Users/pascaltimshel/Dropbox/EGCUT_DATA/geno/Prote_370k_251011.no_mixup.chr_infered.bim"
file_bim = "/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/all_clean/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.bim" # OBS: not maf5 filtered. However, this does not matter since we only look-up the SNPs positions

### Broad
#file_bim = "/cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.maf5.bim" 


###################################### READ bim file ######################################
header_bim = ["chr", "snp", "cm", "pos", "dummy1", "dummy2"]
print "Reading .bim file..."
df_bim = pd.read_csv(file_bim, sep="\t", header=None, names=header_bim)
print "Done!"


###################################### INITIAL Error checks ######################################

################## Check that files exists ##################
if not os.path.exists(file_fastepistasis_lm_combined):
	raise Exception("file_fastepistasis_lm_combined does not exists: %s" % file_SNP2interaction_map)


###################################### Main Loop ######################################

### Initialyze main data structure
epistatic_dict = collections.defaultdict(int)

### Count the number of lines in the file
n_lines_fastepistasis_lm_combined = sum(1 for line in open(file_fastepistasis_lm_combined)) # this is pretty fast (only a bit slower than "wc -l"). 3.5 sec for a 1 GB file

time_start_loop = time.time()

### READ file line-by-line
with open(file_fastepistasis_lm_combined, 'r') as fh_compiled:
	next(fh_compiled) # SKIPPING THE FIRST LINE!
	for line_no, line in enumerate(fh_compiled, start=1):
		count_total += 1

		time_elapsed_loop = time.time() - time_start_loop # <type 'float'>
		if line_no % 5000 == 0: 
			print "Main loop | #{line_no}/#{n_lines} | {pct_complete:.2f} % done | {sec:.2f} sec [{min:.2f} min]".format(line_no=line_no, n_lines=n_lines_fastepistasis_lm_combined, pct_complete=(line_no/float(n_lines_fastepistasis_lm_combined))*100, sec=time_elapsed_loop, min=time_elapsed_loop/60)
		


		#chr1, snp_A, chr2, snp_B, beta, chisq, pvalue, phenotype
		# CHR	SNP_A	CHR	SNP_B	BETA	CHISQ	PVALUE	PHENOTYPE
		# 11	rs1508531	11	rs4939298	-0.13990	33.32995	7.77755E-09	ILMN_1716816
		line = line.strip() # strip() without argument will remove BOTH leading and trailing whitespace.
		fields = line.split() 
		chr_A, snp_A, chr_B, snp_B, pvalue = fields[0], fields[1], fields[2], fields[3], float(fields[6])

		set_A = SNP2interaction_dict[snp_A] # TODO: put this into a try/except or use dict.get() for a default value
		set_B = SNP2interaction_dict[snp_B]

		set_AB = set_A & set_B # intersection
		### REMEMBER: the set_AB contains strings of the following format:
			# experiment_interaction_identifier = {experiment_type}_{experiment_no}_{interaction_no}
		# e.g. 			
			# null_256_22233
			# null_1000_26312
			# ....

		### Mapping SNP-pairs to *UNASSIGNED* group
		if len(set_AB)==0:
			# assign to unassigned group
			count_unassigned += 1
			fh_unassigned.write(line + "\n")

		##############################################################################
		### Checking for *ASSINGED* INTRA chromosomal interactions
		if (chr_A == chr_B) and (len(set_AB) > 0):
			count_intrachromosomal_interactions += 1
			with open(file_epistatic_intrachromosomal, 'a') as fh: # OBS: append mode!
				fh.write( "{}\t{}\t{}\n".format(line, len(set_AB), ";".join(set_AB)) )
			continue # *IMPORTANT*: we skip *ASSINGED* *INTRACHROMOSOMAL* interactions
		##############################################################################

		### Mapping SNP-pairs to their experiments
		if len(set_AB)==1:
			# assign to one group
			count_assigned += 1
			experiment_interaction_identifier = iter(set_AB).next() # e.g. null_256_22233
			tmp_parts_list = experiment_interaction_identifier.split('_')
			experiment_identifier = tmp_parts_list[0] + "_" + tmp_parts_list[1] # e.g. null_256

			if pvalue <= bonferroni_correction_dict[experiment_identifier]:
				flag_significant = True
				epistatic_counts_dict['count_significant'][experiment_identifier] += 1
				epistatic_counts_dict['count_significant_pruned'][experiment_identifier].add(experiment_interaction_identifier)
			else:
				flag_significant = False

			epistatic_counts_dict['count_all'][experiment_identifier] += 1

			#experiment_interaction_identifier = list(set_AB)[0]
			experiment_identifier_fh_dict[experiment_identifier].write( "{}\t{}\t{}\n".format(line, experiment_interaction_identifier, flag_significant) )
			
		elif len(set_AB)>1:
			# "warn"
			# assigned to multiple groups
			# VALIDATION: check that hic_1 is *NOT* co-occurring with ANY null experiments
			count_assigned += 1 # OBS: SNP pairs with multiple assignments will only be counted once!
			#print "Found SNP pair with multiple assignments (#{})".format(len(set_AB))
			flag_null_seen = False
			for experiment_interaction_identifier in set_AB:
				tmp_parts_list = experiment_interaction_identifier.split('_')
				experiment_identifier = tmp_parts_list[0] + "_" + tmp_parts_list[1] # e.g. null_256

				if experiment_identifier.startswith("null"): # hmmm, there is a slight performance loss when doing this testing this 'if statement' many time. However, this will not affect the speed of the code seriously.
					flag_null_seen = True

				if (experiment_identifier == "hic_1") and (flag_null_seen): # checking for co-occurring of hic with ANY null experiments
					print "List of experiment_interaction_identifiers in set_AB:\n{}".format( "\n".join(set_AB) )
					print "DAMN, got hic_1 for a SNP-pair with 'multiple_assigments' that co-occurs with a null experiment. This means that the *NULL CONTAINS AT TRUE POSITIVE*."
					count_null_false_positives += 1
					with open(file_null_false_negatives, 'a') as fh: # OBS: append mode!
						fh.write( "{}\t{}\t{}\n".format(line, len(set_AB), ";".join(set_AB)) )

					# The reason for this 'false negatives' in the null, is that some interactions may be overlaping because of the 'width' of the regions.
					# Remember we use the same 'interationA' for the HiC and Null. Thus if a Null interaction for 'interactionB' falls suficiently close to the 'interactionB' from the HiC, they will have an overlap and the SNPs in this overlap will be in both a Null and HiC experiment.
					# Using a large width or many interaction will cause this 'false negative' problem to be bigger.
					# Remember that the R-script that contructs the null tables *ADVANCED* in the sense that it uses *chr:pos information* to ensure that *NO NULL [in interactionB] OVERLAP WITH THE HIC*.
					# *CONCLUSION*: this is not really a big problem if it happens rarely. We just need to know how often this occur to monitor the 'false negative' rate in the Null

					### IMPORTANT - recent Questions:
					# How is it possible to have two HiC experiments listed in the set_AB for two SNPs? (Should the code check be changed?)
					# How come two SNPs are listed as interchromosomal interactions?
						# --> the null must *NOT HAVE ANY INTRA-CHROMOSOMAL INTERACTIONS*.


					#raise Exception("See above print messages...")

				if pvalue <= bonferroni_correction_dict[experiment_identifier]:
					flag_significant = True
					epistatic_counts_dict['count_significant'][experiment_identifier] += 1
					epistatic_counts_dict['count_significant_pruned'][experiment_identifier].add(experiment_interaction_identifier)
				else:
					flag_significant = False

				epistatic_counts_dict['count_all'][experiment_identifier] += 1

				experiment_identifier_fh_dict[experiment_identifier].write( "{}\t{}\t{}\n".format(line, experiment_interaction_identifier, flag_significant) )
				
					
			### multiple assignment file
			fh_multiple_assignments.write(line + "\t" + ";".join(set_AB) + "\n")

print "Done with main loop!"


###################################### Writing files ######################################

with open(file_epistatic_counts, 'w') as fh:
	#for experiment_identifier in sorted(epistatic_counts_dict, key=epistatic_counts_dict.get, reverse=True):
	for experiment_identifier in bonferroni_correction_dict:
		fh.write( "{}\t{}\t{}\t{}\t{}\n".format(experiment_identifier, 
												epistatic_counts_dict['count_all'][experiment_identifier], 
												epistatic_counts_dict['count_significant'][experiment_identifier], 
												len(epistatic_counts_dict['count_significant_pruned'][experiment_identifier]), 
												bonferroni_correction_dict[experiment_identifier]),
												) # epistatic_counts_dict + epistatic_counts_dict['count_significant']
		#fh.write( "{}\t{}\n".format(experiment_identifier, epistatic_counts_dict['count_all'][experiment_identifier]) ) # 


###################################### CREATING DataFrame: Calculating Empirical P-value ######################################

#df_experiment_identifier_counts = pd.DataFrame()
df_experiment_identifier_counts = pd.DataFrame(columns=epistatic_counts_dict.keys(), index=bonferroni_correction_dict.keys()) # initialize data frame with zeroes.
df_experiment_identifier_counts = df_experiment_identifier_counts.fillna(0) # fill the data with 0s rather than NaNs
### COLUMNS IN DATAFRAME: "keys of epistatic_counts_dict"
	# count_all
	# count_significant
	# count_significant_pruned


for key_master in epistatic_counts_dict:
	for experiment_identifier in epistatic_counts_dict[key_master]:
	#for experiment_identifier in bonferroni_correction_dict:
		#df_experiment_identifier_counts[key_master][]
		if key_master == 'count_significant_pruned':
			df_experiment_identifier_counts.ix[experiment_identifier, key_master] = len(epistatic_counts_dict[key_master][experiment_identifier])
		else:
			df_experiment_identifier_counts.ix[experiment_identifier, key_master] = epistatic_counts_dict[key_master][experiment_identifier]
		#print key_master, experiment_identifier
		#if 'c' in df_experiment_identifier_counts.index: pdb.set_trace()

################## Creating P-value dict ##################
### About p-values: "obtaining a result EQUAL TO or MORE EXTREME than what was actually observed" --> p_val = sum(X >= X_OBS)
p_value_dict = {} # KEYS will be: count_significant_pruned, count_significant, count_all
for key_master in epistatic_counts_dict:
	p_value_dict[key_master] = sum(df_experiment_identifier_counts.ix[:, key_master] >= df_experiment_identifier_counts.ix['hic_1', key_master])/float(len(df_experiment_identifier_counts))
### OLD:
# p_value_count_significant_pruned = sum(df_experiment_identifier_counts.ix[:, 'count_significant_pruned'] >= df_experiment_identifier_counts.ix['hic_1', 'count_significant_pruned'])/float(len(df_experiment_identifier_counts))
# p_value_count_significant = sum(df_experiment_identifier_counts.ix[:, 'count_significant'] >= df_experiment_identifier_counts.ix['hic_1', 'count_significant'])/float(len(df_experiment_identifier_counts))
# p_value_count_all = sum(df_experiment_identifier_counts.ix[:, 'count_all'] >= df_experiment_identifier_counts.ix['hic_1', 'count_all'])/float(len(df_experiment_identifier_counts))

################## Write csv to file ##################
### Sorting, inplace
df_experiment_identifier_counts.sort(['count_significant_pruned', 'count_significant'], ascending=False, inplace=True)
### Writing file
df_experiment_identifier_counts.to_csv(file_epistatic_counts_csv) # sep='\t', index=True, header=True

################## Write stats file ##################

with open(file_epistatic_stats, 'w') as fh:
	fh.write( "count_unassigned: {}\n".format(count_unassigned) )
	fh.write( "count_assigned: {} ({:.2f} %)\n".format(count_assigned, count_assigned/float(count_unassigned+count_assigned)) )
	fh.write( "count_null_false_positives: {}\n".format(count_null_false_positives) )
	fh.write( "count_intrachromosomal_interactions: {} ({:.2f} %)\n".format(count_intrachromosomal_interactions, count_intrachromosomal_interactions/float(count_total)) )

	fh.write( "HIC count_significant_pruned: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant_pruned']) )
	fh.write( "HIC count_significant: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_significant']) )
	fh.write( "HIC count_all: {}\n".format(df_experiment_identifier_counts.ix['hic_1', 'count_all']) )

	fh.write( "p_value_dict[count_significant_pruned]: {}\n".format(p_value_dict['count_significant_pruned']) )
	fh.write( "p_value_dict[count_significant]: {}\n".format(p_value_dict['count_significant']) )
	fh.write( "p_value_dict[count_all]: {}\n".format(p_value_dict['count_all']) )

	### OLD
	# fh.write( "p_value_count_significant_pruned: {}\n".format(p_value_count_significant_pruned) )
	# fh.write( "p_value_count_significant: {}\n".format(p_value_count_significant) )
	# fh.write( "p_value_count_all: {}\n".format(p_value_count_all) )

###################################### Closing file handle ######################################

print "Closing filehandles..."
fh_multiple_assignments.close()
fh_unassigned.close()

#for experiment_interaction_identifier in experiment_identifier_fh_dict:
for fh in experiment_identifier_fh_dict.values():
	fh.close()


###################################### Plotting ######################################

################## Histogram ##################
print "will plot histogram now..."
for column_name, series in df_experiment_identifier_counts.iteritems(): # Iterator over (column, series) pairs
	print "PLOTTING HISTOGRAM FOR COLUMN_NAME: {}".format(column_name)
	fig_filename = "{base}/fig_{plot}.{ext}".format(base=path_main_out, plot=column_name, ext="pdf")

	x = series
	plt.figure()
	ax = plt.gca()  # NEEDED for "transform=ax.transAxes"
	## OR --> fig, ax = plt.subplots(1)
	## OR --> fig = plt.figure(); ax = fig.add_subplot(111)
	_ = plt.hist(x, bins=500) # 100 looks ok # ---> The return value is a tuple (n, bins, patches) 

	### Adding VERTICAL LINE for OBSERVED statistic (HiC)
	xv1 = df_experiment_identifier_counts.ix['hic_1', column_name]
	plt.axvline(xv1, color='blue', linestyle='dashed', linewidth=2)  # SEE --> http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.axhline
	
	### Adding p-value text
	textstr = "p-value={p_value}".format(p_value=p_value_dict[column_name])
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5) # # these are matplotlib.patch.Patch properties
	# place a text box in upper left in axes coords
	plt.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props) # transform=ax.transAxes
	### Our Favorite Recipes - Placing text boxes
	# http://matplotlib.org/users/recipes.html#placing-text-boxes
	# http://matplotlib.org/users/transforms_tutorial.html#transforms-tutorial
	# http://matplotlib.org/api/text_api.html#matplotlib.text.Text

	### Adding titles and labels
	plt.title(column_name, fontsize=18)
	plt.xlabel(column_name)
	plt.ylabel('Counts')
	plt.savefig(fig_filename)
	plt.close() # or plt.close(fig)


################## Subplot histogram ##################
subplot_n_row = len(df_experiment_identifier_counts.columns)
print "subplot_n_row={}".format(subplot_n_row)
for subplot_no, (column_name, series) in enumerate(df_experiment_identifier_counts.iteritems(), start=0): # Iterator over (column, series) pairs
	print "PLOTTING SUBPLOT HISTOGRAM FOR COLUMN_NAME: {}".format(column_name)
	print "subplot_no={}".format(subplot_no)

	x = series
	

	ax = plt.subplot(subplot_n_row, 1, subplot_no) # .subplot(nrow, ncol, plot_idx)
	### matplotlib.pyplot.subplot(*args, **kwargs)
	# Return a subplot axes positioned by the given grid definition.
	# Typical call signature:
	# subplot(nrows, ncols, plot_number)


	### OR ###
	# fig, ax = plt.subplots(10, 10)
	# where ax will contain one hundred axis in a list (of lists).
	#Definition: plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, **fig_kw)
	#Create a figure with a set of subplots already made.

	_ = plt.hist(x, bins=500) # 100 looks ok # ---> The return value is a tuple (n, bins, patches) 
	### CONSIDER: ax.hist()

	### Adding VERTICAL LINE for OBSERVED statistic (HiC)
	xv1 = df_experiment_identifier_counts.ix['hic_1', column_name]
	plt.axvline(xv1, color='blue', linestyle='dashed', linewidth=2)  # SEE --> http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.axhline


	### Adding p-value text
	textstr = "p-value={p_value:.3f}".format(p_value=p_value_dict[column_name])
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5) # # these are matplotlib.patch.Patch properties
	# place a text box in upper left in axes coords
	plt.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', horizontalalignment='right', bbox=props) # transform=ax.transAxes
	# plt.text(x, x, string)
	### Our Favorite Recipes - Placing text boxes
	# http://matplotlib.org/users/recipes.html#placing-text-boxes
	# http://matplotlib.org/users/transforms_tutorial.html#transforms-tutorial
	# http://matplotlib.org/api/text_api.html#matplotlib.text.Text

	### Adding axis labels and title
	plt.xlabel(column_name)
	plt.ylabel('Counts')
	plt.title(column_name, fontsize=18)

plt.tight_layout()
#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
## http://matplotlib.org/users/tight_layout_guide.html
## http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.tight_layout
## http://matplotlib.org/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels

fig_filename = "{base}/fig_{plot}.{ext}".format(base=path_main_out, plot="subplot", ext="pdf")
### saving fig
plt.savefig(fig_filename)
plt.close() # or plt.close(fig)


print "The script is complete!"

###################################### GARBAGE ######################################

	# *Main* path where the FastEpistasis results lives.
	# Expected folder organization in path_main_input:
	# <path_main_input>/
	# 	errors/ 
	# 	LSF_jobs_ids/								[optional]
	# 	cat_epistasis_YYMMDD_HHMMSS/				[optional]
	# 	config.json									[optional]
	# 	epi/										[optional]
	# 	fastEpi_compiled/							[*required*]
	# 		XXX___probe_epistatic_count.txt
	# 		XXX___results.epi.qt.lm.combined.txt	[*required*]
	# 		XXX___stat_file.txt
	# 		SNP2interaction_map/
	# 			SNP2interaction_map.pickle 			[*required*]
	# 			MORE_FILES...
	# 	jobs/										[optional]
	# 	link_bim									[optional]
	# 	link_probes									[optional]
	# 	link_set									[optional]
	# 	logs/										[optional]
	# 	logs_worker/								[optional]

	# The program relies on the *specific names* in fastEpi_compiled to match the specified pattern.

###################################### ABOUT MATPLOTLIB Axes ######################################
### Ref: http://stackoverflow.com/a/11786526

# There is a helpful page here, which provides an overview of the classes in matplotlib.

# Essentially, the process is:

# Create a figure which can hold Axes instances (and other artists)
# Create a canvas for the figure to draw itself to
# Create an Axes instance, ax, to which plotted lines/patches etc can be added. e.g. ax.plot(range(10)) or ax.contourf(array)
# I think your confusion comes from the understanding of what an Axes is. It is "the rectangular area which holds the basic elements" (for rectilinear plots). By default there is only one Axes in a figure, no matter how many times you run the command plt.plot(range(10)), although you may decide to use plt.subplot to have sub-plots in your figure, in which case you would have many Axes instances in your figure.

# HTH,


###################################### ?? ######################################

# def example_plot(ax, fontsize=12):
#      ax.plot([1, 2])
#      ax.locator_params(nbins=3)
#      ax.set_xlabel('x-label', fontsize=fontsize)
#      ax.set_ylabel('y-label', fontsize=fontsize)
#      ax.set_title('Title', fontsize=fontsize)

# plt.close('all')
# fig, ax = plt.subplots()
# example_plot(ax, fontsize=24)

