#!/usr/bin/env python2.7

import sys
import os
import argparse
import collections

import pdb

import datetime

###################################### IMPORTANT ######################################
# Date created: April 2015
# This script is COPIED and MODIFIED from the SNPsnap library

### Modifications
# 1) added writing out the chromosome to "appendfile"
# 2) changed name for "appendfile"
# 3) *IMPORTANT*: now duplicate chr:pos is also detected and written to file

###################################### USAGE ######################################

### Broad
#python get_duplicates.py --input /cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.maf5.bim

### OSX
#python get_duplicates.py --input /Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/all_clean/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.bim

###################################### SCRIPT ######################################

def get_duplicates(inputfile, outputfile, appendfile):
	### rsID
	snps_seen_rsID = {}
	dup_rsID = collections.defaultdict(list)
	
	### chrpos
	snps_seen_chrpos = {}
	dup_chrpos = collections.defaultdict(list)
	
	print "will read inputfile..."
	
	with open(outputfile, 'w') as outfile:
		with open(inputfile, 'r') as infile:
		
		# .bim file == 'extended' MAP file - each line of the MAP file describes a single marker
		
	 	# c1=chromosome (1-22, X, Y or 0 if unplaced)
	    # c2=rs# or snp identifier
	    # c3=Genetic distance (morgans)
	    # c4=Base-pair position (bp units)
	    # c5=allele names (extra col in .bim file comared to MAP)
	    # c6=allele names (extra col in .bim file comared to MAP)
		# 1       rs202029170     0       247916  C       CAGG
		# 1       rs200079338     0       249275  GT      G
		# 1       rs115018998     0       249276  C       T
		# 1       rs72502741      0       251627  A       AC
		# 1       rs199745078     0       255923  GTC     G
		# 1       rs182870673     0       362905  G       T

			lines = infile.readlines()
			for line in lines:
				cols = line.strip().split()
				chromosome = cols[0]
				rsID = cols[1]
				pos = cols[3]

				chrpos = "{}:{}".format(chromosome, pos) # e.g. 3:124322

				### rsID
				if not rsID in snps_seen_rsID:
					snps_seen_rsID[rsID] = chrpos
				else:
					# pdb.set_trace()
					#if len(snps_seen_rsID[rsID]) == 1: # if first time we see duplicate
					if len(dup_rsID[rsID]) == 0: # if first time we see duplicate
						dup_rsID[rsID].extend([snps_seen_rsID[rsID], chrpos]) # save previous seen rs position and duplicate position
					else:
						dup_rsID[rsID].append(chrpos)

				### chrpos
				if not chrpos in snps_seen_chrpos:
					snps_seen_chrpos[chrpos] = rsID
				else:
					if len(dup_chrpos[chrpos]) == 0: # if first time we see duplicate
						dup_chrpos[chrpos].extend([snps_seen_chrpos[chrpos], rsID]) # save previous seen rs position and duplicate position
					else:
						dup_chrpos[chrpos].append(rsID)
		
		print "done read inputfile! will start writing to outfile"
		
		######################### Generating and write unique set of duplicated SNPs ########################
		duplicated_snps = set()
		
		######### populate set #########
		### rsID
		for rs in dup_rsID:
			duplicated_snps.add(rs)

		### chrpos
		for chrpos, rs_list in dup_chrpos.items():
			for rs in rs_list:
				duplicated_snps.add(rs)

		######### write to file #########
		for rs in duplicated_snps:
			outfile.write(rs+"\n")

		###################################### Append file ######################################
		print "will write to appendfile"
		with open(appendfile, 'a') as fappend:
			now_string = datetime.datetime.now().strftime("%a %b %d %Y %H:%M:%S") # e.g. Wed Sep 24 2014 12:51:54'
			fappend.write( "#"*80 + '\n' )
			fappend.write( "INPUT FILE: %s\n" % inputfile )
			fappend.write( "TIME: %s\n" % now_string )
			result_string = "Duplicates Detected!" if dup_rsID else "No Duplicates"
			fappend.write( "RESULT: %s\n" % result_string )
			### rsID
			if dup_rsID:
				fappend.write("### DUPLICATED rsID ###\n")
				fappend.write("rsID\tCount\tPositions\n")
				for rs, pos in dup_rsID.items():
					fappend.write( "%s\t%s\t%s\n" % (rs, len(pos), ";".join(pos)) )
			### chrpos
			if dup_chrpos:
				fappend.write("### DUPLICATED chrpos ###\n")
				fappend.write("rsID\tCount\tPositions\n")
				for chrpos, rs_list in dup_chrpos.items():
					fappend.write( "%s\t%s\t%s\n" % (chrpos, len(rs_list), ";".join(rs_list)) )
			fappend.write( "#"*80 + '\n' )
		
		print "done in function get_duplicates()"



#print "************************************"
#print "%s dublicate in .bim file" % cols[1]
#print "Dublicates are at position: %s and %s" %	(snps_seen_rsID[cols[1]], cols[3])
#print "Retaining entry with lowest chromosomal position: %s" % keep


#Parse Arguments
arg_parser = argparse.ArgumentParser("Finds duplicates in genotype data and write them to filename '<input-prefix>.duplicates.txt' in the input dir. FURTHERMORE the script APPENDS to a file 'get_duplicates_append_file.txt' in the inputdir")
arg_parser.add_argument("--input", help="input .bim file", required=True) # e.g. /home/projects/tp/childrens/snpsnap/data/step1/full_no_pthin/CEU_GBR_TSI_unrelated.phase1.bim
args = arg_parser.parse_args()
print "got some argument..."


inputfile = args.input
(root, ext) = os.path.splitext(inputfile) #Split the pathname path into a pair (root, ext) such that root + ext == path, and ext is empty or begins with a period and contains at most one period.
	### EXAMPLE:
	#os.path.splitext('home/robert/Documents/Workspace/datafile.xlsx')
	#gives --> ('home/robert/Documents/Workspace/datafile', '.xlsx')
outputfile = root + '.duplicates_unique_rsID.txt' # e.g. /home/projects/tp/childrens/snpsnap/data/step1/full_no_pthin/CEU_GBR_TSI_unrelated.phase1.duplicates.txt

# appendfile = os.path.dirname(inputfile) + '/get_duplicates_append_file.txt' # SNPsnap version
appendfile = root + '.get_duplicates_append_file.txt'
print "appendfile={}".format(appendfile)

get_duplicates(inputfile, outputfile, appendfile)


print "script is done"
