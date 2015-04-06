#!/usr/bin/env python

import os
import sys
import glob
import argparse
import datetime
import time

time_start = time.time()
## USAGE:
## 246 probes - full filename
# python extract_csv.py --file_input ./Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.txt --file_out ./Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.with_probes246.extract.txt --columns ./Expression_related_docs_pascal/Pheno_Transposed.ExpressionData.txt.QuantileNormalized.Log2Transformed.sample_list_with_GTypes --rows probeID2arrayAddress_map_hemani_SNPpair-probe_all_501.txt.cut2 > ./Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.with_probes246.extract.out

## 246 probes - shorter filename (NOT USED!)
# python extract_csv.py --file_input ./Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.txt --file_out ./Expression_related_docs_pascal/ETypes.4GWASPCs.Quantile.Log2.Centered.ZTransformed.CovariatesRM.with_GTypes832.with_probes246.extract.txt --columns ./Expression_related_docs_pascal/Pheno_Transposed.ExpressionData.txt.QuantileNormalized.Log2Transformed.sample_list_with_GTypes --rows probeID2arrayAddress_map_hemani_SNPpair-probe_all_501.txt.cut2 > ./Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.with_probes246.extract.out


## All probes - (no rows removed). New paths
# python extract_csv.py --file_input /Users/pascaltimshel/Dropbox/EGCUT_DATA/Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.txt --file_out /Users/pascaltimshel/Dropbox/EGCUT_DATA/Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.extract.txt --columns /Users/pascaltimshel/Dropbox/EGCUT_DATA/Expression_related_docs_pascal/Pheno_Transposed.ExpressionData.txt.QuantileNormalized.Log2Transformed.sample_list_with_GTypes > /Users/pascaltimshel/Dropbox/EGCUT_DATA/Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.extract.out


#python extract_csv.py --file_input /Users/pascaltimshel/Dropbox/Ubuntu/peer_analysis/export_residuals/egcut.peer_residuals_log2_k50.top50_mean_top50_var_refseq.txt --file_output /Users/pascaltimshel/Dropbox/Ubuntu/peer_analysis/export_residuals/egcut.peer_residuals_log2_k50.top50_mean_top50_var_refseq.hemani_102-246probes.txt --columns /Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/hemani_maps_res/hemani_probes_unique_tmp.txt --delim ' '


def ParseArguments():
	arg_parser = argparse.ArgumentParser("""
		This program is a memory efficient way of subsetting a tabular file with rows and/or column names. 
		The script reads and write line-by-line to keep memory usage low.
		The script IS not dependent on pandas package.
		"""
		)
	arg_parser.add_argument("--file_input", help="input file", required=True)
	arg_parser.add_argument("--file_output", help="output file", required=True)
	arg_parser.add_argument("--rows", help="file of list of rows to include. Row file does not have to include 'row label'. It will not have any effect since the header is processed seperately.")
	arg_parser.add_argument("--columns", help="file of list of columns to include. Column file can include 'row label, e.g. 'LAB', but it will not have any effect since the header is processed seperately (row label will always be printed).")
	arg_parser.add_argument("--delim", default='\t', help="delimiter to split on. The argument value must be a string, e.g. ' ' (SINGLE whitespace) or ',' (comma). Default is to split on tab ('\t').")

	args = arg_parser.parse_args()
	return args

def CheckArguments(args):
	if not (args.rows or args.columns):
		raise Exception("either row or columns files must be specified")

def LogArguments(args):
	# PRINT RUNNING DESCRIPTION 
	now = datetime.datetime.now()
	print '# ' + ' '.join(sys.argv)
	print '# ' + now.strftime("%a %b %d %Y %H:%M")
	print '# CWD: ' + os.getcwd()
	print '# COMMAND LINE PARAMETERS SET TO:'
	for arg in dir(args):
		if arg[:1]!='_':
			print '# \t' + "{:<30}".format(arg) + "{:<30}".format(getattr(args, arg))



def read_list(filename):
	with open(filename, 'r') as f:
		lines = f.read().splitlines()
		lines_set = set(lines)
		print "File %s: %s lines | %s unique lines" % ( os.path.basename(filename), len(lines), len(lines_set) )
		return lines_set


args = ParseArguments()
CheckArguments(args)
LogArguments(args)

delim = args.delim
print r"PRINTING RAW STRING | BEFORE decoding. delim=[{}]".format(delim)
### Decoding delim
delim = delim.decode('string_escape')
print r"PRINTING RAW STRING | AFTER decoding. delim=[{}]".format(delim)
### SEE:
# 1) http://stackoverflow.com/questions/5753332/passing-meta-characters-to-python-as-arguments-from-command-line
	# -->	r'\t\n\v\r'.decode('string-escape')
	#		'\t\n\x0b\r'
# 2) http://stackoverflow.com/questions/4020539/process-escape-sequences-in-a-string-in-python


file_input = args.file_input
file_output = args.file_output

## Read files
keep_columns = None
keep_rows = None
if args.columns:
	keep_columns = read_list(args.columns) # set
if args.rows:
	keep_rows = read_list(args.rows) # set

#http://code.activestate.com/recipes/577943-find-multiple-elements-in-a-list/
#find = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem]

keep_columns_not_found = []

rows_seen = []
rows_printed = 0
keep_columns_idx = None
with open(file_input, 'r') as fin, open(file_output, 'w') as fout: # multiple context with statement (Python 2.7 and up!)
	for line_no, line in enumerate(fin):
		line = line.strip('\n')
		fields = line.split(delim) 
		#str.split([sep[, maxsplit]]): 
			#OBS1: Splitting an empty string with a specified separator returns ['']
			#OBS2: 
			# If sep is not specified or is None, a different splitting algorithm is applied: 
			# runs of consecutive whitespace are regarded as a single separator, 
			# and the result will contain no empty strings at the start or end if the string has leading or trailing whitespace

		## Keeping track of rows seen
		rows_seen.append(fields[0])


		### Identify columns to keep
		if line_no == 0 and args.columns: # header line
			keep_columns_idx = [idx for idx, field in enumerate(fields) if field in keep_columns]
			keep_columns_not_found = [keep for keep in keep_columns if keep not in fields]

			if 0 not in keep_columns_idx: # we do not want the first column ('row label') to be printed twice. But we do want to make sure it is printed
				keep_columns_idx.insert(0,0) #list.insert(index, obj) --> inserts index 0 (zero) in the 'front' of the list
			
			if keep_columns_not_found:
				print "Warning: %s columns in your list where not found in the data file" % len(keep_columns_not_found)
			print "Will keep %s columns" % len(keep_columns_idx)
			#fields_set = set(fields) # ---> try to take intersection. [Will cause problem if there are multiple columns with same ID]
			#keep_columns_idx = [fields.index(field) for field in fields if field in keep_columns] # this is NOT tested
			
			## NOT TESTED
			# for keep in keep_columns:
			# 	if keep in fields:
			# 		keep_columns_idx.append(fields.index(keep)) # note: index returns index first item in list
			# 	else:
			# 		keep_columns_not_found.append(keep)

		### Skipping rows (except header, if it is not present in "keep_rows")
		if (args.rows) and (fields[0] not in keep_rows) and (line_no != 0 ):
			continue

		rows_printed += 1
		## Write columns of interest
		cols2_write = [fields[i] for i in keep_columns_idx]
		fout.write( "{}\n".format(delim.join(cols2_write)) )

### Rows that where not found:
if args.rows:
	keep_rows_not_found = [keep for keep in keep_rows if keep not in rows_seen]
	print "Rows not found: %s\n%s" % ( len(keep_rows_not_found), "\n".join(keep_rows_not_found) )
if args.columns:
	print "Columns not found: %s\n%s" % ( len(keep_columns_not_found), "\n".join(keep_columns_not_found) )


print "Number of rows printed (including header): %s" % rows_printed
print "Number of columns printed (including 'row labels') %s" % len(keep_columns_idx)
print "DIMENSION OF FILE: %s X %s" % (rows_printed, len(keep_columns_idx))

time_elapsed = time.time() - time_start
print "RUNTIME: {:.2f} sec ({:.2f} min)".format( time_elapsed, time_elapsed/60 ) # remember: ints can be formated as floats. floats cannot be formatted as ints

