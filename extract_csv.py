#!/usr/bin/env python

import os
import sys
import glob
import argparse
import datetime
import time

time_start = time.time()
## USAGE:
# python extract_csv.py --file_input ./Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.txt.head10 --file_out ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.txt.head10.extract --columns ./Expression_related_docs_pascal/Pheno_Transposed.ExpressionData.txt.QuantileNormalized.Log2Transformed.sample_list_with_GTypes --rows probeID2arrayAddress_map_hemani_SNPpair-probe_all_501.txt.cut2 > extract_csv.out

## full file
# python extract_csv.py --file_input ./Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.txt --file_out ./Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.with_probes246.extract.txt --columns ./Expression_related_docs_pascal/Pheno_Transposed.ExpressionData.txt.QuantileNormalized.Log2Transformed.sample_list_with_GTypes --rows probeID2arrayAddress_map_hemani_SNPpair-probe_all_501.txt.cut2 > ./Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.with_probes246.extract.out

## full file (shorter name)
# python extract_csv.py --file_input ./Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.txt --file_out ./Expression_related_docs_pascal/ETypes.4GWASPCs.Quantile.Log2.Centered.ZTransformed.CovariatesRM.with_GTypes832.with_probes246.extract.txt --columns ./Expression_related_docs_pascal/Pheno_Transposed.ExpressionData.txt.QuantileNormalized.Log2Transformed.sample_list_with_GTypes --rows probeID2arrayAddress_map_hemani_SNPpair-probe_all_501.txt.cut2 > ./Expression_related_docs_pascal/ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.with_probes246.extract.out


def ParseArguments():
	arg_parser = argparse.ArgumentParser()
	arg_parser.add_argument("--file_input", help="input file", required=True)
	arg_parser.add_argument("--file_output", help="output file", required=True)
	arg_parser.add_argument("--rows", help="file to list of rows to include. Row file does not have to include 'row label'. It will not any effect since the header is processed seperately")
	arg_parser.add_argument("--columns", help="file to list of columns to include. Column file CAN include 'row label, e.g. LAB', but the row label will always be printed")
	arg_parser.add_argument("--delim", default='\t', help="delimiter to split on")

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
		if (fields[0] not in keep_rows) and (line_no != 0 ):
			continue

		rows_printed += 1
		## Write columns of interest
		cols2_write = [fields[i] for i in keep_columns_idx]
		fout.write( "{}\n".format(delim.join(cols2_write)) )

### Rows that where not found:
keep_rows_not_found = [keep for keep in keep_rows if keep not in rows_seen]

print "Columns not found: %s\n%s" % ( len(keep_columns_not_found), "\n".join(keep_columns_not_found) )
print "Rows not found: %s\n%s" % ( len(keep_rows_not_found), "\n".join(keep_rows_not_found) )

print "Number of rows printed (including header): %s" % rows_printed
print "Number of columns printed (including 'row labels') %s" % len(keep_columns_idx)
print "DIMENSION OF FILE: %s X %s" % (rows_printed, len(keep_columns_idx))

time_elapsed = time.time() - time_start
print "RUNTIME: {:.2f} sec ({:.2f} min)".format( time_elapsed, time_elapsed/60 ) # remember: ints can be formated as floats. floats cannot be formatted as ints

