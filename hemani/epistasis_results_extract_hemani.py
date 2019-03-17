#!/usr/bin/env python

import os
import sys
import glob
import argparse
import datetime
import time

import re

## USAGE:
# XXX

## EXAMPLE file_hemani_map:
# 5550253 ILMN_1651385    rs7989895       rs4846085
# 4040364 ILMN_1651705    rs872311        rs11032695
# 240097  ILMN_1651886    rs7108734       rs12784396
# 4010021 ILMN_1652333    rs898095        rs9892064

## ### .epi.qt example: /cvar/jhlab/timshel/egcut/epistasis_plink_246probes/990020/990020.epi.qt
# CHR1        SNP1 CHR2        SNP2     BETA_INT         STAT            P
#    1  rs12065581   13   rs7322768      -0.7644         15.9    6.687e-05
#    1   rs1320993   17  rs11150847       0.6257        16.51    4.836e-05
#    1   rs3738725   14  rs11628398        3.163        15.26    9.354e-05
#    1   rs1619856    2   rs1254900       0.6497        17.04    3.652e-05

snps_seen = set()

file_hemani_map = "/cvar/jhlab/timshel/egcut/hemani_SNPpair_probe_map.txt"

f_out = open('epistasis_results.out', 'w')
header = "CHR1\tSNP1\tCHR2\tSNP2\tBETA_INT\tSTAT\tP"
header += "\tprobe\tprobe_ILML"
f_out.write(header+"\n")

f_hemani = open(file_hemani_map, 'r')
for line_no, line in enumerate(f_hemani):
	#if line_no == 20: break
	(probe, probe_ILMN, snp1, snp2) = line.strip().split()
	#file_epi_qt = "/cvar/jhlab/timshel/egcut/epistasis_plink_246probes/{probe}/{probe}.epi.qt".format(probe=probe)
	file_epi_qt = "/cvar/jhlab/timshel/egcut/epistasis_plink_246probes_epil1eq1/{probe}/{probe}.epi.qt".format(probe=probe)
	print "probe %s | %s" % (probe, probe_ILMN)
	f_epi_qt = None
	try:
		f_epi_qt = open(file_epi_qt, 'r')
	except:	
		print "could not find file %s" % file_epi_qt
	else:
		for line in f_epi_qt:
			line = line.strip()
			#if snp1 in line or snp2 in line: print line 
			snps_seen.add(snp1)
			snps_seen.add(snp2)
			pattern = "(.*{snp1}.*{snp2})|(.*{snp2}.*{snp1})".format(snp1=snp1, snp2=snp2)
			if re.match(pattern, line):
				string2write = "\t".join(line.split() + [probe, probe_ILMN])
				f_out.write(string2write+"\n")

print "SNPs seen in .epi.qt files %s" % len(snps_seen)


## remember to write out what probe te pair came from

