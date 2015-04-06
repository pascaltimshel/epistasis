#!/usr/bin/env python2.7

import os
import glob
import datetime
import time

import pdb

###################################### USAGE ######################################

###################################### REMARKS ######################################
### We use the *FIRST* of the position columns in the file_lan_et_al.

###################################### FILE SNIPPETS ######################################
### file_lan_et_al
# chr10   22930819        22930820        interaction11   chr15   75896781        75896782        interaction11   interaction_11
# chr5    111408972       111408973       interaction91   chr1    202167527       202167528       interaction91   interaction_85
# chr4    25934812        25934813        interaction228  chr5    133628590       133628591       interaction228  interaction_214
# chr17   41618048        41618049        interaction276  chr5    133625948       133625949       interaction276  interaction_258

### FORMAT fit-hi-c table [not used - only for FORMAT reference]
### /cvar/jhlab/timshel/egcut/interactome_fit-hi-c/interation_table.fit-hi-c.nosex.interchromosomal.hESC.q_1e-18.txt
# chrA    posA    chrB    posB    contactCount    p.value q.value interactionID   interactionID_final
# 10      100105000       12      55885000        8       4.373243e-29    6.503306e-21    interaction_317 interaction_1
# 2       212405000       10      100535000       8       6.88906e-27     4.923711e-19    interaction_1164        interaction_2
# 10      100905000       11      127055000       7       1.459951e-26    8.740568e-19    interaction_1793        interaction_3
# 3       185305000       10      107435000       21      5.821841e-87    1.578715e-77    interaction_15052       interaction_4

######################################  ######################################


###################################### SCRIPT ######################################
file_lan_et_al = "/cvar/jhlab/timshel/egcut/interactome/lift_findItersection.intersection.paste.clean.nosex.updatedIDs.interchromosomal"
file_out_interaction_table = "/cvar/jhlab/timshel/egcut/interactome/lift_findItersection.intersection.paste.clean.nosex.updatedIDs.interchromosomal.fit-hi-c"
#interation_table.fit-hi-c.nosex.interchromosomal.lan-et-al_K562.q_DUMMY.txt

chromosomes_allowed = range(1,23) # 1, 2, .., 22

#file_out_interaction_table_header = "chrA\tposA\tchrB\tposB\tcontactCount\tp.value\tq.value\tinteractionID\tinteractionID_final"
file_out_interaction_table_header = "chrA\tposA\tchrB\tposB\tcontactCount\tp.value\tq.value"

interactionID_count = 0
with open(file_lan_et_al, 'r') as f_in, open(file_out_interaction_table, 'w') as f_out:
	f_out.write(file_out_interaction_table_header+"\n") # write header
	for line in f_in:

		interactionID_count += 1

		line = line.strip()
		fields = line.split("\t")
		chrA, posA, chrB, posB, interactionID = fields[0], fields[1], fields[4], fields[5], fields[8]

		###################################### VALIDATION ######################################
		### Strip "chr" from chromosome and convert to int.
		# *IN THIS WAY WE KNOW WE ONlY GET CHR 1-22* 
		try:
			chrA_int = int(chrA.lstrip("chr"))
			chrB_int = int(chrB.lstrip("chr"))
		except:
			print "Cannot convert either chrA={} or chrB={} to int by stripping chars".format(chrA, chrB)
			print line
			raise

		if not (chrA_int in chromosomes_allowed):
			raise Exception("chrA_int is not in chromosomes_allowed: {}".format(chrA_int))
		if not (chrB_int in chromosomes_allowed):
			raise Exception("chrB_int is not in chromosomes_allowed: {}".format(chrB_int))

		### Safety check
		if chrA_int==chrB_int:
			raise Exception("ERROR chrA_int==chrB_int. line={}".format(line))
		######################################  ######################################


		#interactionID_final = "interaction_{count}".format(count=interactionID_count)

		### Write out
		# chrA    posA    chrB    posB    contactCount    p.value q.value interactionID   interactionID_final
		# line2write = "{chrA}\t{posA}\t{chrB}\t{posB}\t{contactCount}\t{p_value}\t{q_value}\t{interactionID}\t{interactionID_final}\n".format(chrA=chrA,
		# 																																	posA=posA,
		# 																																	chrB=chrB,
		# 																																	posB=posB,
		# 																																	contactCount="",
		# 																																	p_value="",
		# 																																	q_value="",
		# 																																	interactionID=interactionID,
		# 																																	interactionID_final=interactionID_final)

		line2write = "{chrA}\t{posA}\t{chrB}\t{posB}\t{contactCount}\t{p_value}\t{q_value}\n".format(chrA=chrA,
																									posA=posA,
																									chrB=chrB,
																									posB=posB,
																									contactCount="",
																									p_value="",
																									q_value="")

		f_out.write(line2write)



