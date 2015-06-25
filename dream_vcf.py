#!/usr/bin/env python

import sys, os
import re

try:
	import vcf
except ImportError:
	vcf = None

'''
Written by:
Shadrielle Melijah G. Espiritu, shadrielle.espiritu@oicr.on.ca
Ontario Institute for Cancer Research
'''

# ---------- LOCALS ---------- #
columns = ('#CHROM', 'POS', 'Sample', 'Binary.Cutoff', 'Continuous.Confidence.Score')
chroms = range(1,23)
chroms.extend(('X','Y'))
samples = ['IS1', 'IS2', 'IS3', 'IS4']


'''
Validates submission file by ensuring proper column names and values.
'''
def validate(infile):
	print "Starting Validation.\n"

	try:
		input_fh = open(infile, 'r')
	except IOError:
		sys.stderr.write("ERROR: can't find file or read data from: " + infile + "\n")
		sys.exit(1)
	
	# check list of pipelines
	pipe = re.compile('^##Pipelines=.+')
	pipelines_str = input_fh.readline()
	if pipe.match(pipelines_str) == None:
		sys.stderr.write("ERROR: list of pipelines not found\n")
		sys.stderr.write("\tFound " + pipelines_str)
		sys.stderr.write("\tRequires ##Pipelines=XXX1,XXX2,etc.\n")
		sys.exit(1)
	
	# ensure pipelines in correct format
	pipelines = (pipelines_str.split('=')[1]).split(',')

	skip = re.compile('^##')
	while(1):
		sub_columns = input_fh.readline()
		if skip.match(sub_columns):
			continue
		break
	
	fields = sub_columns.split('\t')
	for i in range(len(columns)):
		try:
			if fields[i].strip() != columns[i]:
				sys.stderr.write("Invalid column name: " + fields[i] + "\n")
				sys.stderr.write("Should be: " + columns[i] + "\n")
				sys.stderr.write("\tCheck whether you intended to have '##' instead of '#'\n")
				sys.exit(1)
		except:
			if columns[i] == 'Continuous.Confidence.Score':
				print "Optional confidence score not provided."
			else:
				sys.stderr.write("Missing column: " + columns[i] + "\n")
				sys.exit(1)
	
	print "Submitted file has proper column names.\n"
	print "Summarizing calls..."

	try:
		unique_sample = None
		positive_calls = 0
		negative_calls = 0

		for rec in input_fh:
			rec_fields = rec.split('\t')
			# check CHROM
			if rec_fields[0] not in str(chroms):
				sys.stderr.write("Invalid CHROM: " + rec_fields[0] + "\n")
				sys.stderr.write("\tValid options: [%s]\n" % ', '.join(map(str,chroms)))
				sys.exit(1)

			# check POS
			try:
				pos = int(rec_fields[1])
			except:
				sys.stderr.write("POS not an int: " + rec_fields[1] + "\n")
				sys.exit(1)

			# check sample name
			if unique_sample == None:
				unique_sample = rec_fields[2]
			
			if rec_fields[2] not in samples:
				sys.stderr.write("Invalid Sample: " + rec_fields[2] + "\n")
				sys.stderr.write("\tValid options: [%s]\n" % ', '.join(map(str,samples)))
				sys.exit(1)

			if rec_fields[2] != unique_sample:
				sys.stderr.write("Sample name not unique: " + rec_fields[2] + "\n")
				sys.stderr.write("\tHave been tracking " + unique_sample + "\n")
				sys.exit(1)

			# check binary cutoff
			try:
				binary = int(rec_fields[3])
			except:
				sys.stderr.write("Binary.Cutoff not an int: " + rec_fields[2] + "\n")
				sys.exit(1)

			if binary == 1:
				positive_calls += 1
			elif binary == 0:
				negative_calls += 1
			else:
				sys.stderr.write("Binary cutoff should be 0 or 1: " + rec_fields[3] + "\n")
				sys.exit(1)
	except:
		sys.stderr.write("\nValidation Failed.\n")
		sys.exit(1)
	
	print "Validation complete.\n"
	print "Total records: ", positive_calls+negative_calls
	print "-" * 60
	print "Positive count: ", positive_calls
	print "Negative count: ", negative_calls
	print "-" * 60
	

'''
Converts output format for IGCG-TCGA DREAM SMC-Meta to simplified vcf format
Generates new output file called sub.vcf
'''
def convert(infile, truth):
	print "Starting Conversion.\n"
	assert truth.endswith('.vcf') or truth.endswith('.vcf.gz')

	# try opening file before proceeding to process data
	try:
		input_fh = open(infile, 'r')
	except IOError:
		sys.stderr.write("ERROR: can't find file or read data!\n")
		sys.exit(1)
	output_fh = open('sub.vcf', 'w')

	# get unique sample name
	skip = re.compile('^#.*')
	with open(infile, 'r') as f:
		while(1):
			line = f.readline()
			if skip.match(line):
				continue
			unique_sample = line.split('\t')[2]
			break

	truvcfh = vcf.Reader(filename=truth)
	# compile list of true records, in case user's output doesn't have calls sorted in order
	trulist = [trurec for trurec in truvcfh]

	vcf_header = "##fileformat=VCFv4.1\n"
	vcf_sample = "##SAMPLE=<ID=" + unique_sample + ">\n"
	vcf_fields = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

	# get pipelines used
	pipelines_str = input_fh.readline()
	while(1):
		sub_columns = input_fh.readline()
		if skip.match(sub_columns):
			next
		break
	
	# writing meta data
	output_fh.write(vcf_header)
	output_fh.write(vcf_sample)
	output_fh.write(pipelines_str)
	output_fh.write(vcf_fields)

	# evaluating each call
	sub_columns = input_fh.readline()
	for line in input_fh:
		fields = line.split('\t')
		chrom = fields[0]
		chrom = chrom.replace("chr","")
		pos = fields[1]
		binary = fields[3]

		if int(binary):
			'''
			Fixing REF and ALT fields:
				For true positive calls, extract REF and ALT from truth vcf
				For false positive calls, set REF=N, ALT=A (does not matter which
				as long as they're not '.' or both the same)
			'''
			ref = "N"
			alt = "A"
			for trurec in trulist:
				if ((chrom == trurec.CHROM) and (int(pos) == trurec.POS)):
					ref = str(trurec.REF)
					alt = str(trurec.ALT[0])
					break

			output_fh.write("\t".join( (chrom,pos,".",ref,alt,".","PASS","SOMATIC") ) )
			output_fh.write("\n")

	input_fh.close()
	output_fh.close()

	print "Conversion complete. Converted VCF file: sub.vcf\n"


'''
Main file
./dream_vcf <submitted file> [<truth file]

Without truth file, performs VALIDATION
With truth file, performs CONVERSION
'''
if __name__ == '__main__':
	if len(sys.argv) == 2:
		sys.stderr.write("Usage: " + sys.argv[0] + " <submitted file> [<truth file>]\n")
		print "No truth file given. Will only be validating submitted file...\n"
		validate(sys.argv[1])
		sys.exit(1)
		
	elif len(sys.argv) == 3:
		print "Truth file given. Converting submitted file to vcf...\n"

		if vcf is None:
			print "Please install PyVCF"
			print ">>> pip install PyVCF"
			sys.exit(1)

		convert(sys.argv[1], sys.argv[2])

	else:
		sys.stderr.write("Usage: " + sys.argv[0] + " <submitted file> [<truth file>]\n")
		sys.exit(1)
