#!/usr/bin/env python

import sys, os
import re
import tempfile

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
columns = ('#CHROM', 'POS', 'Sample', 'Predicted', 'Probability')
chroms = range(1,23)
chroms.extend(('X','Y'))
samples = ['synthetic_1', 'synthetic_2', 'synthetic_3', 'synthetic_4', 'CPCG0100', 'CPCG0183', 'CPCG0184', 'CPCG0196', 'CPCG0235', 'PCSI0023', 'PCSI0044', 'PCSI0046', 'PCSI0048', 'PCSI0072']

sub = 'submitted'
confs = {}
for sample in samples:
	confs[sample] = 'synthetic.challenge.set4.tumour.25pctmasked.truth.vcf.gz'


'''
Validates submission file by
- ensuring proper column names and values.
- ensure ALL tumour samples are included in submission
'''
def validate(infile):
	print "Starting Validation.\n"

	try:
		input_fh = open(infile, 'r')
	except IOError, e:
		raise IOError(e)
	
	unique_samples = [];
	# get unique sample names
	skip = re.compile('^#.*')
	try:
		with open(infile, 'r') as f:
			for line in f:
				if not skip.match(line):	
					unique_sample = line.split('\t')[2]

					assert unique_sample in samples, "\nInvalid sample: " +unique_sample+ "\n"
					if unique_sample not in unique_samples:
						unique_samples.append(unique_sample)

		assert len(samples) == len(unique_samples), "\nMissing tumour samples: [%s]\n" % ', '.join(list(set(samples) - set(unique_samples)))

	except Exception, e:
		print "Valid samples: [%s]\n" % ', '.join(map(str,samples))
		raise Exception(e)
		

	# init num_pipelines to -1
	num_pipelines = -1

	try:
		# check pipelines used for each unique sample
		for sample in unique_samples:
			sample_pipelines_re = re.compile('^##'+sample+'_Pipelines=.*')
			sample_pipelines = input_fh.readline()
			if sample_pipelines_re.match(sample_pipelines) is None:
				raise Exception("\nFound " +sample_pipelines+ "\nRequires ##" +sample+ "_Pipelines={list of pipelines}\n")
	
			# validate correct number of specified pipelines
			pipelines_str = sample_pipelines.split('=')[1]
			pipelines_list = pipelines_str.split(',')
			if (num_pipelines == -1):
				num_pipelines = len(pipelines_list)
			else:
				assert len(pipelines_list) == num_pipelines, "\nExpecting " +str(num_pipelines)+ " pipelines, but seen " +str(len(pipelines_list))+ " for sample " +sample
	except Exception, e:
		raise Exception(e)
	
	skip = re.compile('^##')
	while(1):
		line = input_fh.readline()
		if skip.match(line):
			continue
		break

	# should be at columns 
	fields = line.split('\t')
	for i in range(len(columns)):
		try:
			if fields[i].strip() != columns[i]:
				raise Exception("\nInvalid column name: " +fields[i]+ "\nShould be: " +columns[i])
		except Exception, e:
			if columns[i] == 'Probability':
				print "Optional probability not provided."
			else:
				raise Exception("\nMissing column: " +columns[i])
	
	print "Submitted file has proper column names.\n"
	print "Summarizing calls..."

	try:
		positive_calls = {}
		negative_calls = {}

		for rec in input_fh:
			rec_fields = rec.split('\t')
			# check CHROM
			assert rec_fields[0] in str(chroms), "\nInvalid CHROM: " +rec_fields[0]+ "\nValid options: [%s]\n" % ', '.join(map(str,chroms))
			# check POS
			pos = int(rec_fields[1])
			# get sample
			sample = rec_fields[2]
			# check prediction
			binary = int(rec_fields[3])

			assert binary == 1 or binary == 0, "\nPrediction should be 0 or 1: " +rec_fields[3]

			try:
				positive_calls[sample] += binary
				negative_calls[sample] += (1-binary)
			except KeyError:
				positive_calls[sample] = binary
				negative_calls[sample] = 1 - binary
	except Exception, e:
		raise Exception(e)

	for sample in samples:
		print "Total " +sample+ " records: ", positive_calls[sample]+negative_calls[sample]
		print "-" * 60
		print "Positive count: ", positive_calls[sample]
		print "Negative count: ", negative_calls[sample]
		print "-" * 60

	print "Validation complete.\n"


'''
Splits submitted file by sample
'''
def split(infile):
	print "Starting splitting.\n"
	
	# try opening file before proceeding to process data
	try:
		input_fh = open(infile, 'r')
	except IOError, e:
		raise IOError(e)

	pipelines = re.compile('^##.*_Pipelines=.*')
	header = re.compile('^#[^#].*')

	# generate output file for each sample
	output_fhs = {}
	for sample in samples:
		output_fhs[sample] = tempfile.NamedTemporaryFile()

	with open(infile, 'r') as f:
		for line in f:
			if header.match(line):
				for sample in samples:
					output_fhs[sample].write(line)
			elif pipelines.match(line):
				# get sample name associated with pipeline string
				temp = line.split('_Pipelines')[0]
				sample_name = temp.split('##')[1]
				output_fhs[sample_name].write(line)
			else:
				sample_name = line.split('\t')[2]
				output_fhs[sample_name].write(line)

	input_fh.close()
	print "Split complete.\n"
	return output_fhs


'''
Converts output format for IGCG-TCGA DREAM SMC-Meta to simplified vcf format
Generates new output files called sample_name.vcf
'''
def convert(infile, truth):
	print "Starting Conversion.\n"
	assert truth.endswith('.vcf') or truth.endswith('.vcf.gz')

	infile.seek(0)
	# get unique sample name
	skip = re.compile('^#.*')
	for line in infile:
		if skip.match(line):
			continue
		unique_sample = line.split('\t')[2]
		break

	output_dir = tempfile.mkdtemp()
	output_filepath = os.path.join(output_dir, unique_sample+".vcf")
	with open(output_filepath, 'w') as output_fh:

		truvcfh = vcf.Reader(filename=truth)
		# compile list of true records, in case user's output doesn't have calls sorted in order
		trulist = [trurec for trurec in truvcfh]
	
		vcf_header = "##fileformat=VCFv4.1\n"
		vcf_sample = "##SAMPLE=<ID=" + unique_sample + ">\n"
		vcf_fields = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	
		# writing meta data
		output_fh.write(vcf_header)
		output_fh.write(vcf_sample)
	
		# meta data associated with pipelines 
		pipelines = re.compile('^##.*')
		infile.seek(0)
		for line in infile:
			if pipelines.match(line):
				output_fh.write(line)
			else:
				break
		output_fh.write(vcf_fields)
	
		for line in infile:
			if skip.match(line):
				continue
	
			# evaluating each call
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

	infile.close()
	print "Conversion complete for " + unique_sample + ".\n"
	return output_filepath


'''
Preprocesses submitted file into 14 different VCFs, one for each tumour sample
- splits submission by sample
- converts each sample calls into vcf format
'''
def preprocess(infile, evaluation_confs):
	out_split = split(infile)
	output_list = list()

	for sample in samples:
		vcf = convert(out_split[sample], evaluation_confs[sample])
		output_list.append((vcf, sample))
#		output_list.append((vcf.name, sample))

	print "Preprocessing complete.\n"
	return output_list


'''
Main file
./dream_vcf.py validate <submitted file>
'''
if __name__ == '__main__':
	if len(sys.argv) != 3 or sys.argv[1] != "validate":
		sys.stderr.write("Usage:\n")
		sys.stderr.write("	./dream_vcf.py validate <submitted file>\n")
		print "Exitting...\n";
		sys.exit(1)

	else:
		print "Will only be validating submitted file...\n"
		validate(sys.argv[2])
	
	print "DONE\n"
