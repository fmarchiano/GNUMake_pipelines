#!/usr/bin/env python3
import sys
import os
import gzip

# This script expects a VCF.gz with 2 samples (normal and tumor) as produced by strelka2 in somatic mode.
# It calculates the AD and AF (written as FA, as in Mutect1, to be compatible with the downstream modules).
# It swaps the normal/tumor columns to be compatible with the downstream modules.

if sys.version_info < (3, 0):
    sys.stdout.write("Sorry, requires Python 3, not Python 2\n")
    sys.exit(1)

# Check the number of command line arguments
if not len(sys.argv)==2:
	print("\nError:\tincorrect number of command-line arguments")
	print("Syntax:\tfix_strelka2_vcf_indels.py [Input VCF.gz]\n")
	sys.exit(1)

# Loop through each line in the input file, and calculate AD and FA
with gzip.open(sys.argv[1], 'rt') as fileInput:
	for strLine in fileInput:
		# Strip the endline character from each input line
		strLine = strLine.rstrip("\n")

		# The '#' character in VCF format indicates that the line is a header. Ignore these and just output to the new file
		# Add AD and FA descriptions, and swap NORMAL/TUMOR
		if strLine.startswith("#"):
			if strLine.startswith("##FORMAT=<ID=DP"):
				print("##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"AD computed from tier 1\">")
				print("##FORMAT=<ID=FA,Number=2,Type=Float,Description=\"AF computed from tier 1\">")
				print(strLine)
			elif strLine.startswith("#CHROM"):
				print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL")
			else:
				print(strLine)
		else:
			# Split the tab-delimited line into an array
			strArray = [splits for splits in strLine.split("\t") if splits is not ""]

			# Split the colon-delimited normal column, and extract TAR and TIR
			normalArray = [splits for splits in strArray[9].split(":") if splits is not ""]
			normalTAR = int([splits for splits in normalArray[2].split(",") if splits is not ""][0])
			normalTIR = int([splits for splits in normalArray[3].split(",") if splits is not ""][0])

			# Split the colon-delimited tumor column, and extract TAR and TIR
			tumorArray = [splits for splits in strArray[10].split(":") if splits is not ""]
			tumorTAR = int([splits for splits in tumorArray[2].split(",") if splits is not ""][0])
			tumorTIR = int([splits for splits in tumorArray[3].split(",") if splits is not ""][0])

			# Calcluate AD and FA for normal
			normalAD = normalTIR + normalTAR
			if normalTIR == 0:
				normalFA = 0
			else:
				normalFA = normalTIR / normalAD

			# Calcluate AD and FA for tumor
			tumorAD = tumorTIR + tumorTAR
			if tumorTIR == 0:
				tumorFA = 0
			else:
				tumorFA = tumorTIR / tumorAD

			# Print line to output, add AD (split as REF,ALT) and FA accordingly
			print(strArray[0] + "\t" \
			+ strArray[1] + "\t" \
			+ strArray[2] + "\t" \
			+ strArray[3] + "\t" \
			+ strArray[4] + "\t" \
			+ strArray[5] + "\t" \
			+ strArray[6] + "\t" \
			+ strArray[7] + "\t" \
			+ "AD:FA:" + strArray[8] + "\t" \
			+ str(tumorTAR) + "," + str(tumorTIR) + ":" + str(tumorFA) + ":"  + strArray[10] + "\t" \
			+ str(normalTAR) + "," + str(normalTIR) + ":" + str(normalFA) + ":"  + strArray[9])
