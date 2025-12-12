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
	print("Syntax:\tfix_strelka2_vcf_snvs.py [Input VCF.gz]\n")
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

			# Split the colon-delimited normal column
			normalArray = [splits for splits in strArray[9].split(":") if splits is not ""]

			# Extract values for ALT and REF
			if strArray[3] == "A":
				normalREF = int([splits for splits in normalArray[4].split(",") if splits is not ""][0])
			elif strArray[3] == "C":
				normalREF = int([splits for splits in normalArray[5].split(",") if splits is not ""][0])
			elif strArray[3] == "G":
				normalREF = int([splits for splits in normalArray[6].split(",") if splits is not ""][0])
			elif strArray[3] == "T":
				normalREF = int([splits for splits in normalArray[7].split(",") if splits is not ""][0])
			else:
				print("REF was not A, C, G or T\n")
				sys.exit()

			if strArray[4] == "A":
				normalALT = int([splits for splits in normalArray[4].split(",") if splits is not ""][0])
			elif strArray[4] == "C":
				normalALT = int([splits for splits in normalArray[5].split(",") if splits is not ""][0])
			elif strArray[4] == "G":
				normalALT = int([splits for splits in normalArray[6].split(",") if splits is not ""][0])
			elif strArray[4] == "T":
				normalALT = int([splits for splits in normalArray[7].split(",") if splits is not ""][0])
			else:
				print("ALT was not A, C, G or T\n")
				sys.exit(1)

			# Split the colon-delimited tumor column
			tumorArray = [splits for splits in strArray[10].split(":") if splits is not ""]

			# Extract values for ALT and REF
			if strArray[3] == "A":
				tumorREF = int([splits for splits in tumorArray[4].split(",") if splits is not ""][0])
			elif strArray[3] == "C":
				tumorREF = int([splits for splits in tumorArray[5].split(",") if splits is not ""][0])
			elif strArray[3] == "G":
				tumorREF = int([splits for splits in tumorArray[6].split(",") if splits is not ""][0])
			elif strArray[3] == "T":
				tumorREF = int([splits for splits in tumorArray[7].split(",") if splits is not ""][0])
			else:
				print("REF was not A, C, G or T\n")
				sys.exit(1)

			if strArray[4] == "A":
				tumorALT = int([splits for splits in tumorArray[4].split(",") if splits is not ""][0])
			elif strArray[4] == "C":
				tumorALT = int([splits for splits in tumorArray[5].split(",") if splits is not ""][0])
			elif strArray[4] == "G":
				tumorALT = int([splits for splits in tumorArray[6].split(",") if splits is not ""][0])
			elif strArray[4] == "T":
				tumorALT = int([splits for splits in tumorArray[7].split(",") if splits is not ""][0])
			else:
				print("ALT was not A, C, G or T\n")
				sys.exit(1)

			# Calcluate FA for normal
			normalAD = normalREF + normalALT
			if normalALT == 0:
				normalFA = 0
			else:
				normalFA = normalALT / normalAD

			# Calcluate AD and FA for tumor
			tumorAD = tumorREF + tumorALT
			if tumorALT == 0:
				tumorFA = 0
			else:
				tumorFA = tumorALT / tumorAD

			# Print line to output, add AD (split as REF,ALT) and FA accordingly, swap normal/tumor
			print(strArray[0] + "\t" \
			+ strArray[1] + "\t" \
			+ strArray[2] + "\t" \
			+ strArray[3] + "\t" \
			+ strArray[4] + "\t" \
			+ strArray[5] + "\t" \
			+ strArray[6] + "\t" \
			+ strArray[7] + "\t" \
			+ "AD:FA:" + strArray[8] + "\t" \
			+ str(tumorREF) + "," + str(tumorALT) + ":" + str(tumorFA) + ":"  + strArray[10] + "\t" \
			+ str(normalREF) + "," + str(normalALT) + ":" + str(normalFA) + ":"  + strArray[9])
