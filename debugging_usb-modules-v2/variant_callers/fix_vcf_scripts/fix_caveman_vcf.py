#!/usr/bin/env python3
import sys
import os
import gzip

# This script expects a VCF.gz with NORMAL abd TUMOR produced by CaVEMan.
# It calculates the DP, AD and AF (written as FA to be compatible with the downstream modules).
# It swaps the normal/tumor columns to be compatible with the downstream modules.

if sys.version_info < (3, 0):
    sys.stdout.write("Sorry, requires Python 3, not Python 2\n")
    sys.exit(1)

# Check the number of command line arguments
if not len(sys.argv)==2:
	print("\nError:\tincorrect number of command-line arguments")
	print("Syntax:\tfix_caveman_vcf.py [Input VCF.gz]\n")
	sys.exit(1)

# Loop through each line in the input file
with gzip.open(sys.argv[1], 'rt') as fileInput:
	for strLine in fileInput:
		# Strip the endline character from each input line
		strLine = strLine.rstrip("\n")

		# The '#' character in VCF format indicates that the line is a header. Pass these to the new file
		# Add DP, AD and FA descriptions, and swap NORMAL/TUMOR
		if strLine.startswith("#"):
			if strLine.startswith("##FORMAT=<ID=PM"):
				print(strLine)
				print("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">")
				print("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allelic depths for the ref and alt alleles\">")
				print("##FORMAT=<ID=FA,Number=2,Type=Float,Description=\"Allele fractions of alternate alleles\">")
			elif strLine.startswith("#CHROM"):
				print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL")
			# Omit DS from the header ("DBSnp ID of known SNP") because in conflict with sufam (HaplotypeCaller) where DS means "downsampling". DS will not appear in the variants because we are not using DBSnp during the call
			elif not strLine.startswith("##INFO=<ID=DS"):
				print(strLine)
		else:
			# Split the entire tab-delimited line into an array
			strArray = [splits for splits in strLine.split("\t") if splits is not ""]

			# Split the colon-delimited normal column
			normalArray = [splits for splits in strArray[9].split(":") if splits is not ""]

			# Extract values that match to ALT and REF
			# original format is GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM
			if strArray[3] == "A":
				normalREF = (int([splits for splits in normalArray[1].split(",") if splits is not ""][0])
							+
							int([splits for splits in normalArray[5].split(",") if splits is not ""][0]))
			elif strArray[3] == "C":
				normalREF = (int([splits for splits in normalArray[2].split(",") if splits is not ""][0])
							+
							int([splits for splits in normalArray[6].split(",") if splits is not ""][0]))
			elif strArray[3] == "G":
				normalREF = (int([splits for splits in normalArray[3].split(",") if splits is not ""][0])
							+
							int([splits for splits in normalArray[7].split(",") if splits is not ""][0]))
			elif strArray[3] == "T":
				normalREF = (int([splits for splits in normalArray[4].split(",") if splits is not ""][0])
							+
							int([splits for splits in normalArray[8].split(",") if splits is not ""][0]))
			else:
				print("REF was not A, C, G or T\n")
				sys.exit()

			if strArray[4] == "A":
				normalALT = (int([splits for splits in normalArray[1].split(",") if splits is not ""][0])
							+
							int([splits for splits in normalArray[5].split(",") if splits is not ""][0]))
			elif strArray[4] == "C":
				normalALT = (int([splits for splits in normalArray[2].split(",") if splits is not ""][0])
							+
							int([splits for splits in normalArray[6].split(",") if splits is not ""][0]))
			elif strArray[4] == "G":
				normalALT = (int([splits for splits in normalArray[3].split(",") if splits is not ""][0])
							+
							int([splits for splits in normalArray[7].split(",") if splits is not ""][0]))
			elif strArray[4] == "T":
				normalALT = (int([splits for splits in normalArray[4].split(",") if splits is not ""][0])
							+
							int([splits for splits in normalArray[8].split(",") if splits is not ""][0]))
			else:
				print("ALT was not A, C, G or T\n")
				sys.exit(1)

			# Split the colon-delimited tumor column
			tumorArray = [splits for splits in strArray[10].split(":") if splits is not ""]

			# Extract values that match to ALT and REF
			# original format is GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM
			if strArray[3] == "A":
				tumorREF = (int([splits for splits in tumorArray[1].split(",") if splits is not ""][0])
							+
							int([splits for splits in tumorArray[5].split(",") if splits is not ""][0]))
			elif strArray[3] == "C":
				tumorREF = (int([splits for splits in tumorArray[2].split(",") if splits is not ""][0])
							+
							int([splits for splits in tumorArray[6].split(",") if splits is not ""][0]))
			elif strArray[3] == "G":
				tumorREF = (int([splits for splits in tumorArray[3].split(",") if splits is not ""][0])
							+
							int([splits for splits in tumorArray[7].split(",") if splits is not ""][0]))
			elif strArray[3] == "T":
				tumorREF = (int([splits for splits in tumorArray[4].split(",") if splits is not ""][0])
							+
							int([splits for splits in tumorArray[8].split(",") if splits is not ""][0]))
			else:
				print("REF was not A, C, G or T\n")
				sys.exit()

			if strArray[4] == "A":
				tumorALT = (int([splits for splits in tumorArray[1].split(",") if splits is not ""][0])
							+
							int([splits for splits in tumorArray[5].split(",") if splits is not ""][0]))
			elif strArray[4] == "C":
				tumorALT = (int([splits for splits in tumorArray[2].split(",") if splits is not ""][0])
							+
							int([splits for splits in tumorArray[6].split(",") if splits is not ""][0]))
			elif strArray[4] == "G":
				tumorALT = (int([splits for splits in tumorArray[3].split(",") if splits is not ""][0])
							+
							int([splits for splits in tumorArray[7].split(",") if splits is not ""][0]))
			elif strArray[4] == "T":
				tumorALT = (int([splits for splits in tumorArray[4].split(",") if splits is not ""][0])
							+
							int([splits for splits in tumorArray[8].split(",") if splits is not ""][0]))
			else:
				print("ALT was not A, C, G or T\n")
				sys.exit(1)

			# Calcluate DP for tumor and normal (sum up all Z values)
			normalDP = (int([splits for splits in normalArray[1].split(",") if splits is not ""][0])
						+
						int([splits for splits in normalArray[2].split(",") if splits is not ""][0])
						+
						int([splits for splits in normalArray[3].split(",") if splits is not ""][0])
						+
						int([splits for splits in normalArray[4].split(",") if splits is not ""][0])
						+
						int([splits for splits in normalArray[5].split(",") if splits is not ""][0])
						+
						int([splits for splits in normalArray[6].split(",") if splits is not ""][0])
						+
						int([splits for splits in normalArray[7].split(",") if splits is not ""][0])
						+
						int([splits for splits in normalArray[8].split(",") if splits is not ""][0]))


			tumorDP = (int([splits for splits in tumorArray[1].split(",") if splits is not ""][0])
						+
						int([splits for splits in tumorArray[2].split(",") if splits is not ""][0])
						+
						int([splits for splits in tumorArray[3].split(",") if splits is not ""][0])
						+
						int([splits for splits in tumorArray[4].split(",") if splits is not ""][0])
						+
						int([splits for splits in tumorArray[5].split(",") if splits is not ""][0])
						+
						int([splits for splits in tumorArray[6].split(",") if splits is not ""][0])
						+
						int([splits for splits in tumorArray[7].split(",") if splits is not ""][0])
						+
						int([splits for splits in tumorArray[8].split(",") if splits is not ""][0]))


			# Calcluate AD and FA for normal
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

			# Print to output.
			# Format will be GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM:DP:AD:FA (GT must be first if available, otherwise GATK complains downstream, easiest is to append new tags after the original entry)
			# Also swap normal/tumor columns
			# Also skip rare cases where REF and ALT are the same
			if strArray[3] != strArray[4]:
				print(strArray[0] + "\t" \
				+ strArray[1] + "\t" \
				+ strArray[2] + "\t" \
				+ strArray[3] + "\t" \
				+ strArray[4] + "\t" \
				+ strArray[5] + "\t" \
				+ strArray[6] + "\t" \
				+ strArray[7] + "\t" \
				+ strArray[8] + ":DP:AD:FA" + "\t" \
				+ strArray[10] + ":"  + str(tumorDP) + ":" + str(tumorREF) + "," + str(tumorALT) + ":" + str("{:.3f}".format(tumorFA)) + "\t" \
				+ strArray[9] + ":" + str(normalDP) + ":" + str(normalREF) + "," + str(normalALT) + ":" + str("{:.3f}".format(normalFA)) )
