#!/usr/bin/env python3

import argparse
from SigProfilerAssignment import Analyzer as Analyze

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument('--samples', default=None, help='Path to the input somatic mutations file (if using segmentation file/mutational matrix) or input folder (mutation calling file/s).')
parser.add_argument('--output', default=None, help='Path to the output folder.')
parser.add_argument('--input_type', default='vcf', help='Can be "vcf", "matrix", and "seg:TYPE", where TYPE={"ASCAT", "ASCAT_NGS", "SEQUENZA", "ABSOLUTE", "BATTENBERG", "FACETS", "PURPLE", "TCGA"}')
parser.add_argument('--context_type', default='96', help='Required context type if input_type is "vcf". context_type takes which context type of the input data is considered for assignment. Valid options include "96", "288", "1536", "DINUC", and "ID". The default value is "96".')
parser.add_argument('--cosmic_version', default=3.4, help='version of the COSMIC reference signatures')
parser.add_argument('--exome', action=argparse.BooleanOptionalAction, default=False, help='if the exome renormalized COSMIC signatures will be used')
parser.add_argument('--genome_build', default="GRCh38", help='Supported genomes include "GRCh37", "GRCh38", "mm9", "mm10" and "rn6"')
parser.add_argument('--signature_database', default=None, help='Path to the input set of known mutational signatures (only in case that COSMIC reference signatures are not used)')
parser.add_argument('--exclude_signature_subgroups', default=None, help='Removes the signatures corresponding to specific subtypes to improve refitting (only available when using default COSMIC reference signatures).')
parser.add_argument('--export_probabilities', action=argparse.BooleanOptionalAction, default=True, help='Defines if the probability matrix per mutational context for all samples is created. The default value is True')
parser.add_argument('--export_probabilities_per_mutation', action=argparse.BooleanOptionalAction, default=True, help='Defines if the probability matrices per mutation for all samples are created. Only available when input_type is "vcf".')
parser.add_argument('--make_plots', action=argparse.BooleanOptionalAction, default=True, help='Toggle on and off for making and saving plots.')
parser.add_argument('--sample_reconstruction_plots', default='pdf', help='Select the output format for sample reconstruction plots. Valid inputs are {"pdf", "png", "both", None}')

args = parser.parse_args()


Analyze.cosmic_fit(samples=args.samples,
                   output=args.output,
                   input_type=args.input_type,
                   context_type=args.context_type,
                   collapse_to_SBS96=True,
                   cosmic_version=args.cosmic_version,
                   exome=args.exome,
                   genome_build=args.genome_build,
                   signature_database=args.signature_database,
                   exclude_signature_subgroups=args.exclude_signature_subgroups,
                   export_probabilities=args.export_probabilities,
                   export_probabilities_per_mutation=args.export_probabilities_per_mutation,
                   make_plots=args.make_plots,
                   sample_reconstruction_plots=args.sample_reconstruction_plots,
                   verbose=True)
