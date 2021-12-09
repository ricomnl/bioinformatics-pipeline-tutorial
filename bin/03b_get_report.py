#!/usr/bin/env python3
import argparse
import glob
import os

from bin.io import load_counts, save_report


def get_report(input_files):
	"""
	Generate output report for a given list of input files.
	"""
	count_list = [[
		"Protein",
		"Target Amino Acid",
		"No. of Peptides", 
		"No. of Peptides w/ Target Amino Acid", 
		"Total No. of Amino Acids", 
		"No. of Target Amino Acid"
	]]
	for input_file in input_files:
		counts = load_counts(input_file)
		count_list.append(
			[os.path.splitext(os.path.basename(input_file))[0],] + counts
		)
	return count_list


def main(input_folder, output_file):
	"""
	Get a list of input files from a given folder and create a report.
	"""
	input_files = glob.glob("{}/*.tsv".format(input_folder))
	report = get_report(input_files)
	save_report(output_file, report)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"input_folder", 
		help="A folder of .tsv files containing the peptide counts."
	)
	parser.add_argument(
		"output_file", 
		help="A .tsv output file to write the report to."
	)
	kwargs = {k:v for k,v in vars(parser.parse_args()).items() if v}
	input_folder = kwargs.pop("input_folder")
	output_file = kwargs.pop("output_file")
	main(input_folder, output_file)
	