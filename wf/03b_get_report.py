#!/usr/bin/env python3
import argparse
import glob
import os


def load_counts(filename):
	with open(filename, "r") as count_file:
		counts = [line.split('\t') for line in count_file.read().splitlines()][0]
	return counts


def save_report(filename, report):
	"""
	Save output report to a .tsv with given filename.
	"""
	with open(filename, "w") as output_file:
		for line in report:
			output_file.write("{}\n".format("\t".join(line)))


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
			[os.path.basename(input_file).split('.')[0],] + counts
		)
	return count_list


def main(input_files, output_file):
	"""
	Get a list of input files from a given folder and create a report.
	"""
	report = get_report(input_files)
	save_report(output_file, report)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"input_files",
		nargs="+",
		help="A list of .tsv files containing the peptide counts."
	)
	parser.add_argument(
		"--output_file", 
		help="A .tsv output file to write the report to."
	)
	kwargs = {k:v for k,v in vars(parser.parse_args()).items() if v}
	input_files = kwargs.pop("input_files")
	output_file = kwargs.pop("output_file")
	main(input_files, output_file)
	