#!/usr/bin/env python3
import sys


def load_counts(filename):
	with open(filename, "r") as count_file:
		counts = [line.split('\t') for line in count_file.read().splitlines()][0]
	return counts


def get_report(input_files):
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
		count_list.append([input_file.split("_")[0],] + counts)
	return count_list


def save_report(filename, report):
	with open(filename, "w") as output_file:
		for line in report:
			output_file.write("{}\n".format("\t".join(line)))


if __name__ == "__main__":
	input_files = sys.argv[1:-1]
	output_file = sys.argv[-1]
	report = get_report(input_files)
	save_report(output_file, report)