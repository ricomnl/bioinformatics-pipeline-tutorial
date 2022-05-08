#!/usr/bin/env python3
import argparse

import matplotlib.pyplot as plt


def load_counts(filename):
	with open(filename, "r") as count_file:
		counts = [line.split('\t') for line in count_file.read().splitlines()][0]
	return counts


def plot_counts(counts):
	"""
	Plot the calculated counts.
	"""
	amino_acid, n_peptides, n_peptides_with_aa, total_aa_in_peptides, aa_in_peptides = counts
	labels_n_peptides = ["No. of Peptides", "No. of Peptides w/ {}".format(amino_acid)]
	labels_n_aa = ["Total No. of Amino Acids", "No. of {}'s".format(amino_acid)]
	colors = ["#001425", "#308AAD"]
	x = [1, 2]
	fig, (ax1, ax2) = plt.subplots(1, 2)
	ax1.bar(
		x,
		[int(n_peptides_with_aa), int(n_peptides)], 
		color=colors[0],
	)
	ax1.set_xticks(x)
	ax1.set_xticklabels(labels_n_peptides)
	ax2.bar(
		x,
		[int(aa_in_peptides), int(total_aa_in_peptides)], 
		color=colors[1],
	)
	ax2.set_xticks(x)
	ax2.set_xticklabels(labels_n_aa)
	fig.suptitle("{}'s in Peptides and Amino Acids".format(amino_acid))
	

def main(input_file, output_file):
	"""
	Load the calculated counts and create a plot.
	"""
	counts = load_counts(input_file)
	plot_counts(counts)
	if output_file == "show":
		plt.show()
	else:
		plt.savefig(output_file)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"input_file", 
		help="A .tsv file containing the peptide counts."
	)
	parser.add_argument(
		"output_file", 
		help="A .png output file to write the plot to. Use `show` to just show the output."
	)
	kwargs = {k:v for k,v in vars(parser.parse_args()).items() if v}
	input_file = kwargs.pop("input_file")
	output_file = kwargs.pop("output_file")
	main(input_file, output_file)
