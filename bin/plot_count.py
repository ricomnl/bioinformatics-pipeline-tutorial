#!/usr/bin/env python3
import sys

import matplotlib.pyplot as plt


def load_counts(filename):
	with open(filename, "r") as count_file:
		counts = count_file.read().split('\t')
	return counts


def plot_counts(counts):
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


if __name__ == "__main__":
	input_file = sys.argv[1]
	output_file = sys.argv[2]

	counts = load_counts(input_file)
	plot_counts(counts)
	if output_file == "show":
		plt.show()
	else:
		plt.savefig(output_file)