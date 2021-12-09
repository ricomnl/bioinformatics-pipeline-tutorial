#!/usr/bin/env python3
import sys


def load_peptides(filename):
	with open(filename, "r") as peptide_file:
		lines = peptide_file.read().splitlines()
	return lines


def save_peptide_counts(filename, peptide_counts):
	with open(filename, "w") as output_file:
		output_file.write("{}\n".format("\t".join([str(c) for c in peptide_counts])))


def num_peptides(peptides):
	return len(peptides)


def num_peptides_with_aa(peptides, amino_acid="C"):
	return sum([1 if amino_acid in peptide else 0 for peptide in peptides])


def total_num_aa_in_peptides(peptides):
	return sum([len(peptide) for peptide in peptides])


def num_aa_in_peptides(peptides, amino_acid="C"):
	return sum([peptide.count(amino_acid) for peptide in peptides])


def counts_for_peptides(
		input_file, 
		output_file,
		amino_acid="C"
	):
	peptides = load_peptides(input_file)
	n_peptides = num_peptides(peptides)
	n_peptides_with_aa = num_peptides_with_aa(peptides)
	total_aa_in_peptides = total_num_aa_in_peptides(peptides)
	aa_in_peptides = num_aa_in_peptides(peptides)
	save_peptide_counts(output_file, [amino_acid, n_peptides, n_peptides_with_aa, total_aa_in_peptides, aa_in_peptides])


if __name__ == "__main__":
	input_file = sys.argv[1]
	output_file = sys.argv[2]
	counts_for_peptides(input_file, output_file)