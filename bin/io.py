#!/usr/bin/env python3


def load_fasta(filename):
	"""
	Load a protein with its metadata from a given .fasta file.
	"""
	with open(filename, "r") as fasta_file:
		lines = fasta_file.read().splitlines()
	metadata = lines[0]
	sequence = "".join(lines[1:])
	return metadata, sequence


def load_peptides(filename):
	"""
	Load peptides from a .txt file as a list.
	"""
	with open(filename, "r") as peptide_file:
		lines = peptide_file.read().splitlines()
	return lines


def load_counts(filename):
	with open(filename, "r") as count_file:
		counts = [line.split('\t') for line in count_file.read().splitlines()][0]
	return counts


def save_counts(filename, peptide_counts):
	"""
	Write out the peptide counts to a .tsv file using tabs as a separator.
	"""
	with open(filename, "w") as output_file:
		output_file.write("{}\n".format("\t".join([str(c) for c in peptide_counts])))


def save_peptides(filename, peptides):
	"""
	Write out the list of given peptides to a .txt file. Each line is a different peptide.
	"""
	with open(filename, "w") as output_file:
		for peptide in peptides:
			output_file.write("{}\n".format(peptide))


def save_report(filename, report):
	"""
	Save output report to a .tsv with given filename.
	"""
	with open(filename, "w") as output_file:
		for line in report:
			output_file.write("{}\n".format("\t".join(line)))
