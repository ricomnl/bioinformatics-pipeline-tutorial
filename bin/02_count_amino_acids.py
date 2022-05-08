#!/usr/bin/env python3
import argparse


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


def save_counts(filename, peptide_counts):
	"""
	Write out the peptide counts to a .tsv file using tabs as a separator.
	"""
	with open(filename, "w") as output_file:
		output_file.write("{}\n".format("\t".join([str(c) for c in peptide_counts])))


def num_peptides(peptides):
	"""
	Retrieve the number of peptides in a given list.
	"""
	return len(peptides)


def num_peptides_with_aa(peptides, amino_acid="C"):
	"""
	Count the number of peptides in a given list that contain a given amino acid. 
	Defaults to cysteine.
	"""
	return sum([1 if amino_acid in peptide else 0 for peptide in peptides])


def total_num_aa_in_protein(protein):
	"""
	Count the total number of amino acids in a given protein string.
	"""
	return len(protein)


def num_aa_in_protein(protein, amino_acid="C"):
	"""
	Count the number of times a given amino acid occurs in a given protein.
	Defaults to cysteine.
	"""
	return protein.count(amino_acid)	


def main(input_fasta, input_peptides, output_file, amino_acid="C"):
	"""
	Count the number of times a given amino acid appears in a protein as well
	as its peptides after digestion.
	"""
	_, protein_sequence = load_fasta(input_fasta)
	peptides = load_peptides(input_peptides)
	n_peptides = num_peptides(peptides)
	n_peptides_with_aa = num_peptides_with_aa(peptides, **kwargs)
	total_aa_in_protein = total_num_aa_in_protein(protein_sequence)
	aa_in_protein = num_aa_in_protein(protein_sequence, **kwargs)
	save_counts(
		output_file, 
		[amino_acid, n_peptides, n_peptides_with_aa, total_aa_in_protein, aa_in_protein]
	)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"input_fasta", 
		help="A .fasta file containing a protein."
	)
	parser.add_argument(
		"input_peptides", 
		help="A .txt file containing peptides."
	)
	parser.add_argument(
		"output_file", 
		help="A .tsv output file to write the counts to."
	)
	parser.add_argument(
		"--amino_acid", 
		default="C", 
		help="The one letter code for the amino acid to target. Defaults to C for cysteine."
	)
	kwargs = {k:v for k,v in vars(parser.parse_args()).items() if v}
	input_fasta = kwargs.pop("input_fasta")
	input_peptides = kwargs.pop("input_peptides")
	output_file = kwargs.pop("output_file")
	main(input_fasta, input_peptides, output_file, **kwargs)
