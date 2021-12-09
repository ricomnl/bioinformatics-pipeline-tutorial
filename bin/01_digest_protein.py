#!/usr/bin/env python3
import argparse
import re


def load_fasta(filename):
	"""
	Load a protein with its metadata from a given .fasta file.
	"""
	with open(filename, "r") as fasta_file:
		lines = fasta_file.read().splitlines()
	metadata = lines[0]
	sequence = "".join(lines[1:])
	return metadata, sequence


def save_peptides(filename, peptides):
	"""
	Write out the list of given peptides to a .txt file. Each line is a different peptide.
	"""
	with open(filename, "w") as output_file:
		for peptide in peptides:
			output_file.write("{}\n".format(peptide))


def digest_protein(
	protein_sequence,
	enzyme_regex="[KR]",
	missed_cleavages=0,
	min_length=4,
    max_length=75
):
	"""
	Digest a protein into peptides using a given enzyme. Defaults to trypsin.
	All the code is taken from https://github.com/wfondrie/mokapot/blob/master/mokapot/parsers/fasta.py.
	The only reason why I didn't use mokapot as a package is because I wanted to leave 
	any packaging out of this tutorial.
	"""
	# Find the cleavage sites
	enzyme_regex = re.compile(enzyme_regex)
	sites = (
		[0]
		+ [m.end() for m in enzyme_regex.finditer(protein_sequence)]
		+ [len(protein_sequence)]
	)

	peptides = set()

    # Do the digest
	for start_idx, start_site in enumerate(sites):
		for diff_idx in range(1, missed_cleavages + 2):
			end_idx = start_idx + diff_idx
			if end_idx >= len(sites):
				continue
			end_site = sites[end_idx]
			peptide = protein_sequence[start_site:end_site]
			if len(peptide) < min_length or len(peptide) > max_length:
				continue
			peptides.add(peptide)
	return peptides


def main(input_file, output_file, **kwargs):
	"""
	Digest a protein and save the peptides.
	"""
	_, protein_sequence = load_fasta(input_file)
	peptides = digest_protein(protein_sequence, **kwargs)
	save_peptides(output_file, peptides)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"input_file", 
		help="A .fasta file containing a protein."
	)
	parser.add_argument(
		"output_file", 
		help="A .txt output file to write the peptides to."
	)
	parser.add_argument(
		"--enzyme_regex", 
		help="A regex for the enzyme to use for the digest. E.g. [KR] for trypsin."
	)
	parser.add_argument(
		"--missed_cleavages", 
		type=int, 
		help="Number of missed cleavages for the digest."
	)
	parser.add_argument(
		"--min_length", 
		type=int, 
		help="Minimum length for a peptide to be considered valid."
	)
	parser.add_argument(
		"--max_length", 
		type=int, 
		help="Maximum length for a peptide to be considered valid."
	)
	kwargs = {k:v for k,v in vars(parser.parse_args()).items() if v}
	input_file = kwargs.pop("input_file")
	output_file = kwargs.pop("output_file")
	main(input_file, output_file, **kwargs)
