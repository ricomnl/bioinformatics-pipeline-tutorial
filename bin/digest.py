#!/usr/bin/env python3
# All the code is taken from https://github.com/wfondrie/mokapot/blob/master/mokapot/parsers/fasta.py.
# The only reason why I didn't use mokapot as a package is because I wanted to leave any packaging 
# out of this tutorial.
import re
import sys


def load_fasta(filename):
	with open(filename, "r") as fasta_file:
		lines = fasta_file.read().splitlines()
	metadata = lines[0]
	sequence = "".join(lines[1:])
	return metadata, sequence


def save_peptides(filename, peptides):
	with open(filename, "w") as output_file:
		for peptide in peptides:
			output_file.write("{}\n".format(peptide))


def digest_protein(
	input_file, 
	output_file, 
	enzyme_regex="[KR]",
	missed_cleavages=0,
	min_length=4,
    max_length=75):
	metadata, sequence = load_fasta(input_file)

	# Find the cleavage sites
	enzyme_regex = re.compile(enzyme_regex)
	sites = (
		[0]
		+ [m.end() for m in enzyme_regex.finditer(sequence)]
		+ [len(sequence)]
	)

	peptides = set()

    # Do the digest
	for start_idx, start_site in enumerate(sites):
		for diff_idx in range(1, missed_cleavages + 2):
			end_idx = start_idx + diff_idx
			if end_idx >= len(sites):
				continue
			end_site = sites[end_idx]
			peptide = sequence[start_site:end_site]
			if len(peptide) < min_length or len(peptide) > max_length:
				continue
			peptides.add(peptide)
	save_peptides(output_file, peptides)


if __name__ == "__main__":
	input_file = sys.argv[1]
	output_file = sys.argv[2]
	digest_protein(input_file, output_file)