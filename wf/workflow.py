"""wf/workflow.py"""
from pathlib import Path
import re
from typing import List, Tuple

from latch import small_task, workflow
from latch.types import LatchFile
from flytekit import dynamic, map_task


def load_fasta(filename: str) -> Tuple[str, str]:
	"""
	Load a protein with its metadata from a given .fasta file.
	"""
	with open(filename, "r") as fasta_file:
		lines = fasta_file.read().splitlines()
	metadata = lines[0]
	sequence = "".join(lines[1:])
	return metadata, sequence


def save_peptides(filename: str, peptides: List[str]) -> None:
	"""
	Write out the list of given peptides to a .txt file. Each line is a different peptide.
	"""
	with open(filename, "w") as output_file:
		for peptide in peptides:
			output_file.write("{}\n".format(peptide))


def digest_protein(
	protein_sequence: str,
	enzyme_regex: str ="[KR]",
	missed_cleavages: int = 0,
	min_length: int = 4,
    max_length: int = 75
) -> List[str]:
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


@small_task
def digest_protein_task(input_fasta: LatchFile) -> LatchFile:
    output_folder = f"{Path(input_fasta.remote_path).parts[1]}/data"
    output_path = f"{Path(input_fasta.local_path).stem}.peptides.txt"

    _, protein_sequence = load_fasta(input_fasta.local_path)
    peptides = digest_protein(protein_sequence)
    save_peptides(output_path, peptides)
    return LatchFile(output_path, f"latch:///{output_folder}/{output_path}")


@workflow
def bioinformatics_pipeline_tutorial(input_fasta: List[LatchFile]) -> List[LatchFile]:
    """This is a test workflow.

    Test Protein Workflow
    ----

    # 01 - Digest Protein

    * Digest a protein and save the peptides.

    __metadata__:
        display_name: Digest Protein
        author:
            name: Rico Meinl
            email: dev@rmeinl.com
            github: ricomnl
        repository: https://github.com/ricomnl/bioinformatics-pipeline-tutorial
        license:
            id: MIT

    Args:

        input_fasta:
            Input fasta file containing a single protein.

            __metadata__:
            display_name: Input fasta
    """
    return map_task(digest_protein_task)(input_fasta=input_fasta)



# @dynamic()
# def protein_workflow(input_fasta: List[LatchFile]) -> List[LatchFile]:
#    # return protein_workflow(input_fasta=input_fasta)
#     peptide_list = []
#     for fasta in input_fasta:
#         peptides = digest_protein_task(input_fasta=fasta)
#         peptide_list.append(peptides)
#     return peptide_list
