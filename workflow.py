"""workflow.py"""
from pathlib import Path
import re
from typing import List, Tuple

from redun import task, File, functools


redun_namespace = "bioinformatics_pipeline_tutorial"


def load_fasta(input_file: File) -> Tuple[str, str]:
    """
    Load a protein with its metadata from a given .fasta file.
    """
    with input_file.open("r") as fasta_file:
        lines = fasta_file.read().splitlines()
    metadata = lines[0]
    sequence = "".join(lines[1:])
    return metadata, sequence


def save_peptides(filename: str, peptides: List[str]) -> File:
    """
    Write out the list of given peptides to a .txt file. Each line is a different peptide.
    """
    output_file = File(filename)
    with output_file.open("w") as out:
        for peptide in peptides:
            out.write("{}\n".format(peptide))
    return output_file


def digest_protein(
    protein_sequence: str,
    enzyme_regex: str = "[KR]",
    missed_cleavages: int =0,
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


def load_peptides(input_file: File) -> List[str]:
    """
    Load peptides from a .txt file as a list.
    """
    with input_file.open("r") as peptide_file:
        lines = peptide_file.read().splitlines()
    return lines


def save_counts(filename: str, peptide_counts: List[int]) -> File:
    """
    Write out the peptide counts to a .tsv file using tabs as a separator.
    """
    output_file = File(filename)
    with output_file.open("w") as out:
        out.write("{}\n".format("\t".join([str(c) for c in peptide_counts])))
    return output_file


def num_peptides(peptides: List[str]) -> int:
    """
    Retrieve the number of peptides in a given list.
    """
    return len(peptides)


def num_peptides_with_aa(peptides: List[str], amino_acid: str = "C") -> int:
    """
    Count the number of peptides in a given list that contain a given amino acid. 
    Defaults to cysteine.
    """
    return sum([1 if amino_acid in peptide else 0 for peptide in peptides])


def total_num_aa_in_protein(protein: str) -> int:
    """
    Count the total number of amino acids in a given protein string.
    """
    return len(protein)


def num_aa_in_protein(protein: str, amino_acid: str = "C") -> int:
    """
    Count the number of times a given amino acid occurs in a given protein.
    Defaults to cysteine.
    """
    return protein.count(amino_acid)


@task(version="1")
def digest_protein_task(input_fasta: File) -> File:
    _, protein_sequence = load_fasta(input_fasta)
    peptides = digest_protein(protein_sequence)
    protein = Path(input_fasta.path).stem
    output_path = Path("data").joinpath(f"{protein}.peptides.txt")
    peptides_file = save_peptides(str(output_path), peptides)
    return peptides_file


@task(version="1")
def count_amino_acids_task(input_fasta: File, input_peptides: File, amino_acid: str = "C") -> File:
    """
    Count the number of times a given amino acid appears in a protein as well
    as its peptides after digestion.
    """
    _, protein_sequence = load_fasta(input_fasta)
    peptides = load_peptides(input_peptides)
    n_peptides = num_peptides(peptides)
    n_peptides_with_aa = num_peptides_with_aa(peptides)
    total_aa_in_protein = total_num_aa_in_protein(protein_sequence)
    aa_in_protein = num_aa_in_protein(protein_sequence)
    protein = Path(input_fasta.path).stem
    output_path = Path("data").joinpath(f"{protein}.count.tsv")
    aa_count_file = save_counts(
        str(output_path), 
        [amino_acid, n_peptides, n_peptides_with_aa, total_aa_in_protein, aa_in_protein]
    )
    return aa_count_file


@task(version="1")
def main(input_fastas: List[File]) -> List[File]:
    peptide_files = [digest_protein_task(fasta) for fasta in input_fastas]
    aa_count_files = [count_amino_acids_task(fasta, peptides) for (fasta, peptides) in zip(input_fastas, peptide_files)]
    return aa_count_files