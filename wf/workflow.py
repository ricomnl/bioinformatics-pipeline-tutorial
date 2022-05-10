"""wf/workflow.py"""
import os
import re
from typing import List, Tuple

from latch import small_task, workflow
from latch.types import LatchFile
import matplotlib.pyplot as plt


def load_fasta(filename: str) -> Tuple[str, str]:
    """
    Load a protein with its metadata from a given .fasta file.
    """
    with open(filename, "r") as fasta_file:
        lines = fasta_file.read().splitlines()
    metadata = lines[0]
    sequence = "".join(lines[1:])
    return metadata, sequence


def load_peptides(filename: str) -> List[str]:
    """
    Load peptides from a .txt file as a list.
    """
    with open(filename, "r") as peptide_file:
        lines = peptide_file.read().splitlines()
    return lines


def load_counts(filename: str) -> List[List[str]]:
    with open(filename, "r") as count_file:
        counts = [line.split('\t') for line in count_file.read().splitlines()][0]
    return counts


def save_peptides(filename: str, peptides: List[str]) -> None:
    """
    Write out the list of given peptides to a .txt file. Each line is a different peptide.
    """
    with open(filename, "w") as output_file:
        for peptide in peptides:
            output_file.write("{}\n".format(peptide))


def save_counts(filename: str, peptide_counts: List[int]) -> None:
    """
    Write out the peptide counts to a .tsv file using tabs as a separator.
    """
    with open(filename, "w") as output_file:
        output_file.write("{}\n".format("\t".join([str(c) for c in peptide_counts])))


def plot_counts(filename: str, counts: List[str]) -> None:
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
    fig.savefig(filename)


def save_report(filename: str, report: List[List[str]]) -> None:
	"""
	Save output report to a .tsv with given filename.
	"""
	with open(filename, "w") as output_file:
		for line in report:
			output_file.write("{}\n".format("\t".join(line)))


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


def get_report(input_files: List[LatchFile]) -> List[List[str]]:
    """
    Generate output report for a given list of input files.
    """
    count_list = [
        [
            "Protein",
            "Target Amino Acid",
            "No. of Peptides",
            "No. of Peptides w/ Target Amino Acid",
            "Total No. of Amino Acids",
            "No. of Target Amino Acid",
        ]
    ]
    for input_file in input_files:
        counts = load_counts(input_file.local_path)
        count_list.append(
            [
                os.path.basename(input_file.local_path).split(".")[0],
            ]
            + counts
        )
    return count_list


@small_task
def digest_protein_task(input_fasta: LatchFile) -> LatchFile:
    protein = os.path.basename(input_fasta.local_path).split(".")[0]
    output_path = os.path.join(
        os.path.split(os.path.dirname(input_fasta.local_path))[0], "data", f"{protein}.peptides.txt"
    )
    _, protein_sequence = load_fasta(input_fasta.local_path)
    peptides = digest_protein(protein_sequence)
    save_peptides(output_path, peptides)
    return LatchFile(output_path, f"latch:///{output_path}")


@small_task
def count_amino_acids_task(
    input_fasta: LatchFile, input_peptides: LatchFile, amino_acid: str = "C"
) -> LatchFile:
    """
    Count the number of times a given amino acid appears in a protein as well
    as its peptides after digestion.
    """
    protein = os.path.basename(input_fasta.local_path).split(".")[0]
    output_path = os.path.join(
        os.path.split(os.path.dirname(input_fasta.local_path))[0], "data", f"{protein}.count.tsv"
    )
    _, protein_sequence = load_fasta(input_fasta)
    peptides = load_peptides(input_peptides)
    n_peptides = num_peptides(peptides)
    n_peptides_with_aa = num_peptides_with_aa(peptides, amino_acid=amino_acid)
    total_aa_in_protein = total_num_aa_in_protein(protein_sequence)
    aa_in_protein = num_aa_in_protein(protein_sequence, amino_acid=amino_acid)
    save_counts(
        output_path,
        [
            amino_acid,
            n_peptides,
            n_peptides_with_aa,
            total_aa_in_protein,
            aa_in_protein,
        ],
    )
    return LatchFile(output_path, f"latch:///{output_path}")


@small_task
def plot_count_task(input_count: LatchFile) -> LatchFile:
    """
    Load the calculated counts and create a plot.
    """
    protein = os.path.basename(input_count.local_path).split(".")[0]
    output_path = os.path.join(
        os.path.split(os.path.dirname(input_count.local_path))[0], "data", f"{protein}.plot.png"
    )
    counts = load_counts(input_count)
    plot_counts(output_path, counts)
    return LatchFile(output_path, f"latch:///{output_path}")


@small_task
def get_report_task(input_counts: List[LatchFile]) -> LatchFile:
    """
    Get a list of input files from a given folder and create a report.
    """
    output_path = os.path.join(
        os.path.split(os.path.dirname(input_counts[0].local_path))[0], "data", f"protein_report.tsv"
    )
    report = get_report(input_counts)
    save_report(output_path, report)
    return LatchFile(output_path, f"latch:///{output_path}")


@workflow
def bioinformatics_pipeline_tutorial(input_fasta: LatchFile) -> LatchFile:
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
    peptide_file = digest_protein_task(input_fasta=input_fasta)
    aa_count_file = count_amino_acids_task(input_fasta=input_fasta, input_peptides=peptide_file, amino_acid="C")
    count_plot_file = plot_count_task(input_count=aa_count_file)
    report_file = get_report_task(input_counts=[aa_count_file])
    return report_file


if __name__ == "__main__":
    bioinformatics_pipeline_tutorial(input_fasta=LatchFile("fasta/KLF4.fasta"))

