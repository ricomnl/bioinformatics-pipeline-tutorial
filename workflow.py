"""workflow.py"""
from pathlib import Path
import re
import tarfile
from typing import List, Tuple

from plotly.subplots import make_subplots
import plotly.graph_objects as go
from redun import task, File


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


def load_peptides(input_file: File) -> List[str]:
    """
    Load peptides from a .txt file as a list.
    """
    with input_file.open("r") as peptide_file:
        lines = peptide_file.read().splitlines()
    return lines


def load_counts(input_file: File) -> List[List[str]]:
    """
    Load counts from a .tsv file as a list.
    """
    with input_file.open("r") as count_file:
        counts = [line.split("\t") for line in count_file.read().splitlines()][0]
    return counts


def save_peptides(filename: str, peptides: List[str]) -> File:
    """
    Write out the list of given peptides to a .txt file. Each line is a different peptide.
    """
    output_file = File(filename)
    with output_file.open("w") as out:
        for peptide in peptides:
            out.write("{}\n".format(peptide))
    return output_file


def save_counts(filename: str, peptide_counts: List[int]) -> File:
    """
    Write out the peptide counts to a .tsv file using tabs as a separator.
    """
    output_file = File(filename)
    with output_file.open("w") as out:
        out.write("{}\n".format("\t".join([str(c) for c in peptide_counts])))
    return output_file


def plot_counts(filename: str, counts: List[str]) -> File:
    """
    Plot the calculated counts.
    """
    (
        amino_acid,
        n_peptides,
        n_peptides_with_aa,
        total_aa_in_peptides,
        aa_in_peptides,
    ) = counts
    labels_n_peptides = ["No. of Peptides", "No. of Peptides w/ {}".format(amino_acid)]
    labels_n_aa = ["Total No. of Amino Acids", "No. of {}'s".format(amino_acid)]
    colors = ["#001425", "#308AAD"]
    fig = make_subplots(rows=1, cols=2)
    fig.add_trace(
        go.Bar(
            x=labels_n_peptides,
            y=[int(n_peptides_with_aa), int(n_peptides)],
            marker_color=colors[0],
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Bar(
            x=labels_n_aa,
            y=[int(aa_in_peptides), int(total_aa_in_peptides)],
            marker_color=colors[1],
        ),
        row=1,
        col=2,
    )
    fig.update_layout(
        height=600,
        width=800,
        title_text="{}'s in Peptides and Amino Acids".format(amino_acid),
        showlegend=False,
    )
    output_file = File(filename)
    fig.write_image(filename)
    return output_file


def save_report(filename: str, report: List[List[str]]) -> File:
    """
    Save output report to a .tsv with given filename.
    """
    output_file = File(filename)
    with output_file.open("w") as out:
        for line in report:
            out.write("{}\n".format("\t".join(line)))
    return output_file


def digest_protein(
    protein_sequence: str,
    enzyme_regex: str = "[KR]",
    missed_cleavages: int = 0,
    min_length: int = 4,
    max_length: int = 75,
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


def get_report(input_files: List[File]) -> List[List[str]]:
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
        counts = load_counts(input_file)
        count_list.append(
            [
                Path(input_file.path).stem.split(".")[0],
            ]
            + counts
        )
    return count_list


@task(version="1")
def digest_protein_task(
    input_fasta: File,
    enzyme_regex: str = "[KR]",
    missed_cleavages: int = 0,
    min_length: int = 4,
    max_length: int = 75,
) -> File:
    _, protein_sequence = load_fasta(input_fasta)
    peptides = digest_protein(
        protein_sequence,
        enzyme_regex=enzyme_regex,
        missed_cleavages=missed_cleavages,
        min_length=min_length,
        max_length=max_length,
    )
    protein = Path(input_fasta.path).stem
    output_path = Path("data").joinpath(f"{protein}.peptides.txt")
    peptides_file = save_peptides(str(output_path), peptides)
    return peptides_file


@task(version="1")
def count_amino_acids_task(
    input_fasta: File, input_peptides: File, amino_acid: str = "C"
) -> File:
    """
    Count the number of times a given amino acid appears in a protein as well
    as its peptides after digestion.
    """
    _, protein_sequence = load_fasta(input_fasta)
    peptides = load_peptides(input_peptides)
    n_peptides = num_peptides(peptides)
    n_peptides_with_aa = num_peptides_with_aa(peptides, amino_acid=amino_acid)
    total_aa_in_protein = total_num_aa_in_protein(protein_sequence)
    aa_in_protein = num_aa_in_protein(protein_sequence, amino_acid=amino_acid)
    protein = Path(input_fasta.path).stem
    output_path = Path("data").joinpath(f"{protein}.count.tsv")
    aa_count_file = save_counts(
        str(output_path),
        [
            amino_acid,
            n_peptides,
            n_peptides_with_aa,
            total_aa_in_protein,
            aa_in_protein,
        ],
    )
    return aa_count_file


@task(version="1")
def plot_count_task(input_count: File) -> File:
    """
    Load the calculated counts and create a plot.
    """
    counts = load_counts(input_count)
    protein = Path(input_count.path).stem
    output_path = Path("data").joinpath(f"{protein}.plot.png")
    counts_plot = plot_counts(str(output_path), counts)
    return counts_plot


@task(version="1")
def get_report_task(input_counts: List[File]) -> File:
    """
    Get a list of input files from a given folder and create a report.
    """
    report = get_report(input_counts)
    output_path = Path("data").joinpath("protein_report.tsv")
    report_file = save_report(str(output_path), report)
    return report_file


@task(version="1")
def archive_results_task(inputs_plots: List[File], input_report: File) -> File:
    output_path = Path("data").joinpath("results.tgz")
    tar_file = File(str(output_path))
    with tar_file.open("wb") as out:
        with tarfile.open(fileobj=out, mode="w|gz") as tar:
            for file_path in inputs_plots + [input_report]:
                tar.add(file_path.path)
    return tar_file


@task(version="1")
def main(
    input_fastas: List[File],
    amino_acid: str = "C",
    enzyme_regex: str = "[KR]",
    missed_cleavages: int = 0,
    min_length: int = 4,
    max_length: int = 75,
) -> List[File]:
    peptide_files = [
        digest_protein_task(
            fasta,
            enzyme_regex=enzyme_regex,
            missed_cleavages=missed_cleavages,
            min_length=min_length,
            max_length=max_length,
        )
        for fasta in input_fastas
    ]
    aa_count_files = [
        count_amino_acids_task(fasta, peptides, amino_acid=amino_acid)
        for (fasta, peptides) in zip(input_fastas, peptide_files)
    ]
    count_plots = [plot_count_task(aa_count) for aa_count in aa_count_files]
    report_file = get_report_task(aa_count_files)
    results_archive = archive_results_task(count_plots, report_file)
    return results_archive
