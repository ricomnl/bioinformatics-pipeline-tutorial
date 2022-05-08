"""script_workflow.py"""
import os
from typing import List

from redun import task, File, script
from redun.file import glob_file


redun_namespace = "bioinformatics_pipeline_tutorial.script_workflow"


@task(version="1")
def digest_protein_task(
    input_fasta: File,
    enzyme_regex: str = "[KR]",
    missed_cleavages: int = 0,
    min_length: int = 4,
    max_length: int = 75,
) -> File:
    protein = input_fasta.basename().split(".")[0]
    output_path = os.path.join(
        os.path.split(input_fasta.dirname())[0], "data", f"{protein}.peptides.txt"
    )
    return script(
        f"""
        bin/01_digest_protein.py \
            {input_fasta.path} \
            {output_path} \
            --enzyme_regex {enzyme_regex} \
            --missed_cleavages {missed_cleavages} \
            --min_length {min_length} \
            --max_length {max_length}
        """,
        outputs=File(output_path),
    )


@task(version="1")
def count_amino_acids_task(
    input_fasta: File, input_peptides: File, amino_acid: str = "C"
) -> File:
    """
    Count the number of times a given amino acid appears in a protein as well
    as its peptides after digestion.
    """
    protein = input_fasta.basename().split(".")[0]
    output_path = os.path.join(
        os.path.split(input_fasta.dirname())[0], "data", f"{protein}.count.tsv"
    )
    return script(
        f"""
        bin/02_count_amino_acids.py \
            {input_fasta.path} \
            {input_peptides.path} \
            {output_path} \
            --amino_acid {amino_acid}
        """,
        outputs=File(output_path),
    )


@task(version="1")
def plot_count_task(input_count: File) -> File:
    """
    Load the calculated counts and create a plot.
    """
    protein = input_count.basename().split(".")[0]
    output_path = os.path.join(
        os.path.split(input_count.dirname())[0], "data", f"{protein}.plot.png"
    )
    return script(
        f"""
        bin/03a_plot_count.py \
            {input_count.path} \
            {output_path}
        """,
        outputs=File(output_path),
    )


@task(version="1")
def get_report_task(input_counts: List[File]) -> File:
    """
    Get a list of input files from a given folder and create a report.
    """
    output_path = os.path.join(
        os.path.split(input_counts[0].dirname())[0], "data", f"protein_report.tsv"
    )
    return script(
        f"""
        bin/03b_get_report.py \
            {' '.join([f.path for f in input_counts])} \
            --output_file {output_path}
        """,
        outputs=File(output_path),
    )


@task(version="1")
def archive_results_task(inputs_plots: List[File], input_report: File) -> File:
    output_path = os.path.join(
        os.path.split(input_report.dirname())[0], "data", f"results.tgz"
    )
    return script(
        f"""
        rm -rf results
        mkdir results
        cp {' '.join([f.path for f in inputs_plots] + [input_report.path])} results/
        tar -czf {output_path} results
        rm -r results
        """,
        outputs=File(output_path),
    )


@task(version="1")
def main(
    input_dir: str,
    amino_acid: str = "C",
    enzyme_regex: str = "[KR]",
    missed_cleavages: int = 0,
    min_length: int = 4,
    max_length: int = 75,
) -> List[File]:
    input_fastas = [File(f) for f in glob_file(f"{input_dir}/*.fasta")]
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
        count_amino_acids_task(
            fasta, peptides, amino_acid=amino_acid
        )
        for (fasta, peptides) in zip(input_fastas, peptide_files)
    ]
    count_plots = [
        plot_count_task(aa_count)
        for aa_count in aa_count_files
    ]
    report_file = get_report_task(aa_count_files)
    results_archive = archive_results_task(
        count_plots, report_file
    )
    return results_archive