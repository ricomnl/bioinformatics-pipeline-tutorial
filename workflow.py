"""workflow.py"""
from enum import Enum
from typing import List

from redun import task, File
from redun.file import glob_file

from bioinformatics_pipeline_tutorial.lib import (
    digest_protein_task,
    count_amino_acids_task,
    plot_count_task,
    get_report_task,
    archive_results_task,
)


redun_namespace = "bioinformatics_pipeline_tutorial"


class Executor(Enum):
    default = "default"
    process = "process"
    batch = "batch"
    batch_debug = "batch_debug"


@task(version="1")
def main(
    input_dir: str,
    amino_acid: str = "C",
    enzyme_regex: str = "[KR]",
    missed_cleavages: int = 0,
    min_length: int = 4,
    max_length: int = 75,
    executor: Executor = Executor.default,
) -> List[File]:
    input_fastas = [File(f) for f in glob_file(f"{input_dir}/*.fasta")]
    peptide_files = [
        digest_protein_task.options(executor=executor.value)(
            fasta,
            enzyme_regex=enzyme_regex,
            missed_cleavages=missed_cleavages,
            min_length=min_length,
            max_length=max_length,
        )
        for fasta in input_fastas
    ]
    aa_count_files = [
        count_amino_acids_task.options(executor=executor.value)(
            fasta, peptides, amino_acid=amino_acid
        )
        for (fasta, peptides) in zip(input_fastas, peptide_files)
    ]
    count_plots = [
        plot_count_task.options(executor=executor.value)(aa_count)
        for aa_count in aa_count_files
    ]
    report_file = get_report_task.options(executor=executor.value)(aa_count_files)
    results_archive = archive_results_task.options(executor=executor.value)(
        count_plots, report_file
    )
    return results_archive
