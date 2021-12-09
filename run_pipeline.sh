#!/usr/bin/env bash
# USAGE: bash run_pipeline.sh

# 01. Digest
mkdir -p peptides
bin/01_digest_protein.py fasta/KLF4.fasta peptides/KLF4.txt
bin/01_digest_protein.py fasta/MYC.fasta peptides/MYC.txt
bin/01_digest_protein.py fasta/PO5F1.fasta peptides/PO5F1.txt
bin/01_digest_protein.py fasta/SOX2.fasta peptides/SOX2.txt

# 02. Count
mkdir -p counts
bin/02_count_amino_acids.py fasta/KLF4.fasta peptides/KLF4.txt counts/KLF4.tsv
bin/02_count_amino_acids.py fasta/MYC.fasta peptides/MYC.txt counts/MYC.tsv
bin/02_count_amino_acids.py fasta/PO5F1.fasta peptides/PO5F1.txt counts/PO5F1.tsv
bin/02_count_amino_acids.py fasta/SOX2.fasta peptides/SOX2.txt counts/SOX2.tsv

# 03a. Plot
mkdir -p plots
bin/03a_plot_count.py counts/KLF4.tsv plots/KLF4.png
bin/03a_plot_count.py counts/MYC.tsv plots/MYC.png
bin/03a_plot_count.py counts/PO5F1.tsv plots/PO5F1.png
bin/03a_plot_count.py counts/SOX2.tsv plots/SOX2.png

# 03b. Generate Report
mkdir -p reports
bin/03b_get_report.py counts/ reports/protein_report.txt
