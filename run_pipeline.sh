#!/usr/bin/env bash
# USAGE: bash run_pipeline.sh

./scripts/digest.py fasta/KLF4.fasta peptides/KLF4_peptides.txt
./scripts/digest.py fasta/MYC.fasta peptides/MYC_peptides.txt
./scripts/digest.py fasta/PO5F1.fasta peptides/PO5F1_peptides.txt
./scripts/digest.py fasta/SOX2.fasta peptides/SOX2_peptides.txt

./scripts/count_amino_acids.py peptides/KLF4_peptides.txt counts/KLF4_aa_counts.txt
./scripts/count_amino_acids.py peptides/MYC_peptides.txt counts/MYC_aa_counts.txt
./scripts/count_amino_acids.py peptides/PO5F1_peptides.txt counts/PO5F1_aa_counts.txt
./scripts/count_amino_acids.py peptides/SOX2_peptides.txt counts/SOX2_aa_counts.txt

./scripts/plot_count.py counts/KLF4_aa_counts.txt fig/KLF4_aa_plot.png
./scripts/plot_count.py counts/MYC_aa_counts.txt fig/MYC_aa_plot.png
./scripts/plot_count.py counts/PO5F1_aa_counts.txt fig/PO5F1_aa_plot.png
./scripts/plot_count.py counts/SOX2_aa_counts.txt fig/SOX2_aa_plot.png

./scripts/get_report.py counts/KLF4_aa_counts.txt \
						counts/MYC_aa_counts.txt \
						counts/PO5F1_aa_counts.txt \
						counts/SOX2_aa_counts.txt
