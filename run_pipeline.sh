#!/usr/bin/env bash
# USAGE: bash run_pipeline.sh

mkdir -p data

# 01. Digest
bin/01_digest_protein.py fasta/KLF4.fasta data/KLF4.peptides.txt
bin/01_digest_protein.py fasta/MYC.fasta data/MYC.peptides.txt
bin/01_digest_protein.py fasta/PO5F1.fasta data/PO5F1.peptides.txt
bin/01_digest_protein.py fasta/SOX2.fasta data/SOX2.peptides.txt

# 02. Count
bin/02_count_amino_acids.py fasta/KLF4.fasta data/KLF4.peptides.txt data/KLF4.count.tsv
bin/02_count_amino_acids.py fasta/MYC.fasta data/MYC.peptides.txt data/MYC.count.tsv
bin/02_count_amino_acids.py fasta/PO5F1.fasta data/PO5F1.peptides.txt data/PO5F1.count.tsv
bin/02_count_amino_acids.py fasta/SOX2.fasta data/SOX2.peptides.txt data/SOX2.count.tsv

# 03a. Plot
bin/03a_plot_count.py data/KLF4.count.tsv data/KLF4.plot.png
bin/03a_plot_count.py data/MYC.count.tsv data/MYC.plot.png
bin/03a_plot_count.py data/PO5F1.count.tsv data/PO5F1.plot.png
bin/03a_plot_count.py data/SOX2.count.tsv data/SOX2.plot.png

# 03b. Generate Report
bin/03b_get_report.py data/KLF4.count.tsv \
					  data/MYC.count.tsv \
					  data/PO5F1.count.tsv \
					  data/SOX2.count.tsv \
					  --output_file=data/protein_report.tsv

# 04. Archive the results in a tarball so we can share them with a colleague
rm -rf results
mkdir results
cp data/*plot.png data/protein_report.tsv results/
tar -czf data/results.tgz results
rm -r results
