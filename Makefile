COUNT := data/KLF4.count.tsv data/MYC.count.tsv \
			data/PO5F1.count.tsv data/SOX2.count.tsv
PLOT := data/KLF4.plot.png data/MYC.plot.png \
			data/PO5F1.plot.png data/SOX2.plot.png

# Dummy targets
all: data/results.tgz

clean:
	rm -rf data/*

.PHONY: all clean

# Analysis and plotting
data/%.peptides.txt: bin/01_digest_protein.py fasta/%.fasta
	$^ $@

data/%.count.tsv: bin/02_count_amino_acids.py fasta/%.fasta data/%.peptides.txt
	$^ $@

data/%.plot.png: bin/03a_plot_count.py data/%.count.tsv
	$^ $@

data/protein_report.tsv: bin/03b_get_report.py ${COUNT}
	$^ --output_file=$@

# Archive for sharing
data/results.tgz: ${PLOT} data/protein_report.tsv
	rm -rf results
	mkdir results
	cp $^ results/
	tar -czf $@ results
	rm -r results
