PROTEINS := counts/KLF4.tsv \
			counts/MYC.tsv \
			counts/PO5F1.tsv \
			counts/SOX2.tsv

# Dummy targets
all: reports/protein_report.tsv

clean:
	rm -f counts/* fig/* peptides/* reports/*

.PHONY: all clean clean_nf

# Analysis and plotting
peptides/%.txt: bin/digest.py fasta/%.fasta
	$^ $@

counts/%.tsv: bin/count_amino_acids.py peptides/%.txt
	$^ $@

fig/%.png: bin/plot_count.py counts/%.tsv
	$^ $@

reports/protein_report.tsv: bin/get_report.py ${PROTEINS}
	$^ $@

# Misc
makeviz:
	make -Bnd | make2graph | dot -Tpng -o fig/out.png

clean_nf:
	rm -rf .nextflow.log*
	rm -rf work/*