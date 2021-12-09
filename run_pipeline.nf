#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def getProtein(fileName) {
  fileName.getBaseName().tokenize(".")[0]
}

process digestProtein {
  publishDir "peptides/", mode: "copy"
  
  input:
    path input_fasta

  output:
    path "*.txt"

  script:
    def protein = getProtein(input_fasta)
    """
    digest.py ${input_fasta} ${protein}.txt
    """
}


process countAA {
  publishDir "counts/", mode: "copy"
  
  input:
    path input_peptides

  output:
    path "*.tsv"

  script:
    def protein = getProtein(input_peptides)
    """
    count_amino_acids.py ${input_peptides} ${protein}.tsv
    """
}


process plotCount {
  publishDir "fig/", mode: "copy"
  
  input:
    path input_counts

  output:
    path "*.png"

  script:
    def protein = getProtein(input_counts)
    """
    plot_count.py ${input_counts} ${protein}.png
    """
}


process generateReport {
  publishDir "reports/", mode: "copy"
  
  input:
    path input_counts

  output:
    path "*.tsv"

  script:
    """
    get_report.py ${input_counts} protein_report.tsv
    """
}


workflow {
  fasta = Channel.fromPath("$baseDir/fasta/*.fasta")
  peptides = digestProtein(fasta)
  counts = countAA(peptides)
  plotCount(counts)
  generateReport(counts | collect)
}
