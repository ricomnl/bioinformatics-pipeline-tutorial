#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


// Helper function to extract the protein name from the filename
def getProtein(fileName) {
  fileName.getBaseName().tokenize(".")[0]
}


// Digest a protein and save the peptides
process digestProtein {
  input:
    path input_fasta

  output:
    path "*.txt"

  script:
    def protein = getProtein(input_fasta)
    """
    01_digest_protein.py ${input_fasta} ${protein}.peptides.txt
    """
}


// Count the number of times a given amino acid appears in a protein as well 
// as its peptides after digestion
process countAA {  
  input:
    path input_peptides
    path input_fasta

  output:
    path "*.tsv"

  script:
    def protein = getProtein(input_peptides)
    """
    02_count_amino_acids.py ${input_fasta} ${input_peptides} ${protein}.count.tsv
    """
}


// Load the calculated counts and create a plot
process plotCount {  
  input:
    path input_count

  output:
    path "*.png" 

  script:
    def protein = getProtein(input_count)
    """
    03a_plot_count.py ${input_count} ${protein}.plot.png
    """
}


// Get a list of input files from a given folder and create a report
process generateReport {  
  input:
    path input_count

  output:
    path "*.tsv"

  script:
    """
    03b_get_report.py ${input_count} --output_file=protein_report.tsv
    """
}


// Gather result files and archive them
process archiveResults {  
  input:
    path input_plot
    path input_report

  output:
    path "*.tgz"

  script:
    """
    mkdir results
    cp ${input_plot} ${input_report} results/
    tar -czf results.tgz results
    """
}


workflow {
  // Run workflow for all .fasta files in the fasta directory
  fasta = Channel.fromPath("$baseDir/fasta/*.fasta")
  peptides = digestProtein(fasta)
  aa_count = countAA(peptides, fasta)
  count_plot = plotCount(aa_count)
  protein_report = generateReport(aa_count | collect)
  archive_results = archiveResults(count_plot | collect, protein_report)
}
