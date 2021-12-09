#!/usr/bin/env nextflow

// Run workflow for all .fasta files in the fasta directory
fasta = Channel.fromPath("$baseDir/fasta/*.fasta")
fasta.into { 
  fasta_a
  fasta_b 
}


// Helper function to extract the protein name from the filename
def getProtein(fileName) {
  fileName.getBaseName().tokenize(".")[0]
}


// Digest a protein and save the peptides
process digestProtein {
  input:
    path input_fasta from fasta_a

  output:
    path "*.txt" into peptides

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
    path input_peptides from peptides
    path input_fasta from fasta_b

  output:
    path "*.tsv" into aa_count_a, aa_count_b

  script:
    def protein = getProtein(input_peptides)
    """
    02_count_amino_acids.py ${input_fasta} ${input_peptides} ${protein}.count.tsv
    """
}


// Load the calculated counts and create a plot
process plotCount {  
  input:
    path input_count from aa_count_a

  output:
    path "*.png" into count_plot

  script:
    def protein = getProtein(input_count)
    """
    03a_plot_count.py ${input_count} ${protein}.plot.png
    """
}


// Get a list of input files from a given folder and create a report
process generateReport {  
  input:
    path input_count from aa_count_b.collect()

  output:
    path "*.tsv" into protein_report

  script:
    """
    03b_get_report.py ${input_count} --output_file=protein_report.tsv
    """
}


// Gather result files and archive them
process archiveResults {  
  input:
    path input_plot from count_plot.collect()
    path input_report from protein_report

  output:
    path "*.tgz" into archive_results

  script:
    """
    mkdir results
    cp ${input_plot} ${input_report} results/
    tar -czf results.tgz results
    """
}
