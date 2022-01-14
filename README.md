# Bioinformatics Pipeline Tutorial

This is the accompanying GitHub repository for this blog post: https://ricomnl.com/blog/bottom-up-bioinformatics-pipeline/.

The workflow we're going to wrap in a pipeline looks like this:
1. Take a set of .fasta protein files
2. Split each into peptides using a variable number of missed cleavages
3. Count the number of cysteines in total as well as the number of peptides that contain a cysteine
4. Generate an output report containing this information in a .tsv file
5. Create an archive to share with colleagues
