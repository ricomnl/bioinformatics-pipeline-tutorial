# [Extension: redun] Bioinformatics Pipeline Tutorial

This is the accompanying GitHub repository for this blog post: https://ricomnl.com/blog/bottom-up-bioinformatics-pipeline/.

## Getting Started
Install requirements:
```bash
pip install -r requirements.txt
```

Launch a workflow:
```bash
# Main workflow
redun run workflow.py main --input-dir fasta/

# Script workflow
redun run script_workflow.py main --input-dir fasta/

# Make workflow
redun run make_workflow.py make --target all
```
