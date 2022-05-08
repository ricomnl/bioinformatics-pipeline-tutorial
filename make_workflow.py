"""make_workflow.py"""
import os
from typing import Dict, List, Optional

from redun import task, File
from redun.functools import const


redun_namespace = "bioinformatics_pipeline_tutorial.make_workflow"


COUNT = ["data/KLF4.count.tsv", "data/MYC.count.tsv", "data/PO5F1.count.tsv", "data/SOX2.count.tsv"]
PLOT = ["data/KLF4.plot.png", "data/MYC.plot.png", "data/PO5F1.plot.png", "data/SOX2.plot.png"]


# Custom DSL for describing targets, dependencies (deps), and commands.
rules = {
    "all": {
        "deps": ["data/results.tgz"],
    },
    "clean": {
        "command": "rm -rf data/*",
    },
    "data/%.peptides.txt": {
        "deps": ["fasta/%.fasta"],
        "command": "bin/01_digest_protein.py fasta/%.fasta data/%.peptides.txt"
    },
    "data/%.count.tsv": {
        "deps": ["fasta/%.fasta", "data/%.peptides.txt"],
        "command": "bin/02_count_amino_acids.py fasta/%.fasta data/%.peptides.txt data/%.count.tsv"
    },
    "data/%.plot.png": {
        "deps": ["data/%.count.tsv"],
        "command": "bin/03a_plot_count.py data/%.count.tsv data/%.plot.png"
    },
    "data/protein_report.tsv": {
        "deps": COUNT,
        "command": "bin/03b_get_report.py {COUNT} --output_file=data/protein_report.tsv".format(COUNT=" ".join(COUNT))
    },
    "data/results.tgz": {
        "deps": PLOT + ["data/protein_report.tsv"],
        "command": """rm -rf results
                      mkdir results
                      cp {PLOT} data/protein_report.tsv results/
                      tar -czf data/results.tgz results
                      rm -r results""".format(PLOT=" ".join(PLOT))
    },
}


def match_target(target: str = "all", rules: dict = rules) -> Optional[Dict[str, Dict]]:
    """
    Emulate GNU make pattern matching described here: 
    https://www.gnu.org/software/make/manual/html_node/Pattern-Match.html#Pattern-Match
    """
    rule = rules.get(target)
    if not rule:
        _, tbase = os.path.split(target)
        for rkey, rval in rules.items():
            _, rbase = os.path.split(rkey)
            if not "%" in rbase: continue
            pre, post = rbase.split("%")
            if tbase.startswith(pre) and tbase.endswith(post):
                stem = tbase[len(pre):-len(post)]
                rule = {
                    "deps": [dep.replace("%", stem) for dep in rval.get("deps", [])],
                    "command": rval.get("command", "").replace("%", stem),
                }
                break
    return rule


@task(version="1")
def run_command(command: str, inputs: List[File], output_path: str) -> File:
    """
    Run a shell command to produce a target.
    """
    # Ignore inputs. We pass it as an argument to simply force a dependency.
    assert os.system(command) == 0
    return File(output_path)


@task(version="1")
def make(target: str = "all", rules: dict = rules) -> Optional[File]:
    """
    Make a target (file) using a series of rules.
    """
    rule = match_target(target, rules) if not "%" in target else None
    if not rule:
        # No rule. See if target already exists.
        file = File(target)
        if not file.exists():
            raise ValueError(f"No rule for target: {target}")
        return file

    # Recursively make dependencies.
    inputs = [
        make(dep, rules=rules)
        for dep in rule.get("deps", [])
    ]

    # Run command, if needed.
    if "command" in rule:
        return run_command(rule["command"], inputs, target)
    else:
        return const(None, inputs)
