# MPRA Sequence Comparison Tool

This Python tool compares two FASTA files typically used in MPRA (Massively Parallel Reporter Assays) workflows:

- A **reference sequence FASTA file** (original sequences without adaptors)
- An **MPRA sequence FASTA file** (with adaptors attached)

The tool automatically detects adaptor sequences added to the MPRA sequences and generates a detailed comparison report that includes:

- Detected adaptor sequences (forward and reverse)
- Total mismatches and insertions/deletions excluding adaptor regions
- Number of sequences containing adaptor ends
- Summary statistics useful for quality control in MPRA experiments

---

## Features

- Automatic adaptor detection based on longest common prefixes and suffixes
- Accurate mismatch and indel counts, excluding adaptor regions
- Supports usage both as a standalone Python script and inside Jupyter notebooks
- Built using Python and Biopython for reliable sequence handling

---

## Requirements

- Python 3.x
- [Biopython](https://biopython.org/) package

Install Biopython via pip if needed:

```bash
pip install biopython
