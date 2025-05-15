# MPRA Sequence Comparison Tool

This Python tool compares two FASTA files typically used in MPRA (Massively Parallel Reporter Assays) workflows:

- A **reference sequence FASTA file** (original sequences without adaptors)
- An **MPRA sequence FASTA file** (with adaptors attached)

The tool automatically detects adaptor sequences added to the MPRA sequences and reports:

- The detected adaptor sequences (forward and reverse)
- Total mismatches and insertions/deletions excluding adaptor regions
- Number of sequences containing adaptor ends

---

## Features

- Automatic adaptor detection based on longest common prefixes and suffixes
- Detailed mismatch and indel counts excluding adaptors
- Command-line or Jupyter Notebook usage friendly
- Written in Python using Biopython

---

## Requirements

- Python 3.x
- [Biopython](https://biopython.org/) package

Install Biopython via pip if needed:

```bash
pip install biopython
