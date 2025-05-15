"""
MPRA Sequence Comparison Tool

Author:
Prof. Amir Ali Abbasi
National Center for Bioinformatics, Quaid-i-Azam University Islamabad
Email: abbasiam@qau.edu.pk
"""

from Bio import SeqIO
import re
import os

def clean_sequence(seq):
    """
    Remove bracketed ambiguous bases like [A], [C], etc. from sequences.
    """
    return re.sub(r'\[([ACGT])\]', r'\1', seq)

def longest_common_prefix(strs):
    """
    Find the longest common prefix among a list of strings.
    """
    if not strs:
        return ""
    prefix = strs[0]
    for s in strs[1:]:
        i = 0
        while i < len(prefix) and i < len(s) and prefix[i] == s[i]:
            i += 1
        prefix = prefix[:i]
        if prefix == "":
            break
    return prefix

def longest_common_suffix(strs):
    """
    Find the longest common suffix among a list of strings.
    """
    if not strs:
        return ""
    suffix = strs[0]
    for s in strs[1:]:
        i = 0
        while i < len(suffix) and i < len(s) and suffix[-1 - i] == s[-1 - i]:
            i += 1
        suffix = suffix[-i:] if i > 0 else ""
        if suffix == "":
            break
    return suffix

def detect_and_strip_adaptors_auto(seqs, min_adaptor_len=5):
    """
    Automatically detect common forward and reverse adaptor sequences
    from a list of sequences, then strip them.

    Returns:
        dict with:
            - forward_adaptor (str)
            - reverse_adaptor (str)
            - forward_detected_count (int)
            - reverse_detected_count (int)
            - stripped_sequences (list of str)
    """
    fwd_adaptor = longest_common_prefix(seqs)
    rev_adaptor = longest_common_suffix(seqs)

    if len(fwd_adaptor) < min_adaptor_len:
        fwd_adaptor = ""
    if len(rev_adaptor) < min_adaptor_len:
        rev_adaptor = ""

    stripped_seqs = []
    fwd_detected_count = 0
    rev_detected_count = 0

    for seq in seqs:
        fwd_found = seq.startswith(fwd_adaptor) if fwd_adaptor else False
        rev_found = seq.endswith(rev_adaptor) if rev_adaptor else False

        if fwd_found:
            seq = seq[len(fwd_adaptor):]
            fwd_detected_count += 1
        if rev_found:
            seq = seq[:-len(rev_adaptor)]
            rev_detected_count += 1
        stripped_seqs.append(seq)

    return {
        "forward_adaptor": fwd_adaptor,
        "reverse_adaptor": rev_adaptor,
        "forward_detected_count": fwd_detected_count,
        "reverse_detected_count": rev_detected_count,
        "stripped_sequences": stripped_seqs,
    }

def compare_mpra_sequences(reference_fasta, mpra_fasta, output_filename="mpra_comparison_report.txt"):
    """
    Compare two FASTA files containing reference sequences and MPRA sequences,
    detect adaptor sequences automatically, and report differences.

    Output:
        Prints detailed comparison report and saves to a text file.
    """

    recs_ref = list(SeqIO.parse(reference_fasta, "fasta"))
    recs_mpra = list(SeqIO.parse(mpra_fasta, "fasta"))

    if len(recs_ref) != len(recs_mpra):
        print(f"❌ Sequence count mismatch: {len(recs_ref)} vs {len(recs_mpra)}")
        return

    total_sequences = len(recs_ref)

    # Detect adaptors from MPRA file sequences automatically
    mpra_seqs_raw = [str(rec.seq) for rec in recs_mpra]
    adaptor_info = detect_and_strip_adaptors_auto(mpra_seqs_raw)

    fwd_adaptor = adaptor_info["forward_adaptor"]
    rev_adaptor = adaptor_info["reverse_adaptor"]
    fwd_count = adaptor_info["forward_detected_count"]
    rev_count = adaptor_info["reverse_detected_count"]
    mpra_stripped = adaptor_info["stripped_sequences"]

    total_mismatches = 0
    total_indels_excl_adaptors = 0
    total_adaptor_sequences = fwd_count + rev_count

    for i, rec_ref in enumerate(recs_ref):
        ref_seq = clean_sequence(str(rec_ref.seq))
        mpra_seq_trimmed = mpra_stripped[i]

        len_ref = len(ref_seq)
        len_mpra = len(mpra_seq_trimmed)

        total_indels_excl_adaptors += abs(len_ref - len_mpra)

        # Count mismatches in overlapping length
        min_len = min(len_ref, len_mpra)
        mismatches = sum(1 for j in range(min_len) if ref_seq[j] != mpra_seq_trimmed[j])
        total_mismatches += mismatches

    report_lines = [
        "=== MPRA Sequence Comparison Report ===",
        "",
        f'Total sequences compared: {total_sequences}',
        f'Adaptor sequences detected in MPRA file: {"Yes" if total_adaptor_sequences > 0 else "No"}',
        f'Forward adaptor sequence: {fwd_adaptor if fwd_adaptor else "None"}',
        f'Reverse adaptor sequence: {rev_adaptor if rev_adaptor else "None"}',
        f'Total mismatches (excluding adaptors): {total_mismatches}',
        f'Total insertions/deletions (excluding adaptors): {total_indels_excl_adaptors}',
        f'Total adaptor sequence counts detected: {total_adaptor_sequences} (counts as number of sequences with adaptors)',
    ]

    report_text = "\n".join(report_lines)

    print(report_text)

    with open(output_filename, "w") as f:
        f.write(report_text)

    print(f"\n✅ Report saved as: {os.path.abspath(output_filename)}")
