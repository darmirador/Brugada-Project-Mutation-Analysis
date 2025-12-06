#!/usr/bin/env python3
"""
Pipeline:
1. Parse mutants+WT FASTA and auto-detect:
   - gene symbol (e.g. SCN10A)
   - reference transcript ID from WT header (e.g. ENST00000449082)
   - WT protein sequence
2. Fetch orthologous protein sequences from Ensembl REST for the gene.
3. Append WT sequence as the human reference.
4. Run MAFFT to build an MSA (optional via --run-mafft).
5. Compute amino-acid conservation scores for mutant positions.

Assumes mutant headers like:
  >SCN10A|rs6795970|V1073A|mut_protein|ENST00000449082
  ...
  >SCN10A|WT|ENST00000449082|wildtype_protein

Usage:
  python3 reader_aminoConservation.py \
  -m SCN5A_mutated_protein.fasta \
  -g SCN5A \
  -o SCN5A_aminoConservation.tsv \
  --run-mafft

Please make sure to download MAFFT via https://mafft.cbrc.jp/alignment/software/macstandard.html
"""

import argparse
import math
import sys
from collections import Counter
from typing import Dict, Iterable, List, Tuple, Optional

import requests
import subprocess

ENSEMBL_REST = "https://rest.ensembl.org"


# ---------------- FASTA utils ----------------

def parse_fasta(path: str) -> Iterable[Tuple[str, str]]:
    """Yield (header, sequence) for a FASTA file (header without '>')."""
    header = None
    seq_chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)


def write_fasta(seqs: Dict[str, str], path: str) -> None:
    with open(path, "w") as f:
        for h, s in seqs.items():
            f.write(f">{h}\n{s}\n")


# ---------------- Header parsing ----------------

def parse_mutant_header(header: str) -> Optional[Dict[str, str]]:
    """
    Parse headers like:
      SCN10A|rs6795970|V1073A|mut_protein|ENST00000449082
      SCN10A|WT|ENST00000449082|wildtype_protein
    """
    parts = header.split("|")
    if len(parts) < 3:
        sys.stderr.write(f"[WARN] Cannot parse header: {header}\n")
        return None

    gene = parts[0]
    second = parts[1]

    if "mut_protein" in parts:
        # gene|rsID|V1073A|mut_protein|ENST...
        rsid = second
        aa_change = parts[2]
        transcript = parts[-1]
        rec_type = "mutant"
    elif "wildtype_protein" in parts:
        # gene|WT|ENST...|wildtype_protein
        rsid = second
        aa_change = None
        transcript = parts[2]
        rec_type = "wildtype"
    else:
        rsid = second
        aa_change = None
        transcript = parts[-1]
        rec_type = "unknown"

    return {
        "gene": gene,
        "rsid": rsid,
        "aa_change": aa_change,
        "transcript": transcript,
        "type": rec_type,
        "raw": header,
    }


def parse_aa_change(aa_change: str) -> Tuple[str, int, str]:
    import re
    m = re.match(r"^([A-Z\*])(\d+)([A-Z\*])$", aa_change)
    if not m:
        raise ValueError(f"Invalid AA change: {aa_change}")
    return m.group(1), int(m.group(2)), m.group(3)


def find_reference_info(mutants_fasta: str) -> Tuple[str, str, str]:
    """
    Scan the mutants FASTA to find the WT record and return:
      (gene_symbol, transcript_id, wt_sequence)

    Expects a header containing 'wildtype_protein' as in:
      SCN10A|WT|ENST00000449082|wildtype_protein
    """
    for header, seq in parse_fasta(mutants_fasta):
        meta = parse_mutant_header(header)
        if meta and meta["type"] == "wildtype":
            gene = meta["gene"]
            transcript = meta["transcript"]
            return gene, transcript, seq
    sys.exit("[FATAL] No wildtype_protein record found in mutants FASTA.")


# ---------------- Ortholog fetching ----------------

def fetch_orthologs(gene_symbol: str) -> Dict[str, str]:
    """
    Fetch protein orthologs for a given human gene symbol using Ensembl REST.
    Returns: dict {header: protein_sequence}
    """
    url = f"{ENSEMBL_REST}/homology/symbol/human/{gene_symbol}"
    params = {"type": "orthologues", "sequence": "protein"}

    print("[DEBUG] Homology URL  =", url)
    print("[DEBUG] Homology params =", params)

    r = requests.get(url, headers={"Content-Type": "application/json"}, params=params)
    print("[DEBUG] Final GET URL =", r.url)
    print("[DEBUG] Status        =", r.status_code)
    print("[DEBUG] Response head =", r.text[:300])

    if not r.ok:
        sys.exit(f"[ERROR] Homology lookup failed: {r.text}")

    data = r.json()
    homs = data["data"][0]["homologies"]

    orthos: Dict[str, str] = {}
    for h in homs:
        tgt = h["target"]
        # Ensembl often uses 'align_seq' for the aligned protein sequence
        seq = (
            tgt.get("align_seq")
            or tgt.get("protein_sequence")
            or tgt.get("sequence")
        )
        if seq is None:
            print("[DEBUG] Skipping target with keys:", list(tgt.keys()))
            continue

        header = f"{tgt['id']}|{tgt['species']}"
        orthos[header] = seq

    print(f"[INFO] Retrieved {len(orthos)} orthologous protein sequences")
    return orthos


# ---------------- MSA and conservation ----------------

def load_msa(path: str) -> Dict[str, str]:
    msa: Dict[str, str] = {}
    for header, seq in parse_fasta(path):
        msa_id = header.split()[0]
        msa[msa_id] = seq
    return msa


def map_aa_to_column(aligned_seq: str, aa_pos: int) -> Optional[int]:
    """
    Map 1-based AA position (ungapped) to 0-based column index in an aligned sequence.
    """
    count = 0
    for idx, aa in enumerate(aligned_seq):
        if aa != "-":
            count += 1
        if count == aa_pos:
            return idx
    return None


def column_conservation(msa: Dict[str, str], col_idx: int) -> Tuple[float, Dict[str, int]]:
    """
    Compute amino-acid conservation at a given MSA column using
    1 - (H / Hmax), where H is Shannon entropy over residue frequencies.
    """
    residues: List[str] = []
    for seq in msa.values():
        if col_idx >= len(seq):
            continue
        aa = seq[col_idx]
        if aa in "-X?":
            continue
        residues.append(aa)

    if not residues:
        return float("nan"), {}

    counts = Counter(residues)
    if len(counts) == 1:
        return 1.0, dict(counts)

    n = sum(counts.values())
    H = 0.0
    for c in counts.values():
        p = c / n
        H -= p * math.log2(p)

    Hmax = math.log2(len(counts))
    score = 1.0 - (H / Hmax) if Hmax > 0 else 0.0
    return score, dict(counts)


# ---------------- Main ----------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute amino-acid conservation for mutant positions using ortholog MSA (auto-detect transcript from FASTA)."
    )
    parser.add_argument("-m", "--mutants", required=True,
                        help="Mutants+WT FASTA with headers like SCN10A|rsID|V1073A|mut_protein|ENST... and SCN10A|WT|ENST...|wildtype_protein")
    parser.add_argument("-g", "--gene", required=True,
                        help="Human gene symbol (e.g. SCN10A) for ortholog fetching.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV file.")
    parser.add_argument("--run-mafft", action="store_true",
                        help="Run MAFFT to build MSA (requires 'mafft' in PATH).")
    args = parser.parse_args()

    # 1) Detect reference info from FASTA
    ref_gene, ref_transcript, wt_seq = find_reference_info(args.mutants)
    print(f"[INFO] Reference from FASTA: gene={ref_gene}, transcript={ref_transcript}")

    # Sanity check that provided -g matches FASTA gene
    if ref_gene != args.gene:
        sys.stderr.write(
            f"[WARN] Gene symbol in FASTA ({ref_gene}) "
            f"differs from -g ({args.gene}); using -g for Ensembl lookup.\n"
        )

    # 2) Fetch orthologs for the gene
    print("[INFO] Fetching orthologs from Ensembl...")
    orthos = fetch_orthologs(args.gene)

    # 3) Build ortholog FASTA and append human WT sequence
    orthologs_fa = f"{args.gene}_orthologs.fasta"
    write_fasta(orthos, orthologs_fa)
    human_header = f"{ref_transcript}|Homo_sapiens"
    with open(orthologs_fa, "a") as f:
        f.write(f">{human_header}\n{wt_seq}\n")
    print(f"[INFO] Wrote orthologs + WT to {orthologs_fa}")

    # 4) Build MSA (optional)
    msa_fa = f"{args.gene}_msa.fasta"
    if args.run_mafft:
        print("[INFO] Running MAFFT...")
        with open(msa_fa, "w") as out_f:
            result = subprocess.run(
                ["mafft", "--auto", orthologs_fa],
                stdout=out_f,
                stderr=subprocess.PIPE,
                text=True,
            )
        if result.returncode != 0:
            sys.stderr.write(f"[ERROR] MAFFT failed:\n{result.stderr}\n")
            sys.exit(1)
        print(f"[INFO] MSA written to {msa_fa}")
    else:
        print("[INFO] --run-mafft not used; treating ortholog FASTA as pre-aligned MSA.")
        msa_fa = orthologs_fa

    msa = load_msa(msa_fa)
    if human_header not in msa:
        sys.exit(f"[FATAL] Human reference {human_header} not found in MSA.")
    ref_aln = msa[human_header]

    # 5) Iterate mutants and compute conservation per site
    seen_rsids = set()
    with open(args.output, "w") as out:
        out.write(
            "gene\trsid\taa_change\tref_aa\tpos\talt_aa\t"
            "ref_transcript\taln_col\tcons_score\tresidue_counts\n"
        )

        for header, seq in parse_fasta(args.mutants):
            meta = parse_mutant_header(header)
            if meta is None or meta["type"] != "mutant":
                continue

            rsid = meta["rsid"]
            if rsid in seen_rsids:
                sys.stderr.write(f"[INFO] Skipping duplicate {rsid}\n")
                continue
            seen_rsids.add(rsid)

            aa_change = meta["aa_change"]
            try:
                ref_aa, pos, alt_aa = parse_aa_change(aa_change)
            except ValueError as e:
                sys.stderr.write(f"[WARN] {e}; skipping {header}\n")
                continue

            col_idx = map_aa_to_column(ref_aln, pos)
            if col_idx is None:
                sys.stderr.write(
                    f"[WARN] Cannot map AA position {pos} to alignment column; skipping {aa_change}\n"
                )
                continue

            score, counts = column_conservation(msa, col_idx)
            counts_str = ";".join(f"{aa}:{cnt}" for aa, cnt in sorted(counts.items()))

            out.write(
                f"{meta['gene']}\t{rsid}\t{aa_change}\t{ref_aa}\t{pos}\t{alt_aa}\t"
                f"{ref_transcript}\t{col_idx+1}\t{score:.4f}\t{counts_str}\n"
            )

    print(f"[INFO] Finished. Results written to {args.output}")


if __name__ == "__main__":
    main()
