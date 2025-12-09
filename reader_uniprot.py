#!/usr/bin/env python3

"""
This script takes a FASTA file containing mutant protein sequences plus the
wild-type sequence, extracts amino-acid changes from FASTA headers, and maps
the mutated residue position onto UniProtKB-annotated features via the UniProt
REST API.

Expected FASTA headers (examples):
>SCN10A|rs6795970|V1073A|mut_protein|ENST00000449082
>SCN10A|WT|ENST00000449082|wildtype_protein

Sample usage:
    python3 reader_uniprot_from_fasta.py \
        -m SCN10A_mutated_protein.fasta \
        -g SCN10A \
        -u Q9Y5Y9 \
        -o SCN10A_uniprot_mapping_Q9Y5Y9.tsv

        
Please install the necessary dependencies such as biopython (for processing the fasta files)
by typing 'pip3 install biopython' on terminal.
"""

import argparse
import sys
import requests
from typing import Dict, Any, List, Optional
from Bio import SeqIO

UNIPROT_REST = "https://rest.uniprot.org"


# -------------------------------------------------------------
# UniProt helpers
# -------------------------------------------------------------
def fetch_uniprot_features(uniprot_acc: str) -> List[Dict[str, Any]]:
    """Fetch UniProtKB feature annotations for the given accession."""
    url = f"{UNIPROT_REST}/uniprotkb/{uniprot_acc}.json"
    r = requests.get(url)
    if not r.ok:
        sys.stderr.write(f"[UniProt] Failed to fetch {uniprot_acc}: {r.status_code}\n")
        return []
    return r.json().get("features", [])


def find_uniprot_features_at_position(features: List[Dict[str, Any]], pos: int):
    """Return UniProt feature entries that overlap the given AA position."""
    hits = []
    for feat in features:
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        try:
            s, e = int(start), int(end)
        except ValueError:
            continue
        if s <= pos <= e:
            hits.append(feat)
    return hits


# -------------------------------------------------------------
# FASTA parsing helpers
# -------------------------------------------------------------
def parse_mutation_from_header(header: str):
    """
    Parse header of form:
    SCN10A|rs6795970|V1073A|mut_protein|ENST...
    Returns (rsid, aa_from, pos, aa_to) or None if not a mutant.
    """
    parts = header.split("|")
    if len(parts) < 3:
        return None

    gene = parts[0]
    tag = parts[1]

    # WT entry
    if tag == "WT" or "wild" in header.lower():
        return ("WT", None, None, None)

    # expects something like V1073A
    mut_token = parts[2]
    import re
    m = re.match(r"([A-Z])(\d+)([A-Z])", mut_token)
    if not m:
        return None

    aa_from, pos, aa_to = m.group(1), int(m.group(2)), m.group(3)
    rsid = tag  # rsID or variant name
    return (rsid, aa_from, pos, aa_to)


# -------------------------------------------------------------
# Main
# -------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Map mutant protein FASTA entries to UniProt features.")
    parser.add_argument("-m", "--mutfasta", required=True, help="Mutant+WT protein FASTA file")
    parser.add_argument("-g", "--gene", required=True, help="Gene symbol (e.g., SCN10A)")
    parser.add_argument("-u", "--uniprot", required=True, help="UniProt accession (e.g., Q9Y5Y9)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    # Load UniProt features once
    features = fetch_uniprot_features(args.uniprot)

    out = open(args.output, "w")
    out.write(
        "rsID\tgene\thgvsp\taa_change\tprotein_pos\tuniprot_acc\t"
        "uniprot_feature_type\tuniprot_feature_desc\tuniprot_feature_range\n"
    )

    for rec in SeqIO.parse(args.mutfasta, "fasta"):
        header = rec.id
        parsed = parse_mutation_from_header(header)

        if parsed is None:
            sys.stderr.write(f"[WARN] Could not parse header: {header}\n")
            continue

        rsid, aa_from, pos, aa_to = parsed

        # WT entry → skip
        if rsid == "WT":
            continue

        aa_change_short = f"{aa_from}{pos}{aa_to}"
        hgvsp = f"p.{aa_change_short}"

        sys.stderr.write(f"[INFO] Processing {rsid} at {aa_change_short}\n")

        # Find UniProt features that overlap this position
        hits = find_uniprot_features_at_position(features, pos)

        if not hits:
            out.write(
                f"{rsid}\t{args.gene}\t{hgvsp}\t{aa_change_short}\t{pos}\t"
                f"{args.uniprot}\tNA\tNA\tNA\n"
            )
            continue

        for feat in hits:
            ftype = feat.get("type", "NA")
            fdesc = feat.get("description", "NA")
            loc = feat.get("location", {})
            s = loc.get("start", {}).get("value", "NA")
            e = loc.get("end", {}).get("value", "NA")
            frange = f"{s}-{e}"

            out.write(
                f"{rsid}\t{args.gene}\t{hgvsp}\t{aa_change_short}\t"
                f"{pos}\t{args.uniprot}\t{ftype}\t{fdesc}\t{frange}\n"
            )

    out.close()
    sys.stderr.write(f"[INFO] Wrote UniProt-mapped output to {args.output}\n")


if __name__ == "__main__":
    main()
