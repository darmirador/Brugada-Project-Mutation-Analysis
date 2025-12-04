#!/usr/bin/env python3

"""
This script takes a list of rsIDs for a given human gene, uses Ensembl VEP to
identify amino-acid–changing consequences on that gene’s protein, and then maps
the affected residue positions onto UniProtKB-annotated features (domains,
motifs, PTM sites, etc.) via the UniProt REST API.

Sample input:
    SCN10A_SNP.txt   # one rsID per line; text after ';' is ignored

Sample usage:
    python3 reader_uniprot.py \
        -i SCN10A_SNP.txt \
        -g SCN10A \
        -u Q9Y5Y9 \
        -o SCN10A_uniprot_mapping.tsv

Output:
    SCN10A_uniprot_mapping.tsv  # per-variant mapping to UniProt features
"""

import argparse
import sys
import requests
from typing import Dict, Any, List, Optional

ENSEMBL_VEP = "https://rest.ensembl.org"
UNIPROT_REST = "https://rest.uniprot.org"

AA_CHANGING = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
}

def fetch_uniprot_features(uniprot_acc: str) -> List[Dict[str, Any]]:
    url = f"{UNIPROT_REST}/uniprotkb/{uniprot_acc}.json"
    r = requests.get(url)
    if not r.ok:
        sys.stderr.write(f"[UniProt] Failed to fetch {uniprot_acc}: {r.status_code}\n")
        return []
    return r.json().get("features", [])

def find_uniprot_features_at_position(features: List[Dict[str, Any]], pos: int):
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

def fetch_aa_changing_consequences(rsid: str, gene_symbol: str):
    url = f"{ENSEMBL_VEP}/vep/human/id/{rsid}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    r = requests.get(url, headers=headers)
    if not r.ok:
        sys.stderr.write(f"[VEP] Failed for {rsid}: {r.status_code}\n")
        return []

    data = r.json()
    if not isinstance(data, list) or not data:
        return []

    v = data[0]
    tcs = v.get("transcript_consequences", [])
    gene_tcs = [
        tc for tc in tcs
        if tc.get("gene_symbol") == gene_symbol
        and any(c in AA_CHANGING for c in tc.get("consequence_terms", []))
    ]
    if not gene_tcs:
        return []

    canonical = [tc for tc in gene_tcs if tc.get("canonical") == "YES"]
    return canonical or gene_tcs

def main():
    parser = argparse.ArgumentParser(description="Map rsIDs to UniProt features.")
    parser.add_argument("-i", "--input", required=True, help="Input SNP list file")
    parser.add_argument("-g", "--gene", required=True, help="Gene symbol (e.g., SCN10A)")
    parser.add_argument("-u", "--uniprot", required=True, help="UniProt accession (e.g., Q9Y5Y9)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    features = fetch_uniprot_features(args.uniprot)

    out = open(args.output, "w")
    out.write(
        "rsID\tgene\thgvsp\taa_change\tprotein_pos\tuniprot_acc\t"
        "uniprot_feature_type\tuniprot_feature_desc\tuniprot_feature_range\n"
    )

    seen = set()

    with open(args.input) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            rsid = line.split(";")[0].strip()
            if not rsid:
                continue

            if rsid in seen:
                sys.stderr.write(f"[INFO] Skipping duplicate {rsid}\n")
                continue

            seen.add(rsid)
            sys.stderr.write(f"[INFO] Processing {rsid}\n")

            tcs = fetch_aa_changing_consequences(rsid, args.gene)
            if not tcs:
                out.write(f"{rsid}\t{args.gene}\tNA\tNA\tNA\t{args.uniprot}\tNA\tNA\tNA\n")
                continue

            tc = tcs[0]
            protein_pos = tc.get("protein_start")
            aa = tc.get("amino_acids")
            hgvsp = tc.get("hgvsp")

            if protein_pos is None or aa is None:
                out.write(
                    f"{rsid}\t{args.gene}\t{hgvsp or 'NA'}\tNA\tNA\t"
                    f"{args.uniprot}\tNA\tNA\tNA\n"
                )
                continue

            try:
                pos_int = int(protein_pos)
            except ValueError:
                pos_int = None

            aa_from, aa_to = (aa.split("/") + ["NA", "NA"])[:2]
            aa_change_short = f"{aa_from}{pos_int}{aa_to}" if pos_int else "NA"

            if pos_int is None:
                out.write(
                    f"{rsid}\t{args.gene}\t{hgvsp or 'NA'}\t{aa_change_short}\tNA\t"
                    f"{args.uniprot}\tNA\tNA\tNA\n"
                )
                continue

            hits = find_uniprot_features_at_position(features, pos_int)
            if not hits:
                out.write(
                    f"{rsid}\t{args.gene}\t{hgvsp or 'NA'}\t{aa_change_short}\t"
                    f"{pos_int}\t{args.uniprot}\tNA\tNA\tNA\n"
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
                    f"{rsid}\t{args.gene}\t{hgvsp or 'NA'}\t{aa_change_short}\t"
                    f"{pos_int}\t{args.uniprot}\t{ftype}\t{fdesc}\t{frange}\n"
                )

    out.close()
    sys.stderr.write(f"[INFO] Wrote UniProt-mapped output to {args.output}\n")


if __name__ == "__main__":
    main()