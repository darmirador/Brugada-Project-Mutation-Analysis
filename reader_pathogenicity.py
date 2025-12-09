#!/usr/bin/env python3

"""
Description:
    This script annotates missense mutants from a FASTA file and summarizes
    all pathogenicity-related predictions into a single TSV output. It extracts
    rsIDs and amino-acid changes from FASTA headers, queries NCBI dbSNP/ClinVar
    for clinical significance annotations, and optionally queries Ensembl VEP
    for variant effect predictions (impact, consequence terms, SIFT, PolyPhen).
    Wild-type entries are automatically skipped (the last entry in this workflow). 
    Each mutant is written as a single row in the TSV, with separate columns corresponding 
    to each algorithm.

Input:
    A FASTA file containing wild-type + mutated protein sequences with headers in
    the format:
        >GENE|rsID|AA_CHANGE|mut_protein|TRANSCRIPT
        >GENE|WT|TRANSCRIPT|wildtype_protein

    Example:
        >SCN5A|rs1805124|H558R|mut_protein|ENST00000423572
        MANFLLP...

Output:
    A single TSV file summarizing all algorithms:
        gene
        aa_change
        rsid
        clinvar_class
        clinvar_score
        vep_impact
        vep_consequence_terms
        vep_sift
        vep_polyphen
        clinvar_counts_json

Usage:
    python3 reader_pathogenicity.py \
        -i SCN5A_mutated_protein.fasta \
        -o SCN5A_pathogenicityScore.tsv

Notes:
    - Requires internet access for NCBI and VEP REST API queries.
    - Delay between HTTP calls can be adjusted via --delay (default: 1.0 sec).
    - AlphaMissense and EVE integration removed as requested.
    - This script only produces ONE TSV output containing all available annotations.
"""


import json
import sys
import time
import argparse
from typing import Dict, List, Tuple, Optional
import requests

NCBI_BASE = "https://api.ncbi.nlm.nih.gov/variation/v0"
VEP_BASE = "https://rest.ensembl.org"

SIGNIF_SCORE = {
    "pathogenic": 1.0,
    "likely_pathogenic": 0.8,
    "uncertain_significance": 0.5,
    "likely_benign": 0.2,
    "benign": 0.0,
}

def parse_fasta(path):
    records = []
    header = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line: continue
            if line.startswith(">"):
                if header: records.append((header, seq))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header:
        records.append((header, seq))
    return records

def parse_header(header):
    parts = header.split("|")
    info = {}
    if len(parts) == 5:
        gene, rsid_tag, aa, kind, transcript = parts
        info["gene"] = gene
        info["kind"] = kind
        info["aa"] = aa
        info["rsid"] = rsid_tag if rsid_tag.startswith("rs") else ""
    elif len(parts) == 4:
        gene, rsid_tag, transcript, kind = parts
        info["gene"] = gene
        info["kind"] = kind
        info["aa"] = ""
        info["rsid"] = rsid_tag if rsid_tag.startswith("rs") else ""
    return info

def fetch_refsnp(rsid):
    url = f"{NCBI_BASE}/refsnp/{rsid[2:]}"
    r = requests.get(url, timeout=30)
    if not r.ok: return None
    return r.json()

def extract_clinvar_significance(data):
    sigs = []
    psd = data.get("primary_snapshot_data", {})
    for ann in psd.get("allele_annotations", []):
        for c in ann.get("clinical", []):
            for s in c.get("clinical_significances", []):
                if isinstance(s, str): sigs.append(s.strip())
    return sigs

def normalize(label):
    l = label.lower().replace(" ", "_")
    if l in SIGNIF_SCORE: return l
    if "likely" in l and "pathogenic" in l: return "likely_pathogenic"
    if "likely" in l and "benign" in l: return "likely_benign"
    if "uncertain" in l: return "uncertain_significance"
    if "pathogenic" == l: return "pathogenic"
    if "benign" == l: return "benign"
    return None

def summarize_clinvar(sig_list):
    counts = {}
    scores = []
    for s in sig_list:
        norm = normalize(s)
        if not norm: continue
        counts[norm] = counts.get(norm, 0) + 1
        scores.append(SIGNIF_SCORE[norm])
    if not scores:
        return "no_data", "NA", counts
    mean = sum(scores)/len(scores)
    pred = "pathogenic" if mean>0.5 else "benign" if mean<0.5 else "uncertain"
    return pred, f"{mean:.3f}", counts

def fetch_vep(rsid):
    url = f"{VEP_BASE}/vep/human/id/{rsid}"
    headers = {"Accept": "application/json"}
    r = requests.get(url, headers=headers, timeout=30)
    if not r.ok: return None
    data = r.json()
    return data[0] if isinstance(data, list) and data else None

def summarize_vep(data):
    if not data: return "NA","NA","NA","NA"
    cons = data.get("transcript_consequences") or []
    if not cons: return "NA","NA","NA","NA"
    c = cons[0]
    impact = c.get("impact","NA")
    terms = ",".join(c.get("consequence_terms",[])) or "NA"
    sift = c.get("sift_prediction","NA")
    poly = c.get("polyphen_prediction","NA")
    return impact, terms, sift, poly

def annotate(input_fasta, output_tsv, delay):
    records = parse_fasta(input_fasta)
    with open(output_tsv,"w") as out:
        out.write("gene\taa_change\trsid\tclinvar_class\tclinvar_score\t"
                  "vep_impact\tvep_consequences\tvep_sift\tvep_polyphen\tclinvar_counts\n")
        for header, seq in records:
            info = parse_header(header)
            if info.get("kind","").lower().startswith("wild"):
                continue

            gene = info.get("gene","")
            aa = info.get("aa","")
            rsid = info.get("rsid","")

            # --- ClinVar ---
            if rsid:
                data = fetch_refsnp(rsid)
                sigs = extract_clinvar_significance(data) if data else []
                clin_class, clin_score, counts = summarize_clinvar(sigs)
                time.sleep(delay)
            else:
                clin_class, clin_score, counts = "no_rsid","NA",{}

            # --- VEP ---
            if rsid:
                v = fetch_vep(rsid)
                impact, cons, sift, poly = summarize_vep(v)
                time.sleep(delay)
            else:
                impact=cons=sift=poly="NA"

            out.write(
                f"{gene}\t{aa}\t{rsid or 'NA'}\t"
                f"{clin_class}\t{clin_score}\t"
                f"{impact}\t{cons}\t{sift}\t{poly}\t"
                f"{json.dumps(counts)}\n"
            )

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i","--input",required=True)
    p.add_argument("-o","--output",required=True)
    p.add_argument("--delay",type=float,default=1.0)
    args=p.parse_args()
    annotate(args.input,args.output,args.delay)

if __name__=="__main__":
    main()