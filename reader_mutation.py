#!/usr/bin/env python3

"""
This script takes a list of rsIDs for a given human gene, uses Ensembl REST
endpoints (lookup/sequence/VEP) to identify coding consequences on the
canonical protein-coding transcript, mutates the CDS/cDNA/protein sequences
accordingly, and summarizes the effects.

Sample input:
    SCN5A_SNP.txt   # one rsID per line; text after ';' is ignored

Sample usage:
    python3 reader_mutation.py \
        -g SCN5A \
        -i SCN5A_SNP.txt \
        -c SCN5A_mutated_cdna.fasta \
        -p SCN5A_mutated_protein.fasta

Outputs:
    SCN5A_mutated_cdna.fasta          (mutated + WT cDNA sequences)
    SCN5A_mutated_protein.fasta       (mutated + WT protein sequences)
    SCN5A_variant_classification.tsv  (per-variant consequence summary)
"""

import argparse
import sys
from typing import List, Dict, Tuple, Optional
import requests

ENSEMBL_REST = "https://rest.ensembl.org"

# ---------------- Ensembl helpers ----------------

def get_canonical_transcript(gene_symbol: str) -> Tuple[str, str]:
    url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens/{gene_symbol}"
    params = {"expand": 1}
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, params=params, headers=headers)
    if not r.ok:
        sys.exit(f"Error fetching gene {gene_symbol}: {r.text}")

    data = r.json()
    transcripts = data.get("Transcript", [])
    if not transcripts:
        sys.exit(f"No transcripts found for {gene_symbol}")

    for t in transcripts:
        if t.get("is_canonical") == 1 and t.get("biotype") == "protein_coding":
            print(f"[INFO] Using canonical transcript {t['id']} for {gene_symbol}")
            return t["id"], data["id"]

    pcs = [t for t in transcripts if t.get("biotype") == "protein_coding"]
    if pcs:
        pcs.sort(key=lambda t: t.get("length", 0), reverse=True)
        return pcs[0]["id"], data["id"]

    sys.exit("No usable protein-coding transcript found.")


def get_sequence(id_: str, seq_type: str) -> str:
    if id_ is None:
        sys.exit(f"No ID provided for {seq_type} sequence.")
    url = f"{ENSEMBL_REST}/sequence/id/{id_}"
    params = {"type": seq_type}
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, params=params, headers=headers)
    if not r.ok:
        sys.exit(f"Error fetching {seq_type} for {id_}: {r.text}")
    return r.json().get("seq").upper()


def get_translation_id(tid: str) -> Optional[str]:
    url = f"{ENSEMBL_REST}/lookup/id/{tid}"
    params = {"expand": 1}
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, params=params, headers=headers)
    if not r.ok:
        return None
    trans = r.json().get("Translation")
    return trans.get("id") if trans else None


def query_vep(ids: List[str]) -> List[Dict]:
    """
    GET /vep/human/id/:id without conservation plugins.
    """
    vid = ids[0]
    url = f"{ENSEMBL_REST}/vep/human/id/{vid}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    r = requests.get(url, headers=headers)
    if not r.ok:
        return []
    return r.json()


def pick_tc(records: List[Dict], transcript_id: str, gene: str) -> Optional[Dict]:
    for rec in records:
        for tc in rec.get("transcript_consequences", []):
            if tc.get("transcript_id") == transcript_id:
                return tc
    for rec in records:
        for tc in rec.get("transcript_consequences", []):
            if tc.get("gene_symbol") == gene:
                return tc
    return None

# ---------------- Mutators ----------------

def mutate_cds(cds: str, protein_pos: int, codons: str) -> Optional[str]:
    if protein_pos is None or "/" not in codons:
        return None
    ref_codon, alt_codon = codons.upper().split("/")
    idx0 = (protein_pos - 1) * 3
    if idx0 + 3 > len(cds) or cds[idx0:idx0+3] != ref_codon:
        return None
    return cds[:idx0] + alt_codon + cds[idx0+3:]


def embed_cds(cdna: str, ref_cds: str, new_cds: str) -> Optional[str]:
    pos = cdna.find(ref_cds)
    if pos == -1:
        return None
    return cdna[:pos] + new_cds + cdna[pos+len(ref_cds):]


def mutate_protein(prot: str, protein_pos: int, aa_str: str) -> Optional[str]:
    if protein_pos is None or "/" not in aa_str:
        return None
    ref, alt = aa_str.split("/")
    idx = protein_pos - 1
    if idx < 0 or idx >= len(prot) or prot[idx] != ref.upper():
        return None
    return prot[:idx] + alt.upper() + prot[idx+1:]


# ---------------- I/O ----------------

def read_variants(path: str) -> List[str]:
    out = []
    with open(path) as f:
        for line in f:
            x = line.strip()
            if x:
                out.append(x.split(";")[0])
    return out


def write_fasta(entries, path):
    with open(path, "w") as f:
        for hdr, seq in entries:
            f.write(f">{hdr}\n{seq}\n")


def write_tsv(rows, path):
    with open(path, "w") as f:
        f.write("variant\tgene\tcategory\tamino_acid_change\tnote\n")
        for r in rows:
            f.write(
                f"{r['variant']}\t{r['gene']}\t{r['category']}\t"
                f"{r['aa_change']}\t{r['note']}\n"
            )


# ---------------- MAIN ----------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-g", "--gene", required=True)
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-c", "--cdna_out", default="mutated_cdna.fasta")
    p.add_argument("-p", "--protein_out", default="mutated_protein.fasta")
    args = p.parse_args()

    gene = args.gene
    variants = read_variants(args.input)

    transcript_id, gene_id = get_canonical_transcript(gene)
    cdna = get_sequence(transcript_id, "cdna")
    cds = get_sequence(transcript_id, "cds")
    protein_id = get_translation_id(transcript_id)
    prot = get_sequence(protein_id, "protein")

    cdna_fastas, prot_fastas, rows = [], [], []
    coding_terms = {"missense_variant", "stop_gained", "synonymous_variant"}

    seen = set()

    for rs in variants:
        if rs in seen:
            print(f"[INFO] Skipping duplicate {rs}")
            continue
        seen.add(rs)

        print(f"[INFO] Processing {rs}")
        records = query_vep([rs])

        if not records:
            rows.append({
                "variant": rs, "gene": gene,
                "category": "NA", "aa_change": "NA",
                "note": "does not map to CDS"
            })
            continue

        rec0 = records[0]
        tc = pick_tc(records, transcript_id, gene)
        top_terms = rec0.get("consequence_terms", [])

        if tc is None:
            rows.append({
                "variant": rs, "gene": gene,
                "category": ",".join(top_terms),
                "aa_change": "NA",
                "note": "does not map to CDS"
            })
            continue

        terms = tc.get("consequence_terms", [])
        protein_pos = tc.get("protein_start")
        aa_str = tc.get("amino_acids", "")
        codons = tc.get("codons", "")

        aa_change = "NA"
        dna_change = "NA"

        if "/" in aa_str and protein_pos:
            raa, aaa = aa_str.split("/")
            aa_change = f"{raa}{protein_pos}{aaa}"

        if "/" in codons and protein_pos:
            rc, ac = codons.upper().split("/")
            codon_start = (protein_pos - 1) * 3 + 1
            codon_end = codon_start + 2
            dna_change = f"c.{codon_start}_{codon_end} {rc}>{ac}"

        if not any(t in coding_terms for t in terms):
            rows.append({
                "variant": rs, "gene": gene,
                "category": ",".join(terms),
                "aa_change": aa_change,
                "note": "does not map to CDS"
            })
            continue

        # Attempt CDS mutation
        new_cds = mutate_cds(cds, protein_pos, codons)
        if new_cds is None:
            rows.append({
                "variant": rs, "gene": gene,
                "category": ",".join(terms),
                "aa_change": aa_change,
                "note": "no change"
            })
            continue

        # Mutate cDNA & protein
        new_cdna = embed_cds(cdna, cds, new_cds) or new_cds
        new_prot = mutate_protein(prot, protein_pos, aa_str) or prot

        cdna_changed = new_cdna != cdna
        prot_changed = new_prot != prot

        if not cdna_changed and not prot_changed:
            note = "no change"

        elif cdna_changed and not prot_changed:
            note = f"cDNA changed only - {dna_change}"
            cdna_fastas.append((f"{gene}|{rs}|{aa_change}|mut_cDNA|{transcript_id}", new_cdna))

        elif cdna_changed and prot_changed:
            note = f"cDNA and protein changed - DNA:{dna_change}; protein:{aa_change}"
            cdna_fastas.append((f"{gene}|{rs}|{aa_change}|mut_cDNA|{transcript_id}", new_cdna))
            prot_fastas.append((f"{gene}|{rs}|{aa_change}|mut_protein|{transcript_id}", new_prot))

        else:
            note = f"protein changed only - {aa_change}"
            prot_fastas.append((f"{gene}|{rs}|{aa_change}|mut_protein|{transcript_id}", new_prot))

        rows.append({
            "variant": rs, "gene": gene,
            "category": ",".join(terms),
            "aa_change": aa_change,
            "note": note
        })

    # Append WT once
    cdna_fastas.append((f"{gene}|WT|{transcript_id}|wildtype_cDNA", cdna))
    prot_fastas.append((f"{gene}|WT|{transcript_id}|wildtype_protein", prot))

    write_fasta(cdna_fastas, args.cdna_out)
    write_fasta(prot_fastas, args.protein_out)
    write_tsv(rows, f"{gene}_variant_classification.tsv")

    print("[INFO] Completed successfully.")

if __name__ == "__main__":
    main()