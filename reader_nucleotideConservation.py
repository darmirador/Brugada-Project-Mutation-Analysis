#!/usr/bin/env python3
# TO DOWNLOAD .bw file

"""
reader_phyloP.py

This script takes a list of human SNP rsIDs, resolves their GRCh38 genomic
coordinates using the Ensembl REST API, and then retrieves phyloP100way
conservation scores from a local bigWig file via pyBigWig. Results are written
to a tab-delimited (.tsv) output file.

Sample usage:
    python3 reader_phyloP.py \
        -i SCN10A_SNP.txt \
        -b hg38.phyloP100way.bw \
        -o SCN10A_phyloP.tsv

Arguments:
    -i / --input   : Text file with rsIDs (rsID at start of each line; text after ';' ignored)
    -b / --bigwig  : Path to hg38.phyloP100way.bw (PLEASE MAKE SURE TO DOWNLOAD THIS VIA UCSC
    AND PUT IN THE SAME WORKING DIRECTORY. See 'hg38.phyloP100way.bw'. It is intentionally excluded from
    version control via `.gitignore`, so it will not appear in the Git repository.)
    -o / --output  : Output TSV file containing rsID, genomic_location, and phyloP100way

Please install the necessary dependencies such as pyBigWig (for reading the .bw file) 
by typing 'pip3 install pyBigWig' on terminal.
"""

import sys
import math
import argparse
import requests
import pyBigWig

ENSEMBL_REST = "https://rest.ensembl.org"


def get_grch38_coords_from_rs(rsid: str):
    """Return (chrom, pos) in GRCh38 (UCSC style chrN) using Ensembl REST."""
    url = f"{ENSEMBL_REST}/variation/human/{rsid}"
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, headers=headers)
    if not r.ok:
        sys.stderr.write(f"[Ensembl] Failed for {rsid}: {r.status_code} {r.text}\n")
        return None

    data = r.json()
    mappings = data.get("mappings", [])
    if not mappings:
        sys.stderr.write(f"[Ensembl] No mappings for {rsid}\n")
        return None

    for m in mappings:
        asm = (m.get("assembly_name") or "")
        if asm.startswith("GRCh38"):
            chrom = "chr" + str(m["seq_region_name"])
            pos = int(m["start"])  # 1-based
            return chrom, pos

    sys.stderr.write(f"[Ensembl] No GRCh38 mapping for {rsid}\n")
    return None


def get_phylop_from_bigwig(bw, chrom: str, pos_1based: int):
    """Return phyloP score from bigWig at 1-based genomic position, or None."""
    start0 = pos_1based - 1  # bigWig is 0-based, half-open
    end = pos_1based
    try:
        vals = bw.values(chrom, start0, end)
    except RuntimeError as e:
        sys.stderr.write(f"[pyBigWig] Error at {chrom}:{pos_1based}: {e}\n")
        return None

    if not vals:
        return None
    v = vals[0]
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return None
    return float(v)


def main():
    parser = argparse.ArgumentParser(
        description="Map rsIDs to GRCh38 coordinates and phyloP100way scores from a local bigWig."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input SNP list file (rsID at start of each line; text after ';' ignored)"
    )
    parser.add_argument(
        "-b", "--bigwig", required=True,
        help="Path to hg38.100way.phyloP100way.bw (or other compatible phyloP bigWig)"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file"
    )
    args = parser.parse_args()

    snp_file = args.input
    bw_path = args.bigwig
    out_path = args.output

    bw = pyBigWig.open(bw_path)

    with open(out_path, "w") as out:
        out.write("rsID\tgenomic_location\tphyloP100way\n")

        seen = set()
        with open(snp_file) as f:
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

                coords = get_grch38_coords_from_rs(rsid)
                if coords is None:
                    out.write(f"{rsid}\tNA\tNA\n")
                    continue

                chrom, pos = coords
                loc_str = f"{chrom}:{pos}"

                score = get_phylop_from_bigwig(bw, chrom, pos)
                if score is None:
                    out.write(f"{rsid}\t{loc_str}\tNA\n")
                else:
                    out.write(f"{rsid}\t{loc_str}\t{score:.6f}\n")

    bw.close()
    sys.stderr.write(f"[INFO] Wrote phyloP output to {out_path}\n")


if __name__ == "__main__":
    main()
