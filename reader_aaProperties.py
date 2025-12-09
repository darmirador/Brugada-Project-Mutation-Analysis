#!/usr/bin/env python3
"""
reader_aaProperties.py reads a FASTA file of mutant + WT protein sequences
(e.g. SCN5A_mutated_protein.fasta) and computes per-variant amino-acid
physicochemical properties and changes (charge, polarity, volume, hydropathy,
helix propensity, Grantham distance). It writes a tab-separated table
(e.g. SCN5A_aaProperties.tsv).

Usage:
    python3 reader_aaProperties.py \
        -i SCN5A_mutated_protein.fasta \
        -o SCN5A_aaProperties.tsv
"""

import argparse
import csv
import re
from typing import Dict, Tuple, List

MUT_CODE_RE = re.compile(r"\|([A-Z]\d+[A-Z])\|")

# Kyte–Doolittle hydropathy
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Chou–Fasman helix propensity
HELIX = {
    "A": 1.45, "R": 0.79, "N": 0.73, "D": 1.01, "C": 0.77,
    "Q": 1.17, "E": 1.51, "G": 0.53, "H": 1.00, "I": 1.08,
    "L": 1.34, "K": 1.23, "M": 1.20, "F": 1.12, "P": 0.59,
    "S": 0.79, "T": 0.82, "W": 1.14, "Y": 0.61, "V": 1.06,
}

# Simple charge and polarity
AA_CHARGE = {
    "D": -1, "E": -1,
    "K": +1, "R": +1, "H": +1,
}

AA_POLARITY = {
    "A": "nonpolar", "V": "nonpolar", "L": "nonpolar", "I": "nonpolar",
    "M": "nonpolar", "F": "nonpolar", "W": "nonpolar", "P": "nonpolar",
    "G": "nonpolar",
    "S": "polar", "T": "polar", "Y": "polar", "C": "polar", "N": "polar",
    "Q": "polar",
    "K": "basic", "R": "basic", "H": "basic",
    "D": "acidic", "E": "acidic",
}

# Side-chain volumes (Å^3)
AA_VOLUME = {
    "A":  88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5,
    "Q": 143.8, "E": 138.4, "G":  60.1, "H": 153.2, "I": 166.7,
    "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7,
    "S":  89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0,
}

# Grantham distances (compact table)
_GRANTHAM_TABLE = """
A R 112  A N 111  A D 126  A C 195  A Q 91   A E 107  A G 60   A H 86   A I 94   A L 96
A K 106  A M 84   A F 113  A P 27   A S 99   A T 58   A W 148  A Y 112  A V 64
R N 86   R D 96   R C 180  R Q 43   R E 54   R G 125  R H 29   R I 97   R L 102
R K 26   R M 91   R F 97   R P 103  R S 110  R T 71   R W 101  R Y 77   R V 96
N D 23   N C 139  N Q 46   N E 42   N G 80   N H 68   N I 149  N L 153  N K 94
N M 142  N F 158  N P 91   N S 46   N T 65   N W 174  N Y 143  N V 133
D C 154  D Q 61   D E 45   D G 94   D H 81   D I 168  D L 172  D K 101
D M 160  D F 177  D P 108  D S 65   D T 85   D W 181  D Y 160  D V 152
C Q 154  C E 170  C G 159  C H 174  C I 198  C L 198  C K 202
C M 196  C F 205  C P 169  C S 112  C T 149  C W 215  C Y 194  C V 192
Q E 29   Q G 87   Q H 24   Q I 109  Q L 113  Q K 53   Q M 101  Q F 116
Q P 76   Q S 68   Q T 42   Q W 130  Q Y 99   Q V 96
E G 98   E H 40   E I 134  E L 138  E K 56   E M 126  E F 140
E P 93   E S 80   E T 65   E W 152  E Y 122  E V 121
G H 98   G I 135  G L 138  G K 127  G M 127  G F 153
G P 42   G S 56   G T 59   G W 184  G Y 147  G V 109
H I 94   H L 99   H K 32   H M 87   H F 100
H P 77   H S 89   H T 47   H W 115  H Y 83   H V 84
I L 5    I K 102  I M 10   I F 21   I P 95   I S 142  I T 89
I W 61   I Y 33   I V 29
L K 107  L M 15   L F 22  L P 98   L S 145  L T 92
L W 61   L Y 36   L V 32
K M 95   K F 102  K P 103  K S 121  K T 78
K W 110  K Y 85   K V 97
M F 28   M P 87   M S 135  M T 81
M W 67   M Y 36   M V 21
F P 114  F S 155  F T 103
F W 40   F Y 22   F V 50
P S 74   P T 38   P W 147  P Y 110  P V 68
S T 58   S W 177  S Y 144  S V 124
T W 128  T Y 92   T V 69
W Y 37   W V 88
Y V 55
"""

def build_grantham() -> Dict[Tuple[str, str], int]:
    g: Dict[Tuple[str, str], int] = {}
    tokens = _GRANTHAM_TABLE.split()
    for i in range(0, len(tokens), 3):
        aa1, aa2, val = tokens[i], tokens[i + 1], int(tokens[i + 2])
        g[(aa1, aa2)] = val
        g[(aa2, aa1)] = val
    for aa in KD.keys():
        g[(aa, aa)] = 0
    return g

GRANTHAM = build_grantham()

def parse_fasta(path: str) -> List[Tuple[str, str]]:
    records = []
    header = None
    seq_chunks: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            records.append((header, "".join(seq_chunks)))
    return records

def main():
    parser = argparse.ArgumentParser(
        description="Compute physicochemical changes for mutant proteins vs WT."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input FASTA with mutant proteins + WT"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file (e.g. SCN5A_aaProperties.tsv)"
    )
    args = parser.parse_args()

    records = parse_fasta(args.input)

    wt_seq = None
    wt_header = None
    for h, s in records:
        if "wildtype_protein" in h:
            wt_seq = s
            wt_header = h
            break
    if wt_seq is None:
        raise SystemExit("WT sequence (wildtype_protein) not found in FASTA.")

    with open(args.output, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow([
            "gene", "rsid", "mutation", "position",
            "wt_aa", "mut_aa",
            "charge_wt", "charge_mut", "delta_charge",
            "polarity_wt", "polarity_mut",
            "volume_wt", "volume_mut", "delta_volume",
            "kd_wt", "kd_mut", "delta_kd",
            "helix_wt", "helix_mut", "delta_helix",
            "grantham"
        ])

        for header, seq in records:
            if header == wt_header:
                continue

            m = MUT_CODE_RE.search(header)
            if not m:
                continue
            mut_code = m.group(1)  # e.g. H558R
            wt_aa = mut_code[0]
            mut_aa = mut_code[-1]
            pos = int(mut_code[1:-1])

            if wt_seq[pos - 1] != wt_aa:
                print(
                    f"[WARN] WT residue mismatch at {mut_code} "
                    f"(WT has {wt_seq[pos - 1]})."
                )

            gene, rsid = "", ""
            parts = header.split("|")
            if len(parts) >= 1:
                gene = parts[0]
            if len(parts) >= 2:
                rsid = parts[1]

            charge_wt = AA_CHARGE.get(wt_aa, 0)
            charge_mut = AA_CHARGE.get(mut_aa, 0)
            delta_charge = charge_mut - charge_wt

            pol_wt = AA_POLARITY.get(wt_aa, "unknown")
            pol_mut = AA_POLARITY.get(mut_aa, "unknown")

            vol_wt = AA_VOLUME[wt_aa]
            vol_mut = AA_VOLUME[mut_aa]
            delta_vol = vol_mut - vol_wt

            kd_wt = KD[wt_aa]
            kd_mut = KD[mut_aa]
            delta_kd = kd_mut - kd_wt

            helix_wt = HELIX[wt_aa]
            helix_mut = HELIX[mut_aa]
            delta_helix = helix_mut - helix_wt

            grantham = GRANTHAM.get((wt_aa, mut_aa), -1)

            writer.writerow([
                gene, rsid, mut_code, pos,
                wt_aa, mut_aa,
                charge_wt, charge_mut, delta_charge,
                pol_wt, pol_mut,
                vol_wt, vol_mut, delta_vol,
                kd_wt, kd_mut, delta_kd,
                helix_wt, helix_mut, delta_helix,
                grantham,
            ])

if __name__ == "__main__":
    main()
