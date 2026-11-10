#!/usr/bin/env python3
"""
Build the UV hotspot position table for the GMS-Solid panel (hg38).

Combines the UV sites in the panel design (GMS_Solid_target_annotated_hg38.UV.bed) with
the ranked list of recurrent ETS promoter hotspots in melanoma
(recurrent_ETS_top500_dbsnp_GMS.xlsx) and writes one row per position, hg38 1-based:

  Y1, Y2  cytosines in the YY dipyrimidine of the YYTTCCG motif, the main hotspot
  C5      the cytosine after TT, a weaker hotspot (not checked against dbSNP)
  flank   every C/G in the +/-15 bp window around YY outside the septamer, used as
          an in-sample background for C>T errors

ref/alt are on the reference strand: for motifs on the minus strand the C>T change is
G>A on the reference. CpG is True for a C followed by G, or a G preceded by C.

Usage:
  python3 make_uv_hotspot_positions.py \
      --bed GMS_Solid_target_annotated_hg38.UV.bed \
      --xlsx recurrent_ETS_top500_dbsnp_GMS.xlsx \
      --out GMS_Solid_UV_hotspot_positions_20261110_hg38.tsv
"""

import argparse
import csv

import pandas as pd

COMPLEMENT = str.maketrans("ACGT", "TGCA")
CT_ALT = {"C": "T", "G": "A"}
COLUMNS = ["site", "gene", "rank", "motif", "motif_strand", "cls", "chr", "pos", "ref", "alt", "CpG"]


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def site_positions(site, r):
    """
    Position rows for one site. r is a row of the hotspot list, where sequence[i] is
    1-based position start_adj + i and region_start/region_end are 0-based half-open.
    """
    offset = r.start - r.start_adj
    septamer = r.sequence[offset:offset + 7]
    if septamer != r.Septamer:
        raise ValueError(f"{site}: septamer {r.Septamer} does not match the sequence ({septamer})")
    plus = r.plusstr == 1
    motif = septamer if plus else reverse_complement(septamer)
    # Genomic position of motif positions 1-7, read 5'->3' on the pyrimidine strand
    motif_pos = [r.start + i for i in range(7)] if plus else [r.end - i for i in range(7)]

    def is_cpg(pos, ref):
        i = pos - r.start_adj
        if ref == "C":
            return i + 1 < len(r.sequence) and r.sequence[i + 1] == "G"
        return i > 0 and r.sequence[i - 1] == "C"

    def row(cls, pos, ref):
        return {
            "site": site, "gene": r.gene, "rank": r["rank"], "motif": motif, "motif_strand": "+" if plus else "-",
            "cls": cls, "chr": r.chr, "pos": pos, "ref": ref, "alt": CT_ALT[ref], "CpG": is_cpg(pos, ref),
        }

    rows = []
    for index, cls in ((0, "Y1"), (1, "Y2"), (4, "C5")):
        if motif[index] == "C":
            rows.append(row(cls, motif_pos[index], "C" if plus else "G"))

    septamer_span = range(r.start, r.end + 1)
    for pos in range(r.region_start + 1, r.region_end + 1):
        ref = r.sequence[pos - r.start_adj]
        if pos not in septamer_span and ref in CT_ALT:
            rows.append(row("flank", pos, ref))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bed", required=True, help="UV sites in the panel design (septamer coordinates)")
    parser.add_argument("--xlsx", required=True, help="ranked recurrent ETS promoter hotspot list")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    hotspots = pd.read_excel(args.xlsx)
    hotspots["rank"] = range(1, len(hotspots) + 1)
    bed = pd.read_csv(args.bed, sep="\t", header=None, names=["chr", "bed_start", "bed_end", "site"])
    merged = bed.merge(hotspots, left_on=["chr", "bed_start"], right_on=["chr", hotspots["start"] - 1], how="left")
    missing = merged[merged["gene"].isna()]
    if len(missing):
        raise ValueError(f"BED sites not in the hotspot list: {', '.join(missing.site)}")

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for _, r in merged.iterrows():
            writer.writerows(site_positions(r.site, r))


if __name__ == "__main__":
    main()
