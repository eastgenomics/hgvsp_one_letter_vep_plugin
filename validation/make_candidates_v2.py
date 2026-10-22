#!/usr/bin/env python3
"""Generate synthetic GRCh38 variants covering every HGVSp pattern VEP emits.

Uses only the public Ensembl REST API (gene models + reference sequence), no
patient data. For every coding base, plus a few bases either side of each
coding exon edge (start/stop codons, splice boundaries), it writes:
  - all three SNVs
  - 2 and 3 bp substitutions (delins)
  - 1, 2, 3, 4 and 6 bp deletions
  - 1 bp insertions of each base, 3 and 6 bp random insertions
  - in-frame insertion of each stop codon
  - 3 and 6 bp duplications

Usage: python3 make_candidates_v2.py <out.vcf> GENE [GENE ...]
"""
import json
import random
import sys
import urllib.request

REST = "https://rest.ensembl.org"
FLANK = 3
random.seed(115)


def get(path):
    req = urllib.request.Request(REST + path, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as r:
        return json.load(r)


def coding_positions(gene):
    g = get(f"/lookup/symbol/homo_sapiens/{gene}?expand=1")
    tx = next(t for t in g["Transcript"] if t.get("is_canonical"))
    tl = tx["Translation"]
    cds_lo, cds_hi = tl["start"], tl["end"]
    pos = set()
    for ex in tx["Exon"]:
        lo, hi = max(ex["start"], cds_lo), min(ex["end"], cds_hi)
        if lo > hi:
            continue
        pos.update(range(lo - FLANK, hi + FLANK + 1))
    return g["seq_region_name"], g["strand"], sorted(pos), tx["id"]


def rand_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def main():
    out, genes = sys.argv[1], sys.argv[2:]
    records = set()
    for gene in genes:
        chrom, strand, pos, tx = coding_positions(gene)
        lo, hi = pos[0] - 10, pos[-1] + 10
        seq = get(f"/sequence/region/human/{chrom}:{lo}..{hi}:1?coord_system_version=GRCh38")["seq"].upper()
        at = lambda p, n=1: seq[p - lo:p - lo + n]
        stops = ["TGA", "TAA", "TAG"] if strand == 1 else ["TCA", "TTA", "CTA"]
        for p in pos:
            r = at(p)
            if r not in "ACGT":
                continue
            # SNVs
            for b in "ACGT":
                if b != r:
                    records.add((chrom, p, r, b))
            # 2 and 3 bp substitutions
            for n in (2, 3):
                ref = at(p, n)
                alt = "".join(random.choice([b for b in "ACGT" if b != x]) for x in ref)
                records.add((chrom, p, ref, alt))
            # deletions (VCF anchor base + n deleted)
            for n in (1, 2, 3, 4, 6):
                records.add((chrom, p, at(p, n + 1), r))
            # insertions
            for b in "ACGT":
                records.add((chrom, p, r, r + b))
            for n in (3, 3, 6):
                records.add((chrom, p, r, r + rand_seq(n)))
            for s in stops:
                records.add((chrom, p, r, r + s))
            # duplications of the n bases ending at p
            for n in (3, 6):
                records.add((chrom, p, r, r + seq[p - lo - n + 1:p - lo + 1]))
        print(f"{gene} {tx} chr{chrom} strand {strand}: {len(pos)} positions", file=sys.stderr)

    order = {c: i for i, c in enumerate(sorted({c for c, *_ in records}, key=lambda c: (len(c), c)))}
    with open(out, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=make_candidates_v2.py (synthetic variants, no patient data)\n")
        for c in order:
            fh.write(f"##contig=<ID={c}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIMULATED\n")
        for c, p, ref, alt in sorted(records, key=lambda x: (order[x[0]], x[1], x[2], x[3])):
            fh.write(f"{c}\t{p}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:AD\t0/1:20,10\n")
    print(f"{len(records)} records", file=sys.stderr)


if __name__ == "__main__":
    main()
