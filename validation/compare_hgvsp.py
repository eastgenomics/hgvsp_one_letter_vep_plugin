#!/usr/bin/env python3
"""Compare HGVSp between a VEP run without and with a one-letter plugin.

Joins the two annotated VCFs on variant + transcript. Checks that the plugin
run's HGVSp is unchanged from the base run, and that its HGVSp_OneLetter equals
an independent code-by-code conversion of HGVSp. Also catalogues every distinct
HGVSp pattern seen.

Usage: python3 compare_hgvsp.py <base.vcf.gz> <plugin.vcf.gz> <out_prefix>
Writes <out_prefix>.rows.tsv and <out_prefix>.patterns.tsv
"""
import collections
import re
import subprocess
import sys

AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Xaa": "X", "Ter": "*",
}
AA3_RE = re.compile("|".join(AA3))
FMT = "%CHROM:%POS:%REF:%ALT\t%Feature\t%SYMBOL\t%Consequence\t%HGVSc\t%HGVSp\t{one}\n"


def rows(vcf, one_field):
    fmt = FMT.format(one="%HGVSp_OneLetter" if one_field else ".")
    out = subprocess.run(["bcftools", "+split-vep", "-d", "-f", fmt, vcf],
                         check=True, capture_output=True, text=True).stdout
    d = {}
    for line in out.splitlines():
        var, tx, sym, csq, hgvsc, hgvsp, one = line.split("\t")
        d[(var, tx)] = (sym, csq, hgvsc, hgvsp, one)
    return d


def desc(hgvs):
    return hgvs.split(":", 1)[1] if ":" in hgvs else hgvs


def pattern(p):
    p = AA3_RE.sub("Aaa", p)
    p = re.sub(r"(Aaa)+", "Aaa", p)
    return re.sub(r"\d+", "#", p)


def main():
    base_vcf, one_vcf, prefix = sys.argv[1:4]
    base, one = rows(base_vcf, False), rows(one_vcf, True)
    status = collections.Counter()
    pats = collections.defaultdict(lambda: collections.Counter())
    examples = {}

    with open(prefix + ".rows.tsv", "w") as fh:
        fh.write("variant\ttranscript\tsymbol\tconsequence\tHGVSc\tHGVSp_3letter\tHGVSp_plugin\texpected\tstatus\tpattern\n")
        for key in sorted(set(base) | set(one)):
            if key not in base or key not in one:
                status["missing_row"] += 1
                continue
            sym, csq, hgvsc, b = base[key][:4]
            p = one[key][4]
            if one[key][3] != b:
                status["HGVSp_changed_by_plugin"] += 1
                continue
            if b == ".":
                st = "ok_both_empty" if p == "." else "plugin_added_hgvsp"
                status[st] += 1
                continue
            if b.split(":", 1)[0] != p.split(":", 1)[0]:
                st = "protein_id_changed"
            else:
                exp = AA3_RE.sub(lambda m: AA3[m.group(0)], desc(b))
                got = desc(p)
                if got != exp:
                    st = "mismatch"
                elif AA3_RE.search(got):
                    st = "residual_3letter"
                else:
                    st = "ok"
            pat = pattern(desc(b))
            status[st] += 1
            pats[pat][st] += 1
            examples.setdefault(pat, (sym, desc(hgvsc), desc(b), desc(p)))
            fh.write("\t".join([key[0], key[1], sym, csq, hgvsc, b, p,
                                AA3_RE.sub(lambda m: AA3[m.group(0)], desc(b)), st, pat]) + "\n")

    with open(prefix + ".patterns.tsv", "w") as fh:
        fh.write("pattern\trows\tok\tnot_ok\texample_symbol\texample_HGVSc\texample_3letter\texample_plugin\n")
        for pat, c in sorted(pats.items(), key=lambda x: -sum(x[1].values())):
            n = sum(c.values())
            fh.write("\t".join([pat, str(n), str(c["ok"]), str(n - c["ok"]), *examples[pat]]) + "\n")

    for k, v in sorted(status.items()):
        print(f"{k}\t{v}")
    print(f"distinct_patterns\t{len(pats)}")


if __name__ == "__main__":
    main()
