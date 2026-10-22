#!/usr/bin/env python3
"""Cross-check HGVSpOneLetter conversion with the biocommons hgvs package.

Parses every distinct three-letter HGVSp in compare_hgvsp.py row files with
hgvs, formats it with one-letter codes, and compares with a plain code swap.
Reports strings hgvs rejects as invalid HGVS.

Needs: hgvs (install --no-deps with bioutils, Parsley, attrs).
Usage: python3 hgvs_parse_check.py <rows.tsv> [<rows.tsv> ...]
"""
import csv, sys, collections, re
import hgvs.parser
from hgvs.config import global_config
AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Sec":"U","Pyl":"O","Xaa":"X","Ter":"*"}
RX = re.compile("|".join(AA3))
hp = hgvs.parser.Parser()
uniq = {}
for f in sys.argv[1:]:
    for r in csv.DictReader(open(f), delimiter="\t"):
        uniq.setdefault(r["HGVSp_3letter"], r["pattern"])
res = collections.Counter(); ex = collections.defaultdict(list)
for h, pat in uniq.items():
    simple = h.split(":",1)[0] + ":" + RX.sub(lambda m: AA3[m.group(0)], h.split(":",1)[1])
    try:
        v = hp.parse_hgvs_variant(h.replace("%3D", "="))
        out = v.format(conf={"p_3_letter": False})
    except Exception as e:
        k = "parse_error"; res[k] += 1
        if len(ex[k+pat]) < 2: ex[k+pat].append((h, type(e).__name__ + ": " + str(e)[:90]))
        continue
    k = "same_as_simple" if out == simple.replace("%3D", "=") else "differs"
    res[k] += 1
    if k == "differs" and len(ex[k+pat]) < 2: ex[k+pat].append((h, out + "   vs simple " + simple))
print(len(uniq), "distinct HGVSp strings"); [print(k, v) for k, v in res.items()]
for k, v in sorted(ex.items()):
    print("==", k); [print("  ", a, "->", b) for a, b in v]
