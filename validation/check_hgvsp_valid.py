#!/usr/bin/env python3
"""Check that HGVSp_Valid agrees with the biocommons hgvs parser.

Runs the plugin's valid_hgvsp() and hgvs parse_hgvs_variant() on the same
strings and reports every disagreement. Strings are distinct HGVSp_3letter
values from compare_hgvsp.py row files, plus (with --fuzz N) N random
protein-like strings built from residue codes, positions and edit keywords.

Needs: hgvs (install --no-deps with bioutils, Parsley, attrs) and perl.
Usage: python3 check_hgvsp_valid.py [--fuzz N] [<rows.tsv> ...]
"""
import argparse, csv, os, random, subprocess, collections
import hgvs.parser

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R3 = "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Sec Pyl Xaa Ter Asx Glx".split()
R1 = list("ACDEFGHIKLMNPQRSTVWYBZXUO*")
TOK = R3 * 3 + R1 + ["0", "1", "5", "175", "-", "-5", "_", "del", "ins", "dup", "delins", "fs", "ext", "=", "%3D", "?", "(", ")"]
IDS = ["NP_000537.3", "NP_000537", "ENSP00000269305.4", "NP_1.1(TP53)", "NP-1", "1NP", "NP_1.", ""]

def fuzz(n, seed=1):
    rnd = random.Random(seed)
    pos = lambda: rnd.choice(R3 + R1[:5]) + rnd.choice(["1", "5", "175", "0", ""])
    out = set()
    for _ in range(n // 2):  # interval + edit
        iv = pos() + ("_" + pos() if rnd.random() < 0.4 else "")
        e = rnd.choice(["", rnd.choice(R3 + R1), "?"]) + rnd.choice(["", "fs", "ext"]) + rnd.choice(["", rnd.choice(["Ter", "*", "X", "Met", "-"]) + rnd.choice(["", "?", "5", "-5"])])
        e = rnd.choice([e, e, "del", "dup", "=", "%3D", "delins" + "".join(rnd.choices(R3 + R1, k=rnd.randint(0, 3))),
                        "ins" + "".join(rnd.choices(R3 + R1, k=rnd.randint(0, 3))), "delext" + rnd.choice(["Ter7", ""])])
        d = iv + e
        out.add(rnd.choice(IDS[:4]) + ":p." + ("(" + d + ")" if rnd.random() < 0.2 else d))
    for _ in range(n - n // 2):  # random tokens
        out.add(rnd.choice(IDS) + rnd.choice([":p.", ":p.", ":p.", "p.", ":p"]) + "".join(rnd.choices(TOK, k=rnd.randint(1, 7))))
    return out

ap = argparse.ArgumentParser()
ap.add_argument("--fuzz", type=int, default=0)
ap.add_argument("rows", nargs="*")
a = ap.parse_args()
strings = fuzz(a.fuzz)
for f in a.rows:
    strings.update(r["HGVSp_3letter"] for r in csv.DictReader(open(f), delimiter="\t") if r["HGVSp_3letter"])
strings = sorted(strings)

perl = ("use lib '%s/t/lib', '%s'; use HGVSpOneLetter; "
        "while (<STDIN>) { chomp; print HGVSpOneLetter::valid_hgvsp($_) ? \"true\\n\" : \"false\\n\" }") % (REPO, REPO)
ours = subprocess.run(["perl", "-e", perl], input="\n".join(strings) + "\n", capture_output=True, text=True, check=True).stdout.split()

hp = hgvs.parser.Parser()
res = collections.Counter(); diff = []
for s, v in zip(strings, ours):
    try: hp.parse_hgvs_variant(s.replace("%3D", "=")); h = "true"
    except Exception: h = "false"
    res[(h, v)] += 1
    if h != v: diff.append((h, v, s))
print(len(strings), "distinct strings")
for (h, v), n in sorted(res.items()): print(f"hgvs={h} plugin={v}: {n}")
for h, v, s in diff[:50]: print(f"  DIFF hgvs={h} plugin={v} {s}")
