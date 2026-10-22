#!/usr/bin/env python3
"""Check plugin one-letter HGVSp against VariantValidator for a sample of rows.

SYNTHETIC VARIANTS ONLY - this sends variants to the public VariantValidator
API. Never run it on patient-derived data.

For each sampled row (variant + RefSeq transcript) it queries VariantValidator
(GRCh38, mane_select), then:
  1. compares VEP's three-letter HGVSp with VV's tlr
  2. where those agree, compares the plugin's one-letter HGVSp with VV's slr

Usage: python3 check_variantvalidator.py <sample.tsv> <cache_dir> <out.tsv>
sample.tsv columns: variant (chrom:pos:ref:alt), transcript, HGVSp_3letter, HGVSp_plugin
"""
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

API = "https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/{}/mane_select?content-type=application%2Fjson"


def query(variant, cache_dir):
    path = os.path.join(cache_dir, variant.replace(":", "_") + ".json")
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    url = API.format(urllib.parse.quote(variant.replace(":", "-")))
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                data = json.load(r)
            break
        except Exception as e:  # rate limit / timeout: back off and retry
            print(f"retry {variant}: {e}", file=sys.stderr)
            time.sleep(10 * (attempt + 1))
    else:
        return None
    with open(path, "w") as fh:
        json.dump(data, fh)
    time.sleep(1.5)
    return data


def norm(p):
    """Drop protein ID and prediction brackets; unescape '='."""
    if p in (None, "", "."):
        return ""
    p = p.split(":", 1)[-1].replace("%3D", "=")
    return re.sub(r"^p\.\((.*)\)$", r"p.\1", p)


def main():
    sample, cache_dir, out = sys.argv[1:4]
    os.makedirs(cache_dir, exist_ok=True)
    with open(sample) as fh:
        todo = list(csv.DictReader(fh, delimiter="\t"))

    with open(out, "w") as fh:
        fh.write("variant\ttranscript\tVEP_3letter\tVV_tlr\tplugin_1letter\tVV_slr\tVV_warnings\tthree_letter_agree\tone_letter_agree\n")
        for i, row in enumerate(todo, 1):
            data = query(row["variant"], cache_dir)
            vv = None
            for key, val in (data or {}).items():
                if isinstance(val, dict) and key.split(":")[0] == row["transcript"]:
                    vv = val
            if vv is None:
                fh.write("\t".join([row["variant"], row["transcript"], row["HGVSp_3letter"], "", row["HGVSp_plugin"], "", "no_VV_result_for_transcript", "", ""]) + "\n")
                continue
            pc = vv.get("hgvs_predicted_protein_consequence") or {}
            tlr, slr = norm(pc.get("tlr")), norm(pc.get("slr"))
            vep3, one = norm(row["HGVSp_3letter"]), norm(row["HGVSp_plugin"])
            agree3 = vep3 == tlr
            agree1 = (one == slr) if agree3 else ""
            warn = " | ".join(vv.get("validation_warnings") or []).replace("\t", " ")
            fh.write("\t".join([row["variant"], row["transcript"], vep3, tlr, one, slr, warn, str(agree3), str(agree1)]) + "\n")
            print(f"{i}/{len(todo)} {row['variant']} 3L:{agree3} 1L:{agree1}", file=sys.stderr)


if __name__ == "__main__":
    main()
