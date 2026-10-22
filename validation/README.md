# Validation

All validation used VEP 115 via eggd_vep 1.3.0 (`--refseq --hgvs --shift_3prime 1`, GRCh38) in DNAnexus project `003_261021_vep-single-letter-aa` (`project-JFG9V704JJvggp43gJfXXQ7g`), plugin config `additional_flags: ["dir_plugins=/data/Plugins"]`. Test data in this repository is synthetic only.

## Results

| Test | Rows with HGVSp | `HGVSp_OneLetter` correct |
|---|---|---|
| Synthetic TP53/PTEN/GPX1, default flags | 531,982 | 531,982 |
| Synthetic, `--hgvsp_use_prediction --no_escape` | 531,982 | 531,982 |
| Reproducer VCF (20 cases) | 164 | 164 |
| Atlas sample (full Atlas config) | 2,323 | 2,323 |

"Correct" means equal to an independent code-by-code conversion of the same row's `HGVSp` (`compare_hgvsp.py`). In the plugin runs `HGVSp` itself was unchanged from runs without the plugin. On the Atlas sample, every other CSQ field and every non-INFO column was identical to production.

The synthetic runs produced 19 distinct HGVSp patterns (`data/simulated_hgvsp_patterns.tsv`), including selenocysteine (`p.Sec49Arg` → `p.U49R`), `?` endings, `ext`, and `%3D`.

### biocommons hgvs

`hgvs_parse_check.py` parsed all 764,890 distinct HGVSp strings from the synthetic runs with `hgvs` 1.5.7 and formatted them with `p_3_letter=False`:

- 764,122 parsed, and every one matched the plain code swap
- 768 were rejected as invalid HGVS; all are VEP output forms listed under "VEP limitations" below

### VariantValidator

`check_variantvalidator.py` queried the public VariantValidator API (GRCh38, `mane_select`) for 137 synthetic variants covering 18 patterns (`data/vv_sample.tsv`, results in `data/vv_results.tsv`):

| Result | Variants |
|---|---|
| VEP three-letter = VV `tlr`, and HGVSp_OneLetter = VV `slr` | 90 |
| VEP three-letter = VV `tlr`, one-letter differs | 0 |
| VEP and VV describe the change differently | 44 |
| No VV result | 3 |

The 44 differences are in VEP's three-letter `HGVSp`, so they are carried into `HGVSp_OneLetter` (see below).

## VEP limitations (inherited)

| Situation | VEP 115 | VariantValidator / HGVS |
|---|---|---|
| Insertion creating an immediate stop | `p.Ala3_Ile4insTer` | `p.Ile4Ter` |
| Delins ending in a stop | `p.Tyr315delinsTer` | `p.Tyr315Ter` |
| Deletion across the start codon | `p.Met1_?2`, `p.MetThr1_?2` | `p.Met1?` |
| Deletion of the stop codon | `p.Ter404delextTer7` | `p.Ter404IleextTer7` |
| Synonymous substitution across two codons | `p.LeuSer112%3D` | not valid HGVS |
| Delins within the stop codon | `p.Ter404delinsLeuArg` | missing `ext` |
| Selenoproteins (GPX1) | `extTer?`, `fsTer?` | VV reads TGA as selenocysteine; review manually |

The malformed forms do occur in practice: in 178 Atlas validation VCFs (pre-filtering), 272 of 579,740 HGVSp values were rejected by `hgvs`: synonymous multi-codon substitutions (174), start-codon deletions (71) and stop-codon deletions (27).

## SingleLetterAA reproducer

`data/SingleLetterAA_reproducer_GRCh38.vcf`: 20 synthetic variants on TP53 (NM_000546.6) and PTEN (NM_000314.8), labelled in `INFO/CASE` (`bug_*`, `diff_*`, `control_*`). `data/SingleLetterAA_reproducer_GRCh38_expected.tsv` gives VEP's three-letter HGVSp and Ensembl SingleLetterAA's output for each.

## Scripts

| Script | Purpose |
|---|---|
| `make_candidates_v2.py` | Generate synthetic variants (SNVs, MNVs, indels, dups, stop insertions) across the CDS of given genes, from the Ensembl REST API |
| `compare_hgvsp.py` | Compare a VEP run without and with the plugin, row by row; writes per-row results and a pattern catalogue |
| `hgvs_parse_check.py` | Cross-check with biocommons `hgvs` |
| `check_variantvalidator.py` | Compare a sample against VariantValidator. **Synthetic variants only** - it sends variants to a public API |

`make_candidates_v2.py` writes contig `MT` variants if asked for mitochondrial genes; leave them out for eggd_vep, whose bcftools reference has no `MT`.
