# vep_plugin_hgvsp_one_letter

A VEP plugin, `HGVSpOneLetter`, that adds two fields:

- `HGVSp_OneLetter`: VEP's HGVSp with single-letter amino acid codes
- `HGVSp_SyntaxFlag`: `hgvs_unparseable` when VEP's three-letter HGVSp cannot be parsed by the biocommons `hgvs` package, empty otherwise (see [HGVSp_SyntaxFlag](#hgvsp_syntaxflag))

VEP's own three-letter `HGVSp` field is left unchanged, so both forms are available.

```
HGVSp             NP_000537.3:p.Arg175His   NP_000305.3:p.Lys6Ter   NP_000305.3:p.Met1_?2
HGVSp_OneLetter   NP_000537.3:p.R175H       NP_000305.3:p.K6*       NP_000305.3:p.M1_?2
HGVSp_SyntaxFlag                                                    hgvs_unparseable
```

## Why not Ensembl's SingleLetterAA?

VEP has no native option for one-letter HGVSp (VEP 115 `OutputFactory.pm` always calls `hgvs_protein()` in three-letter mode). Ensembl's `SingleLetterAA` plugin asks ensembl-variation to rebuild HGVSp in one-letter mode. That code path has bugs in release 115 and is unchanged in 116:

| Case | HGVSpOneLetter | SingleLetterAA |
|---|---|---|
| Frameshift whose first changed residue is a stop | `p.K6*` | `p.K6*fs*1` |
| Insertion/delins creating a stop | `p.Y27delins*` | `p.Y27delins*D` |
| Single-residue duplication | `p.D393dup` | `p.Asp393dup` |
| Stop-retained delins | `p.*394delins*` | `p.*394delins**` |

These come from `TranscriptVariationAllele.pm` (release/115) lines 1933, 1897 and 2400, which only handle three-letter codes. In one Atlas sample, 41 of 2,323 HGVSp values were wrong.

`HGVSpOneLetter` never uses that code path. VEP runs plugins after it has built each output line (`OutputFactory.pm` `run_plugins`), so the plugin reads the finished three-letter `HGVSp` and swaps each amino acid code for its one-letter equivalent (`Ter` → `*`, `Sec` → `U`, `Pyl` → `O`, `Xaa` → `X`). Everything else is passed through unchanged, which means:

- `--hgvsp_use_prediction` (`p.(R175H)`), `--remove_hgvsp_version` and `--no_escape` are honoured
- `HGVSp_OneLetter` always describes exactly the same change as `HGVSp`

## Usage

Requires `--hgvs`.

```bash
cp HGVSpOneLetter.pm ~/.vep/Plugins/
vep -i input.vcf --hgvs --plugin HGVSpOneLetter --fields "...,HGVSp,HGVSp_OneLetter,HGVSp_SyntaxFlag"
```

### With eggd_vep (1.3.0)

Add a plugin entry and `dir_plugins` to the VEP config:

```json
"plugins": [
  {
    "name": "HGVSpOneLetter",
    "pm_file": "file-xxxx",
    "required_fields": "HGVSp_OneLetter,HGVSp_SyntaxFlag",
    "resource_files": []
  }
],
"additional_flags": ["dir_plugins=/data/Plugins"]
```

- `"resource_files": []` is required: eggd_vep iterates over it with jq and fails if it is missing.
- `required_fields` adds `HGVSp_OneLetter` and `HGVSp_SyntaxFlag` to `--fields` (they are appended after the other fields).
- **`dir_plugins` is required.** eggd_vep downloads each `pm_file` to `/home/dnanexus/Plugins` (`/data/Plugins` in the container), but the Ensembl image defaults to `--dir_plugins /plugins`, its bundled copy of VEP_plugins. Without the flag, a plugin that is not bundled in the image is skipped with only a warning (`Failed to compile plugin HGVSpOneLetter`) and the job still completes. VEP puts `dir_plugins` first on `@INC` (`Runner.pm`), so other plugins in the config are then also loaded from their `pm_file`.

Check that the output CSQ header contains `HGVSp_OneLetter`; a missing field means the plugin did not load.

## Limitations

The plugin converts VEP's HGVSp faithfully, so anything VEP gets wrong in three-letter form is also in `HGVSp_OneLetter`. Compared with VariantValidator, VEP 115:

- describes insertions or delins that create an immediate stop as `p.Ala3_Ile4insTer` / `p.Tyr315delinsTer` rather than `p.Ile4Ter` / `p.Tyr315Ter`
- writes deletions across the start codon as `p.Met1_?2` rather than `p.Met1?`
- writes stop-codon deletions as `p.Ter404delextTer7` rather than `p.Ter404IleextTer7`
- writes synonymous multi-codon substitutions as `p.LeuSer112%3D` (not valid HGVS)
- omits `ext` for delins within the stop codon (`p.Ter404delinsLeuArg`)

Of these, `hgvs` rejects (and `HGVSp_SyntaxFlag` flags) only `p.Met1_?2`, `p.Ter404delextTer7` and `p.LeuSer112%3D`; the others parse, they are just not the preferred description.

## HGVSp_SyntaxFlag

`HGVSp_SyntaxFlag` is `hgvs_unparseable` when the biocommons [`hgvs`](https://github.com/biocommons/hgvs) package (1.5.7) cannot parse VEP's three-letter HGVSp, and empty otherwise. The plugin does not call `hgvs`: it uses a Perl port of the protein rules of its grammar, so it flags exactly what `hgvs` rejects and nothing more. VEP's `%3D` escaping of `=` is unescaped first. Two deliberate differences: the string must be a protein description (`:p.`; `hgvs` would also parse `c.` and other types), and accessions and gene symbols are matched as ASCII only.

It is a format check, not a validity check. An empty flag does not mean the HGVSp is correct: like `hgvs`, the check does not look at positions, ranges, where a stop may appear, or the change against the protein sequence, so `p.Arg0His`, `p.Arg175_Gly173del` and `p.Ter394CysfsTer77` are not flagged. The flag is also empty when there is no HGVSp.

VEP forms that are flagged:

- `p.Met1_?2`, `p.MetThr1_?2` (start-codon deletion)
- `p.Ter404delextTer7` (stop-codon deletion)
- `p.LeuSer112%3D` (synonymous change across codons)

Residue codes `hgvs` does not support are also flagged: `Pyl`, `Asx` and `Glx` anywhere, and `Xaa` as the residue at a position (`p.Xaa10Ala`; `hgvs` reads the `X` as a stop). `Xaa` as a new residue (`p.Ala10Xaa`) is not flagged.

`HGVSp_OneLetter` is still written for flagged rows, converted code for code.

## Tests

```bash
prove -v t/
```

The tests use a stub of VEP's plugin base class (`t/lib`), so no VEP installation is needed.
