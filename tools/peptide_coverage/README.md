# Peptide reference-coverage audit

A working, dependency-free Python 3 command-line tool for checking whether unmodified peptide strings occur in a protein search FASTA. It prevents a processed search non-hit from being treated as evidence of normal-tissue safety when the candidate sequence is absent from the search reference.

This is research tooling built around published targets. It is not a new biomarker, detection assay, therapeutic recommendation, clinical safety test or measured benchmark win.

## Run

From the repository root:

```sh
python3 tools/peptide_coverage/audit.py --fasta YOUR_REFERENCE.fasta \
  --peptide KLHGILVEA --peptide VVHLIKNAY \
  --expected-sha256 EXPECTED_SHA256 --expected-bytes EXPECTED_BYTES
python3 -m unittest discover -s tools/peptide_coverage/tests -v
```

Omit the expected hash/byte arguments only when file identity is unknown. Without an expected hash, negative findings stay explicitly unverified. Read the expected hash from a trusted source inventory, not a possibly incomplete current download. The hash is over file bytes, including gzip compression if used.

The tool prints JSON to stdout. Redirect it to your chosen output path. Input errors or expected hash/byte mismatches stop with exit code 2 and no JSON conclusions. A successful audit returns 0 even if a candidate has no exact match. Repeated `--peptide` arguments are deduplicated. `--max-hits` limits displayed matches per peptide, not the exact occurrence count. Default: five. Locations are 1-based inclusive.

No network requests or assay reads occur. The full file is parsed before any report is printed. Wrapped FASTA sequences and gzip input are supported. Matches never cross record boundaries, stop symbols or unknown residues. Duplicate headers remain distinguishable through record indices. Exact matching keeps I and L distinct and refuses modified, lowercase or ambiguous query peptides. Only protein FASTA input is appropriate; the parser cannot establish biological molecule type from letters alone.

## Status meanings

- `exact_reference_match`: candidate string occurs in at least one reference record. This establishes reference inclusion only.
- `not_in_verified_reference`: zero exact matches, expected file hash verified, and no noncanonical reference symbols found. This is file-level exact-string absence only.
- `no_exact_match_in_unverified_or_ambiguous_reference`: zero exact matches but no expected hash was supplied or at least one record contains symbols outside the 20 canonical amino acids plus stop. X, U and other noncanonical symbols are conservatively flagged; they do not necessarily invalidate the entire search, but they prevent this tool's strongest negative label.

None establishes whether the experimental pipeline actually used that reference, whether its search settings could recover the peptide, whether the peptide is displayed, whether a source is healthy, or whether it is a safe target. Source-to-model/HLA/condition mapping and participant independence remain separate gates.

## Measured reference-only examples

The bundled examples were generated October 1, 2026 from previously downloaded reference FASTAs. No PSM, peptide assay, intensity, clinical label or held-out outcome files were opened.

| Deposit | Bytes | SHA-256 | Parsed records | Noncanonical-symbol records |
|---|---:|---|---:|---:|
| PXD056070 | 34,539,419 | 36879acd229bf1515f974892cda7210a874d7fb48a75252c5ec15c9ec532cbbd | 76,094 | 8,233 |
| PXD055862 | 52,334,907 | 6f83f3353f33aff67a7796cd01a55084b92a89440765537e4e3295cf8cb37464 | 102,372 | 8,233 |

Both files have two exact occurrences of ovarian ITGB2 VVHLIKNAY and one of CTCFL KLHGILVEA. Both have zero exact occurrences of published NU11 KLFLWPYKV and LUAD ILMILQPQL. Each file contains 8,209 X characters and 36 U characters; the zero-match statuses therefore retain the ambiguity qualifier. Existing processed non-hits cannot support safety claims for the latter candidates merely by comparison with this reference.

The PXD056070 parsed record count corrects an earlier manual inventory count of 76,113. Hash and byte identity are unchanged. Example JSON contains reference headers and exact match coordinates, with all interpretation limits alongside the result. This is a deterministic software/reference check, not an independent scientific validation.

## Sources and provenance limits

- PXD056070: https://www.ebi.ac.uk/pride/archive/projects/PXD056070
- PXD055862: https://www.ebi.ac.uk/pride/archive/projects/PXD055862
- Beta-cell study: https://pmc.ncbi.nlm.nih.gov/articles/PMC12645175/
- Ovarian CTCFL nomination: https://pmc.ncbi.nlm.nih.gov/articles/PMC10070997/
- Ovarian ITGB2 nomination: https://www.nature.com/articles/s41541-025-01234-6
- NU11 nomination: https://pmc.ncbi.nlm.nih.gov/articles/PMC12163983/
- LUAD nomination: https://www.nature.com/articles/s43018-025-00979-2

Accession association and expected hashes come from the preserved source audit. They are not signatures issued by the publisher. The original reference binaries are not redistributed here. The archive landing pages are source pointers, not direct download URLs. A fresh retrieval of PXD055862's landing page failed during packaging; this package does not claim to have re-downloaded the reference today. Raw-file-to-condition/model/HLA mapping remains unresolved for these sources.
