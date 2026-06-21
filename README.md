# bamsteg

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20778947.svg)](https://doi.org/10.5281/zenodo.20778947)

Format-preserving steganography and forensic steganalysis for BAM (genomic sequence alignment) files.

bamsteg is a paired offense/defense research tool. It includes two working methods for hiding encrypted, arbitrary-content payloads inside a standard BAM file without breaking its structure, alongside a forensic detector purpose-built to catch them.

## Why this exists

Genomic data pipelines move enormous binary files between labs, sequencing cores, cloud storage, and collaborators with comparatively little content inspection. Quality scores and auxiliary tags are treated as inert sequencing byproduct, and nobody scans them for hidden payloads. That assumption is the vulnerability this project documents.

bamsteg was built to answer one question: **if someone hid data inside genomic files moving through a real bioinformatics pipeline, would anyone notice?** The embedder proves the channel exists. The detector exists so the answer can be "yes."

This is a research and forensic-tooling project. It is not intended, and should not be used, for unauthorized data exfiltration from systems you do not own or have explicit permission to test.

## How it works

- **AES-256-GCM** encryption (PBKDF2-HMAC-SHA256 key derivation, per-message salt and nonce) wraps every payload before it touches the BAM file. Nothing is embedded in plaintext.
- **Reed-Solomon error correction** wraps the ciphertext so the hidden payload can survive realistic BAM processing (e.g. re-sorting, deduplication) without becoming unrecoverable.
- **Passphrase-seeded read selection**: a PRNG seeded from the passphrase chooses which reads carry payload bits, so the carrier set is not sequential or guessable without the key.
- **Two embedding modes:**
  - `aux`: payload bytes encoded across a custom auxiliary tag (`ZS`) on selected reads. Zero impact on SEQ, QUAL, CIGAR, or POS -- variant callers see an identical file.
  - `lsb`: payload bits hidden in the least-significant bit of Phred quality scores. A ±1 quality delta (e.g. Q30 vs Q31) is invisible to variant calling pipelines.

## Detection (`bamsteg detect`)

The detector runs four independent steganalysis methods and combines them into a single threat level (CLEAN / SUSPICIOUS / HIGH / CRITICAL):

1. **Header signature scan** -- checks `@CO` and `@PG` fields for known steganography tool signatures.
2. **Generic aux tag scan** -- flags any `X*/Y*/Z*` tag (SAM spec reserved namespace) carrying a hex-pattern covert channel value, not just bamsteg-specific tag names.
3. **Adaptive LSB chi-squared** -- splits each read into a head zone (potential payload) and tail zone (used as instrument baseline), then chi-squares their LSB parity distributions. Instrument-independent: no hardcoded quality score assumptions.
4. **RS-SPA (Sample Pairs Analysis)** -- measures the P1/P2 adjacent-score pair ratio. In clean Illumina data P1 < P2 naturally; LSB flipping inflates P1 toward and above P2. Detects low-rate embedding (<3% of reads) where chi-squared loses statistical power.

`detect` exits with code 0 (clean) or 1 (anomaly found), making it pipeline-friendly as a BAM intake check.

## Installation

```bash
git clone https://github.com/MarkBarsoumMarkarian/bamsteg.git
cd bamsteg
pip install -e .
```

Requires `pysam`, `pycryptodome`, and `reedsolo`.

## Usage

```bash
# Embed a payload (auxiliary-tag mode)
bamsteg embed --mode aux --input in.bam --output out.bam --payload secret.txt --passphrase "key"

# Embed a payload (quality-score LSB mode)
bamsteg embed --mode lsb --input in.bam --output out.bam --payload secret.txt --passphrase "key"

# Extract a payload
bamsteg extract --mode aux --input out.bam --passphrase "key" --out-payload recovered.txt

# Run forensic steganalysis on an unknown BAM file
bamsteg detect --input target.bam

# Scan with a custom sample size (default: 50,000 reads)
bamsteg detect --input target.bam --sample 10000
```

## Detection limitations

No single heuristic generalizes to every possible encoding scheme. `bamsteg detect` is tuned against the methods this project ships and may miss novel embedding strategies or non-standard implementations. Treat a clean scan as "no known signature found," not as a guarantee of integrity.

## Responsible use

This project is published for security research, bioinformatics forensics, and education. Do not use it to hide or move data through systems, pipelines, or organizations without explicit authorization. If you work in genomics security and this exposes a gap in your own infrastructure, that is the point.

## Citation

If you use bamsteg in research, please cite the Zenodo record:

```
Markarian, M. B. (2026). bamsteg: Format-preserving steganography and
steganalysis for NGS BAM files. Zenodo. https://doi.org/10.5281/zenodo.20778947
```

## License

MIT, see [LICENSE](LICENSE).

---

*To anyone using technology to hide what should be seen, or to do harm in the dark: this is built to drag it back into the light.*

> "Whose hatred is covered by deceit, his wickedness shall be shewed before the whole congregation." Proverbs 26:26, KJV
