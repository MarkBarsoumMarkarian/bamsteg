# bamsteg

Format-preserving steganography and forensic steganalysis for BAM (genomic sequence alignment) files.

bamsteg is a paired offense/defense research tool. It includes two working methods for hiding encrypted, arbitrary-content payloads inside a standard BAM file without breaking its structure, alongside a forensic detector purpose-built to catch them.

## Why this exists

Genomic data pipelines move enormous binary files between labs, sequencing cores, cloud storage, and collaborators with comparatively little content inspection. Quality scores and auxiliary tags are treated as inert sequencing byproduct, and nobody scans them for hidden payloads. That assumption is the vulnerability this project documents.

bamsteg was built to answer one question: **if someone hid data inside genomic files moving through a real bioinformatics pipeline, would anyone notice?** The embedder proves the channel exists. The detector exists so the answer can be "yes."

This is a research and forensic-tooling project. It is not intended, and should not be used, for unauthorized data exfiltration from systems you don't own or have explicit permission to test.

## How it works

- **AES-256-GCM** encryption (PBKDF2-HMAC-SHA256 key derivation, per-message salt and nonce) wraps every payload before it touches the BAM file. Nothing is embedded in plaintext.
- **Error-correcting code** wraps the ciphertext so the hidden payload can survive realistic processing (e.g. re-sorting) without becoming unrecoverable.
- **Passphrase-seeded read selection**: a PRNG seeded from the passphrase chooses which reads carry payload bits, so the carrier set isn't sequential or guessable without the key.
- **Two embedding modes:**
  - `aux`: payload bytes encoded across a custom auxiliary tag on selected reads.
  - `lsb`: payload bits hidden in the least-significant-bit of Phred quality scores.
- **`bamsteg detect`**: forensic steganalysis combining BAM header signature scanning, auxiliary tag scanning, and per-read Shannon entropy analysis on quality-score LSBs to flag anomalous files.

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
bamsteg detect --input target.bam --sample-size 10000 --entropy-threshold 0.95
```

`detect` exits non-zero and prints a threat level (LOW / HIGH / CRITICAL) when it flags anomalies, so it's pipeline-friendly. Drop it into a CI step or an intake script for sequencing data.

## Detection limitations

No single heuristic generalizes to every possible encoding scheme. `bamsteg detect` is tuned against the methods this project ships, not as a universal genomic-steganalysis solution. Treat a clean scan as "no known signature found," not as a guarantee of integrity. This is an active area the project intends to keep improving.

## Responsible use

This project is published for security research, bioinformatics forensics, and education. Do not use it to hide or move data through systems, pipelines, or organizations without explicit authorization. If you work in genomics security and this exposes a gap in your own infrastructure, that's the point. Go fix it.

## License

MIT, see [LICENSE](LICENSE).

---

*To anyone using technology to hide what should be seen, or to do harm in the dark: this is built to drag it back into the light.*

> "Whose hatred is covered by deceit, his wickedness shall be shewed before the whole congregation." Proverbs 26:26, KJV
