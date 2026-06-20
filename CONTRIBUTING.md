# Contributing to bamsteg

Contributions are welcome, especially on the detection side. That's the part of this project that matters most.

## Setup

```bash
git clone https://github.com/MarkBarsoumMarkarian/bamsteg.git
cd bamsteg
python3 -m venv env && source env/bin/activate
pip install -e .
```

Run the integration test before opening a PR:

```bash
python3 run_integration_test.py
```

## What's useful

- **New detection heuristics**: additional statistical or structural signals that catch embedding methods the current detector misses.
- **False-positive reduction**: real-world BAM files from varied sequencers/pipelines that trip `detect` incorrectly.
- **Bug reports** against the crypto, ECC, or read-selection logic.

## On new embedding/evasion techniques

If you find a way to embed data that bypasses `bamsteg detect`, please don't open a public PR or issue demonstrating it outright. Report it the same way you'd report a security vulnerability (see `SECURITY.md`) so the detector can be updated first. Once a fix ships, the technique and the fix can both be documented openly. That's the same coordinated-disclosure norm the rest of security research follows, and it keeps this project useful as a detector rather than turning into a how-to.

## Code style

Keep it readable over clever. This is a forensics tool, and code that's hard to audit defeats the purpose.
