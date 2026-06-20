# Security Policy

bamsteg ships both a steganography embedder and a forensic detector for genomic (BAM) files. Two categories of report are relevant here, and they're handled differently.

## Detector bypass / new evasion technique

If you've found a way to hide a payload that `bamsteg detect` does not flag, **please do not open a public issue or PR demonstrating it.** Report it privately first (see contact below) so a detection update can ship before the bypass is public. Once addressed, it'll be credited and documented openly in the changelog.

## Implementation vulnerabilities

Bugs in the crypto (`crypto.py`), error correction (`ecc.py`), or anything that weakens confidentiality of an embedded payload beyond its intended threat model: report the same way.

## How to report

Open a private report via GitHub's **Security Advisories** tab on this repo (Security → Report a vulnerability), rather than a public issue. If you'd prefer email, add a contact address here before publishing.

## Scope

This project assumes the embedder is used by someone with legitimate access to the BAM files they're modifying. Reports about misuse of the tool by unauthorized third parties are a policy/responsible-use matter, not a vulnerability. See the README's Responsible Use section.
