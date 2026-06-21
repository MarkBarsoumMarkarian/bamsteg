"""
bamsteg CLI -- Format-preserving steganography for NGS BAM files.

Commands:
  bamsteg embed   --mode [aux|lsb] --input in.bam --output out.bam --payload file --passphrase key
  bamsteg extract --mode [aux|lsb] --input out.bam --passphrase key [--out-payload file]
  bamsteg detect  --input target.bam [--sample N]
"""

import argparse
import sys
import os
import getpass

from .embed_aux import embed_aux, extract_aux
from .embed_lsb import embed_lsb, extract_lsb
from .bamsteg_detect import run_forensics, SAMPLE_SIZE


def cmd_embed(args):
    passphrase = args.passphrase or getpass.getpass("Passphrase: ")
    if not os.path.exists(args.input):
        print(f"[ERROR] Input BAM not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    with open(args.payload, "rb") as f:
        payload = f.read()
    print(f"[*] Payload: {len(payload)} bytes | Mode: {args.mode}")
    print(f"[*] {args.input} -> {args.output}")
    if args.mode == "aux":
        stats = embed_aux(args.input, args.output, payload, passphrase)
    else:
        stats = embed_lsb(args.input, args.output, payload, passphrase)
    print("[+] Done.")
    for k, v in stats.items():
        print(f"    {k}: {v}")


def cmd_extract(args):
    passphrase = args.passphrase or getpass.getpass("Passphrase: ")
    if not os.path.exists(args.input):
        print(f"[ERROR] Input BAM not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    if args.mode == "aux":
        payload = extract_aux(args.input, passphrase)
    else:
        payload = extract_lsb(args.input, passphrase)
    if args.out_payload:
        with open(args.out_payload, "wb") as f:
            f.write(payload)
        print(f"[+] Payload written to: {args.out_payload} ({len(payload)} bytes)")
    else:
        try:
            print(payload.decode("utf-8"))
        except UnicodeDecodeError:
            print(f"<binary payload, {len(payload)} bytes -- use --out-payload to save>")


def cmd_detect(args):
    if not os.path.exists(args.input):
        print(f"[ERROR] Target BAM not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    level, alerts, _ = run_forensics(args.input, sample_size=args.sample)
    # exit code 0 = clean, 1 = anything suspicious or above
    # makes bamsteg detect usable in CI/pipeline scripts
    sys.exit(0 if level == "CLEAN" else 1)


def main():
    parser = argparse.ArgumentParser(
        prog="bamsteg",
        description="Format-preserving steganography + steganalysis for NGS BAM files.",
    )
    sub = parser.add_subparsers(dest="command")

    # --- embed ---
    pe = sub.add_parser("embed", help="Embed encrypted payload into BAM")
    pe.add_argument("--mode", choices=["aux", "lsb"], default="aux",
                    help="aux: custom ZS tags (zero bio impact) | lsb: quality score LSBs (higher capacity)")
    pe.add_argument("--input",      required=True, help="Input BAM")
    pe.add_argument("--output",     required=True, help="Output BAM with payload")
    pe.add_argument("--payload",    required=True, help="File to embed")
    pe.add_argument("--passphrase", default=None,  help="Encryption passphrase (prompted if omitted)")

    # --- extract ---
    px = sub.add_parser("extract", help="Extract payload from embedded BAM")
    px.add_argument("--mode", choices=["aux", "lsb"], default="aux")
    px.add_argument("--input",       required=True, help="BAM containing payload")
    px.add_argument("--passphrase",  default=None)
    px.add_argument("--out-payload", default=None,  help="Write payload to file (else print to stdout)")

    # --- detect ---
    pd = sub.add_parser("detect", help="Forensic steganalysis: scan BAM for hidden payloads")
    pd.add_argument("--input",  required=True,          help="BAM file to scan")
    pd.add_argument("--sample", type=int, default=SAMPLE_SIZE,
                    help=f"Reads to sample (default: {SAMPLE_SIZE})")

    args = parser.parse_args()

    if args.command == "embed":
        cmd_embed(args)
    elif args.command == "extract":
        cmd_extract(args)
    elif args.command == "detect":
        cmd_detect(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
