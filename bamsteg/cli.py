"""
bamsteg - Format-preserving steganography for BAM files.
Usage:
  bamsteg embed --mode aux --input in.bam --output out.bam --payload secret.txt --passphrase "key"
  bamsteg extract --mode aux --input out.bam --passphrase "key" --out-payload recovered.txt
  bamsteg detect --input target.bam --sample-size 10000 --entropy-threshold 0.95
"""

import argparse
import sys
import os
import getpass

from .embed_aux import embed_aux, extract_aux
from .embed_lsb import embed_lsb, extract_lsb
from .bamsteg_detect import scan_headers, analyze_reads

def cmd_embed(args):
    passphrase = args.passphrase or getpass.getpass("Passphrase: ")

    if not os.path.exists(args.input):
        print(f"[ERROR] Input BAM not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    with open(args.payload, "rb") as f:
        payload = f.read()

    print(f"[*] Payload size: {len(payload)} bytes")
    print(f"[*] Mode: {args.mode}")
    print(f"[*] Input: {args.input}")
    print(f"[*] Output: {args.output}")

    if args.mode == "aux":
        stats = embed_aux(args.input, args.output, payload, passphrase)
    elif args.mode == "lsb":
        stats = embed_lsb(args.input, args.output, payload, passphrase)
    else:
        print(f"[ERROR] Unknown mode: {args.mode}", file=sys.stderr)
        sys.exit(1)

    print("[+] Embedding complete.")
    for k, v in stats.items():
        print(f"    {k}: {v}")

def cmd_extract(args):
    passphrase = args.passphrase or getpass.getpass("Passphrase: ")

    if not os.path.exists(args.input):
        print(f"[ERROR] Input BAM not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"[*] Mode: {args.mode}")
    print(f"[*] Extracting from: {args.input}")

    if args.mode == "aux":
        payload = extract_aux(args.input, passphrase)
    elif args.mode == "lsb":
        if not args.payload_size:
            print("[ERROR] LSB mode requires --payload-size (bytes)", file=sys.stderr)
            sys.exit(1)
        payload = extract_lsb(args.input, passphrase, int(args.payload_size))
    else:
        print(f"[ERROR] Unknown mode: {args.mode}", file=sys.stderr)
        sys.exit(1)

    if args.out_payload:
        with open(args.out_payload, "wb") as f:
            f.write(payload)
        print(f"[+] Payload written to: {args.out_payload} ({len(payload)} bytes)")
    else:
        print(f"[+] Payload ({len(payload)} bytes):")
        try:
            print(payload.decode("utf-8"))
        except UnicodeDecodeError:
            print(f"  <binary data, {len(payload)} bytes>")

def cmd_detect(args):
    if not os.path.exists(args.input):
        print(f"[ERROR] Target BAM not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print("=" * 60)
    print(f"  BAMSTEG DETECT - FORENSIC ANALYSIS")
    print(f"  Target: {args.input}")
    print(f"  Sample Size: {args.sample_size} reads")
    print(f"  Entropy Threshold: {args.entropy_threshold}")
    print("=" * 60)

    header_alerts = scan_headers(args.input)

    try:
        read_metrics = analyze_reads(args.input, sample_size=args.sample_size, entropy_threshold=args.entropy_threshold)
    except TypeError:
        # Fallback if analyze_reads hasn't been updated with the entropy_threshold parameter yet
        read_metrics = analyze_reads(args.input, sample_size=args.sample_size)

    threat_level = "LOW (CLEAN)"
    color_code = "\033[92m" # Green
    alerts = []

    if len(header_alerts) > 0:
        threat_level = "CRITICAL"
        alerts.extend(header_alerts)

    if read_metrics.get("suspicious_tags", 0) > 0:
        threat_level = "HIGH"
        alerts.append(f"Found {read_metrics['suspicious_tags']} reads with suspicious 'ZS' carrier tags.")

    if read_metrics.get("high_entropy_reads", 0) > 0:
        if threat_level != "CRITICAL":
            threat_level = "HIGH"
        alerts.append(f"Found {read_metrics['high_entropy_reads']} reads exhibiting localized high entropy (>{args.entropy_threshold}). Indicates LSB ciphertext.")

    if threat_level == "HIGH":
        color_code = "\033[93m" # Yellow
    elif threat_level == "CRITICAL":
        color_code = "\033[91m" # Red

    reset_code = "\033[0m"

    print("\n--- RESULTS ---")
    print(f"Reads Analyzed       : {read_metrics.get('reads_scanned', 0)}")
    print(f"High-Entropy Reads   : {read_metrics.get('high_entropy_reads', 0)}")
    print(f"Suspicious AUX Tags  : {read_metrics.get('suspicious_tags', 0)}")

    print(f"\nTHREAT LEVEL: {color_code}{threat_level}{reset_code}")

    if alerts:
        print("\nForensic Flags Raised:")
        for alert in alerts:
            print(f" -> {alert}")
        print("=" * 60)
        sys.exit(1) # Return error code for pipeline automation
    else:
        print("\nNo steganographic anomalies detected. File appears biologically standard.")
        print("=" * 60)
        sys.exit(0) # Return success code

def main():
    parser = argparse.ArgumentParser(
        prog="bamsteg",
        description="Format-preserving steganography for NGS BAM files.",
    )
    sub = parser.add_subparsers(dest="command")

    # embed
    p_embed = sub.add_parser("embed", help="Embed payload into BAM file")
    p_embed.add_argument("--mode", choices=["aux", "lsb"], default="aux")
    p_embed.add_argument("--input", required=True, help="Input BAM file")
    p_embed.add_argument("--output", required=True, help="Output BAM file")
    p_embed.add_argument("--payload", required=True, help="File to embed")
    p_embed.add_argument("--passphrase", default=None, help="Encryption passphrase")

    # extract
    p_extract = sub.add_parser("extract", help="Extract payload from BAM file")
    p_extract.add_argument("--mode", choices=["aux", "lsb"], default="aux")
    p_extract.add_argument("--input", required=True, help="BAM file with payload")
    p_extract.add_argument("--passphrase", default=None)
    p_extract.add_argument("--out-payload", default=None, help="Output file for payload")
    p_extract.add_argument("--payload-size", default=None, help="Expected payload bytes (LSB mode only)")

    # detect (Forensics)
    p_detect = sub.add_parser("detect", help="Run forensic steganalysis on a BAM file")
    p_detect.add_argument("--input", required=True, help="Target BAM file to analyze")
    p_detect.add_argument("--sample-size", type=int, default=10000, help="Number of reads to scan (default: 10000)")
    p_detect.add_argument("--entropy-threshold", type=float, default=0.95, help="Threshold for LSB Shannon entropy (default: 0.95)")

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
