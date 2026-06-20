#!/usr/bin/env python3
"""
bamsteg_detect.py -- Forensic Steganalysis Engine for BAM Files
Detects structural and statistical anomalies indicative of genomic steganography,
including localized high-entropy scanning for LSB ciphertext.
"""

import sys
import pysam
import math
from collections import Counter

def calculate_shannon_entropy(bit_list):
    """Calculates the Shannon Entropy of a list of bits."""
    if not bit_list:
        return 0.0
    counts = Counter(bit_list)
    entropy = 0.0
    total = len(bit_list)
    for bit, count in counts.items():
        probability = count / total
        if probability > 0:
            entropy -= probability * math.log2(probability)
    return entropy

def scan_headers(bam_path):
    """Checks the BAM header for known steganography signatures."""
    print("[*] Scanning BAM Headers...")
    suspicious_flags = []
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            header = bam.header.to_dict()
            if 'CO' in header:
                for comment in header['CO']:
                    if 'bamsteg' in comment:
                        suspicious_flags.append(f"CRITICAL: Found explicit bamsteg signature -> {comment}")
    except Exception as e:
        print(f"[!] Header read error: {e}")

    return suspicious_flags

def analyze_reads(bam_path, sample_size=10000):
    """Scans structural tags and performs localized LSB analysis per read."""
    print(f"[*] Analyzing up to {sample_size} reads for structural and statistical anomalies...")

    suspicious_tags = 0
    total_reads_scanned = 0
    high_entropy_reads = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            total_reads_scanned += 1

            # 1. Structural Scan (AUX Tags)
            if read.has_tag('ZS'):
                suspicious_tags += 1

            # 2. Localized Statistical Scan (Per-Read Entropy)
            if read.query_qualities is not None:
                read_lsbs = [q & 1 for q in read.query_qualities]
                read_entropy = calculate_shannon_entropy(read_lsbs)

                # A natural read rarely hits perfect 1.0 entropy unless it's pure noise.
                if read_entropy > 0.95:
                    high_entropy_reads += 1

            if total_reads_scanned >= sample_size:
                break

    return {
        "reads_scanned": total_reads_scanned,
        "suspicious_tags": suspicious_tags,
        "high_entropy_reads": high_entropy_reads
    }

def run_forensics(bam_path):
    print("=" * 60)
    print(f"  BAMSTEG DETECT - FORENSIC ANALYSIS")
    print(f"  Target: {bam_path}")
    print("=" * 60)

    header_alerts = scan_headers(bam_path)
    read_metrics = analyze_reads(bam_path)

    # --- THREAT SCORING LOGIC ---
    threat_level = "LOW (CLEAN)"
    color_code = "\033[92m" # Green

    alerts = []

    if len(header_alerts) > 0:
        threat_level = "CRITICAL"
        alerts.extend(header_alerts)

    if read_metrics["suspicious_tags"] > 0:
        threat_level = "HIGH"
        alerts.append(f"Found {read_metrics['suspicious_tags']} reads with suspicious 'ZS' carrier tags.")

    if read_metrics["high_entropy_reads"] > 0:
        if threat_level != "CRITICAL":
            threat_level = "HIGH"
        alerts.append(f"Found {read_metrics['high_entropy_reads']} reads exhibiting localized high entropy (>0.95). Indicates LSB ciphertext.")

    if threat_level == "HIGH":
        color_code = "\033[93m" # Yellow
    elif threat_level == "CRITICAL":
        color_code = "\033[91m" # Red

    reset_code = "\033[0m"

    print("\n--- RESULTS ---")
    print(f"Reads Analyzed       : {read_metrics['reads_scanned']}")
    print(f"High-Entropy Reads   : {read_metrics['high_entropy_reads']}")
    print(f"Suspicious AUX Tags  : {read_metrics['suspicious_tags']}")

    print(f"\nTHREAT LEVEL: {color_code}{threat_level}{reset_code}")

    if alerts:
        print("\nForensic Flags Raised:")
        for alert in alerts:
            print(f" -> {alert}")
    else:
        print("\nNo steganographic anomalies detected. File appears biologically standard.")
    print("=" * 60)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 bamsteg_detect.py <path_to_bam_file>")
        sys.exit(1)

    target_file = sys.argv[1]
    run_forensics(target_file)
