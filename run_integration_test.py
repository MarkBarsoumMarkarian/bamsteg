#!/usr/bin/env python3
import os
import sys
import pysam

# Ensure the local directory takes import priority
sys.path.insert(0, '.')
from bamsteg.embed_aux import embed_aux, extract_aux
from bamsteg.embed_lsb import embed_lsb, extract_lsb

def generate_sample_bam(filepath, num_reads=500, read_length=100):
    """Generates a valid, clean sample BAM file for localized testing."""
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'SN': 'chr1', 'LN': 10000}]
    }
    with pysam.AlignmentFile(filepath, "wb", header=header) as out:
        for i in range(num_reads):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i}"
            a.query_sequence = "A" * read_length
            a.query_qualities = [30] * read_length  # Flat Q30 baseline
            a.reference_id = 0
            a.reference_start = i * 10
            a.cigar = ((0, read_length),)
            out.write(a)
    print(f"[+] Generated baseline sample file: {filepath}")

def main():
    # Setup paths inside your current workspace
    workspace = "./test_sandbox"
    os.makedirs(workspace, exist_ok=True)

    input_bam = os.path.join(workspace, "baseline.bam")
    aux_out   = os.path.join(workspace, "output_aux.bam")
    lsb_out   = os.path.join(workspace, "output_lsb.bam")

    passphrase = "MelinoeOS_Crypto_Pass_2026"

    # Generate the baseline 500-read file
    generate_sample_bam(input_bam)

    print("\n" + "="*50)
    print("  RUNNING BAMSTEG INTEGRATION SUITE")
    print("="*50)

    # --- AUX MODE TEST ---
    payload_aux = b"AUX-PATIENT-ID-WATERMARK-99X"
    print("[*] Testing AUX Mode...")
    stats_aux = embed_aux(input_bam, aux_out, payload_aux, passphrase)
    recovered_aux = extract_aux(aux_out, passphrase)
    aux_pass = recovered_aux == payload_aux
    print(f"    Round-trip: {'PASS' if aux_pass else 'FAIL'}")
    print(f"    Metrics: {stats_aux}")

    # --- LSB MODE TEST ---
    payload_lsb = b"LSB-COVERT-CHANNEL-DATA-007"
    print("\n[*] Testing LSB Mode...")
    stats_lsb = embed_lsb(input_bam, lsb_out, payload_lsb, passphrase)
    recovered_lsb = extract_lsb(lsb_out, passphrase)
    lsb_pass = recovered_lsb == payload_lsb
    print(f"    Round-trip: {'PASS' if lsb_pass else 'FAIL'}")
    print(f"    Metrics: {stats_lsb}")

    # --- QUALITY CHECK ---
    with pysam.AlignmentFile(input_bam, 'rb') as fi, pysam.AlignmentFile(lsb_out, 'rb') as fo:
        ri = list(fi.fetch(until_eof=True))
        ro = list(fo.fetch(until_eof=True))
    max_d = max(abs(int(a)-int(b)) for r1,r2 in zip(ri,ro)
                if r1.query_qualities and r2.query_qualities
                for a,b in zip(r1.query_qualities, r2.query_qualities))
    print(f"\n[+] Maximum LSB Phred score variance: {max_d} (Expected <= 1)")

    print("="*50)
    if aux_pass and lsb_pass and max_d <= 1:
        print("  ALL SYSTEM CHECKS: PASSED")
    else:
        print("  SYSTEM VERIFICATION: FAILED")
    print("="*50)

if __name__ == "__main__":
    main()
