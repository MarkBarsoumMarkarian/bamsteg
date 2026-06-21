#!/usr/bin/env python3
import os, sys
import pysam

sys.path.insert(0, '.')
from bamsteg.embed_aux import embed_aux, extract_aux
from bamsteg.embed_lsb import embed_lsb, extract_lsb
from bamsteg.bamsteg_detect import run_forensics

PASS = "\033[92mPASS\033[0m"
FAIL = "\033[91mFAIL\033[0m"
results = []

def check(label, condition):
    results.append((label, condition))
    print(f"  {'[+]' if condition else '[!]'} {label}: {PASS if condition else FAIL}")

def generate_sample_bam(filepath, num_reads=500, read_length=100):
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'SN': 'chr1', 'LN': 10000}]
    }
    with pysam.AlignmentFile(filepath, "wb", header=header) as out:
        for i in range(num_reads):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i}"
            a.query_sequence = "A" * read_length
            a.query_qualities = [30] * read_length
            a.reference_id = 0
            a.reference_start = i * 10
            a.cigar = ((0, read_length),)
            out.write(a)
    print(f"[+] Generated baseline: {filepath}")

def main():
    workspace = "./test_sandbox"
    os.makedirs(workspace, exist_ok=True)

    input_bam = os.path.join(workspace, "baseline.bam")
    aux_out   = os.path.join(workspace, "output_aux.bam")
    lsb_out   = os.path.join(workspace, "output_lsb.bam")
    passphrase = "MelinoeOS_Crypto_Pass_2026"

    generate_sample_bam(input_bam)

    print("\n" + "="*56)
    print("  BAMSTEG INTEGRATION SUITE")
    print("="*56)

    # ------------------------------------------------------------------
    # 1. AUX mode embed/extract round-trip
    # ------------------------------------------------------------------
    print("\n[*] AUX mode...")
    payload_aux = b"AUX-PATIENT-ID-WATERMARK-99X"
    stats_aux = embed_aux(input_bam, aux_out, payload_aux, passphrase)
    recovered_aux = extract_aux(aux_out, passphrase)
    check("AUX round-trip", recovered_aux == payload_aux)

    # biological invariance: non-carrier reads unchanged
    with pysam.AlignmentFile(input_bam, 'rb') as fi, \
         pysam.AlignmentFile(aux_out, 'rb') as fo:
        ri = list(fi.fetch(until_eof=True))
        ro = list(fo.fetch(until_eof=True))
    dirty = sum(1 for a, b in zip(ri, ro)
                if not b.has_tag('ZS') and
                list(a.query_qualities) != list(b.query_qualities))
    check("AUX biological invariance (0 non-carrier reads modified)", dirty == 0)

    # wrong passphrase must raise
    try:
        extract_aux(aux_out, "wrong")
        check("AUX wrong passphrase rejected", False)
    except Exception:
        check("AUX wrong passphrase rejected", True)

    # ------------------------------------------------------------------
    # 2. LSB mode embed/extract round-trip
    # ------------------------------------------------------------------
    print("\n[*] LSB mode...")
    payload_lsb = b"LSB-COVERT-CHANNEL-DATA-007"
    stats_lsb = embed_lsb(input_bam, lsb_out, payload_lsb, passphrase)
    recovered_lsb = extract_lsb(lsb_out, passphrase)
    check("LSB round-trip", recovered_lsb == payload_lsb)

    with pysam.AlignmentFile(input_bam, 'rb') as fi, \
         pysam.AlignmentFile(lsb_out, 'rb') as fo:
        ri = list(fi.fetch(until_eof=True))
        ro = list(fo.fetch(until_eof=True))
    max_d = max(abs(int(a) - int(b))
                for r1, r2 in zip(ri, ro)
                if r1.query_qualities and r2.query_qualities
                for a, b in zip(r1.query_qualities, r2.query_qualities))
    check(f"LSB max quality delta <= 1 (got {max_d})", max_d <= 1)

    try:
        extract_lsb(lsb_out, "wrong")
        check("LSB wrong passphrase rejected", False)
    except Exception:
        check("LSB wrong passphrase rejected", True)

    # ------------------------------------------------------------------
    # 3. Detector assertions
    # ------------------------------------------------------------------
    print("\n[*] Detector...")

    level_clean, _, _ = run_forensics(input_bam, sample_size=500)
    check("Detector: clean baseline -> CLEAN", level_clean == "CLEAN")

    level_aux, alerts_aux, _ = run_forensics(aux_out, sample_size=500)
    check("Detector: AUX embedded -> CRITICAL",
          level_aux == "CRITICAL")
    check("Detector: AUX flags covert AUX tags",
          any("AUX" in a for a in alerts_aux))

    level_lsb, alerts_lsb, _ = run_forensics(lsb_out, sample_size=500)
    check("Detector: LSB embedded -> CRITICAL or HIGH",
          level_lsb in ("CRITICAL", "HIGH"))
    check("Detector: LSB flags header or LSB anomaly",
          any("HEADER" in a or "LSB" in a for a in alerts_lsb))

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "="*56)
    passed = sum(1 for _, ok in results if ok)
    total  = len(results)
    print(f"  {passed}/{total} checks passed")
    if passed == total:
        print(f"  ALL CHECKS: {PASS}")
    else:
        print(f"  SOME CHECKS: {FAIL}")
        for label, ok in results:
            if not ok:
                print(f"    FAILED: {label}")
    print("="*56)
    sys.exit(0 if passed == total else 1)

if __name__ == "__main__":
    main()
