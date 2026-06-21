#!/usr/bin/env python3
"""
bamsteg_detect.py -- Forensic Steganalysis Engine for BAM Files
v0.3 -- generic steganalysis, not just bamsteg-specific signatures

Detection methods:
  1. Header signature scan        -- catches naive/lazy embedders
  2. Generic suspicious aux scan  -- any X*/Y*/Z* tag with hex-pattern values
  3. Adaptive LSB chi-squared     -- head vs tail parity drift, instrument-independent
  4. RS-SPA (Sample Pairs Analysis) -- detects LSB embedding at low rates (<3% reads)
     where chi-squared loses power. Based on Fridrich et al. 2001.
  5. Per-score parity entropy     -- detects flattening of natural score clusters

RS-SPA background:
  In a clean quality score stream, adjacent score pairs (q_i, q_i+1) have a
  natural asymmetry: values differing by 1 occur less often than values
  differing by 2, because the instrument produces clustered distributions.
  LSB flipping destroys this asymmetry. SPA measures it via the ratio of
  same-parity pairs to cross-parity pairs across the stream.
  Unlike chi-squared it does not require a baseline calibration and is
  effective even when only 1-2% of bases are modified.
"""

import sys
import math
import pysam
from collections import Counter

CHI2_P005 = 3.841
CHI2_P001 = 6.635
SAMPLE_SIZE = 50_000
TAIL_FRAC = 0.15

# Generic aux tag scan: any two-letter tag starting with X, Y, Z
# (reserved for local use per SAM spec) with a suspicious hex-format string value
SUSPICIOUS_TAG_PREFIXES = ("X", "Y", "Z")
# bamsteg AUX pattern: exactly "XXXXXX:YY" where both parts are valid hex
def _is_hex_covert_pattern(val):
    if not isinstance(val, str):
        return False
    if len(val) == 9 and val[6] == ":":
        try:
            int(val[:6], 16); int(val[7:], 16)
            return True
        except ValueError:
            pass
    # also catch longer patterns: any string that's pure hex or hex:hex
    parts = val.split(":")
    if len(parts) == 2:
        try:
            int(parts[0], 16); int(parts[1], 16)
            return len(parts[0]) >= 4 and len(parts[1]) >= 2
        except ValueError:
            pass
    return False


def _entropy(counts, total):
    h = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            h -= p * math.log2(p)
    return h


# ---------------------------------------------------------------------------
# Test 1: header signature
# ---------------------------------------------------------------------------
def scan_headers(bam_path):
    alerts = []
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            hdr = bam.header.to_dict()
            for comment in hdr.get("CO", []):
                if "bamsteg" in comment.lower() or "steg" in comment.lower():
                    alerts.append(f"Steganography signature in @CO: '{comment}'")
            # also check for unusual @PG entries
            for pg in hdr.get("PG", []):
                if "steg" in str(pg).lower():
                    alerts.append(f"Steganography-related @PG entry: {pg}")
    except Exception as e:
        alerts.append(f"Header read error: {e}")
    return alerts


# ---------------------------------------------------------------------------
# Test 2: generic suspicious aux tag scan
# ---------------------------------------------------------------------------
def scan_aux_tags(bam_path, sample_size=SAMPLE_SIZE):
    """
    Scans all X*/Y*/Z* aux tags (SAM spec reserved-for-local-use namespace).
    Flags any that match a hex covert channel pattern.
    Does not require knowing the specific tag name (ZS, XV, etc.).
    """
    tag_counter = Counter()   # tag_name -> count
    covert_counter = Counter()  # tag_name -> count with suspicious pattern
    reads = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            reads += 1
            for tag, val in (read.get_tags() if read.get_tags() else []):
                if tag[0] in SUSPICIOUS_TAG_PREFIXES:
                    tag_counter[tag] += 1
                    if _is_hex_covert_pattern(val):
                        covert_counter[tag] += 1
            if reads >= sample_size:
                break

    return {
        "reads_scanned": reads,
        "suspicious_tags": dict(tag_counter),
        "covert_pattern_tags": dict(covert_counter),
        "total_covert_hits": sum(covert_counter.values()),
    }


# ---------------------------------------------------------------------------
# Test 3: adaptive LSB chi-squared
# ---------------------------------------------------------------------------
def scan_lsb_distribution(bam_path, sample_size=SAMPLE_SIZE):
    head_lsb = Counter()
    tail_lsb = Counter()
    score_parity = Counter()
    reads = 0
    total_bases = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            reads += 1
            quals = read.query_qualities
            if quals is None:
                continue
            n = len(quals)
            tail_start = int(n * (1 - TAIL_FRAC))
            for i, q in enumerate(quals):
                lsb = q & 1
                total_bases += 1
                score_parity[(q >> 1, lsb)] += 1
                if i < tail_start:
                    head_lsb[lsb] += 1
                else:
                    tail_lsb[lsb] += 1
            if reads >= sample_size:
                break

    if total_bases == 0 or sum(tail_lsb.values()) == 0:
        return None

    tail_total = sum(tail_lsb.values())
    head_total = sum(head_lsb.values())
    tail_lsb0_rate = tail_lsb[0] / tail_total

    exp_lsb0 = head_total * tail_lsb0_rate
    exp_lsb1 = head_total * (1 - tail_lsb0_rate)
    chi2 = 0.0
    if exp_lsb0 > 0:
        chi2 += (head_lsb[0] - exp_lsb0) ** 2 / exp_lsb0
    if exp_lsb1 > 0:
        chi2 += (head_lsb[1] - exp_lsb1) ** 2 / exp_lsb1

    bucket_entropies = []
    for score_half in range(5, 30):
        c0 = score_parity.get((score_half, 0), 0)
        c1 = score_parity.get((score_half, 1), 0)
        bt = c0 + c1
        if bt >= 50:
            bucket_entropies.append(_entropy({0: c0, 1: c1}, bt))

    mean_bucket_h = (sum(bucket_entropies) / len(bucket_entropies)
                     if bucket_entropies else 0.0)

    return {
        "total_bases": total_bases,
        "reads_scanned": reads,
        "tail_lsb0_rate": round(tail_lsb0_rate, 5),
        "head_lsb0_rate": round(head_lsb[0] / head_total, 5),
        "chi2_statistic": round(chi2, 4),
        "chi2_p005": chi2 > CHI2_P005,
        "chi2_p001": chi2 > CHI2_P001,
        "mean_bucket_entropy": round(mean_bucket_h, 5),
        "buckets_checked": len(bucket_entropies),
    }


# ---------------------------------------------------------------------------
# Test 4: RS-SPA (Sample Pairs Analysis)
# ---------------------------------------------------------------------------
def scan_spa(bam_path, sample_size=SAMPLE_SIZE):
    """
    Sample Pairs Analysis for LSB steganography detection.

    For a sequence of quality scores q_0, q_1, ..., q_n consider pairs
    (q_i, q_i+1). Define:
      P_m = number of pairs where |q_i - q_{i+1}| == m

    In a natural quality score stream, P_1 < P_2 because the instrument
    produces smooth distributions. LSB flipping randomizes the LSB,
    which specifically inflates P_1 relative to P_2 (adjacent values
    that differ only in their LSB become more common).

    The SPA statistic is: R = P_1 / P_2
    Clean file:   R < 1.0 (P_1 < P_2 naturally)
    Embedded file: R approaches or exceeds 1.0

    We also compute the asymmetry coefficient:
      A = (C_even - C_odd) / total_pairs
    where C_even = pairs where both scores have same parity,
          C_odd  = pairs where scores have different parity.
    In clean data A deviates from 0; after LSB embedding A -> 0.
    """
    p_counts = Counter()   # |delta| -> count
    same_parity = 0
    diff_parity = 0
    total_pairs = 0
    reads = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            reads += 1
            quals = read.query_qualities
            if quals is None or len(quals) < 2:
                continue
            for i in range(len(quals) - 1):
                q0, q1 = int(quals[i]), int(quals[i+1])
                delta = abs(q0 - q1)
                p_counts[delta] += 1
                total_pairs += 1
                if (q0 & 1) == (q1 & 1):
                    same_parity += 1
                else:
                    diff_parity += 1
            if reads >= sample_size:
                break

    if total_pairs == 0 or p_counts[2] == 0:
        return None

    r_statistic = p_counts[1] / p_counts[2]
    asymmetry = (same_parity - diff_parity) / total_pairs

    # In clean Illumina data: R is typically 0.3-0.7 (P1 << P2)
    # After LSB embedding:   R climbs toward and above 1.0
    # Asymmetry: clean data has A != 0; embedded data has A -> 0
    return {
        "total_pairs": total_pairs,
        "p1": p_counts[1],
        "p2": p_counts[2],
        "r_statistic": round(r_statistic, 5),
        "r_suspicious": r_statistic > 0.85,   # within 15% of 1.0
        "r_high": r_statistic > 1.0,
        "asymmetry_coeff": round(asymmetry, 5),
        "asymmetry_flat": abs(asymmetry) < 0.02,  # near-zero = suspicious
    }


# ---------------------------------------------------------------------------
# Threat scoring
# ---------------------------------------------------------------------------
def _score_threat(header_alerts, aux, lsb, spa):
    alerts = []
    level = "CLEAN"
    rank = ["CLEAN", "SUSPICIOUS", "HIGH", "CRITICAL"]

    def up(l):
        nonlocal level
        if rank.index(l) > rank.index(level):
            level = l

    # Header
    for a in header_alerts:
        alerts.append(f"[HEADER] {a}"); up("CRITICAL")

    # AUX
    if aux["total_covert_hits"] > 0:
        tags = ", ".join(
            f"{t}:{n}" for t, n in aux["covert_pattern_tags"].items()
        )
        alerts.append(
            f"[AUX] {aux['total_covert_hits']} reads carry hex-pattern "
            f"covert channel tags: {tags}."
        ); up("CRITICAL")
    elif aux["suspicious_tags"]:
        tags = ", ".join(
            f"{t}:{n}" for t, n in aux["suspicious_tags"].items()
        )
        alerts.append(
            f"[AUX] Non-standard X*/Y*/Z* tags present: {tags} "
            f"(no covert pattern matched, may be benign)."
        ); up("SUSPICIOUS")

    # Chi-squared
    if lsb:
        chi2 = lsb["chi2_statistic"]
        drift = abs(lsb["head_lsb0_rate"] - lsb["tail_lsb0_rate"])
        if lsb["chi2_p001"]:
            alerts.append(
                f"[LSB-CHI2] Statistic={chi2:.2f} >> {CHI2_P001} (p<0.001). "
                f"Head/tail parity drift={drift:.4f}."
            ); up("HIGH")
        elif lsb["chi2_p005"]:
            alerts.append(
                f"[LSB-CHI2] Statistic={chi2:.2f} > {CHI2_P005} (p<0.05). "
                f"Borderline parity anomaly (drift={drift:.4f})."
            ); up("SUSPICIOUS")
        if lsb["mean_bucket_entropy"] > 0.90 and lsb["buckets_checked"] >= 3:
            alerts.append(
                f"[LSB-ENTROPY] Per-score-bucket entropy={lsb['mean_bucket_entropy']:.4f} "
                f"(>0.90 indicates LSB ciphertext noise)."
            ); up("SUSPICIOUS")

    # SPA
    if spa:
        if spa["r_high"]:
            alerts.append(
                f"[SPA] R-statistic={spa['r_statistic']:.4f} > 1.0. "
                f"P1/P2 pair ratio inverted -- strong LSB embedding signal. "
                f"Parity asymmetry={spa['asymmetry_coeff']:.4f}."
            ); up("HIGH")
        elif spa["r_suspicious"]:
            alerts.append(
                f"[SPA] R-statistic={spa['r_statistic']:.4f} > 0.85. "
                f"Approaching P1/P2 inversion -- possible low-rate LSB embedding."
            ); up("SUSPICIOUS")
        if spa["asymmetry_flat"] and not spa["r_suspicious"]:
            alerts.append(
                f"[SPA] Parity asymmetry={spa['asymmetry_coeff']:.4f} near zero. "
                f"Natural score structure flattened -- consistent with LSB manipulation."
            ); up("SUSPICIOUS")

    return level, alerts


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run_forensics(bam_path, sample_size=SAMPLE_SIZE):
    COLOR = {"CLEAN": "\033[92m", "SUSPICIOUS": "\033[93m",
             "HIGH": "\033[91m", "CRITICAL": "\033[95m"}
    RESET = "\033[0m"

    print("=" * 66)
    print("  BAMSTEG-DETECT  --  Genomic Steganalysis Engine v0.3")
    print(f"  Target : {bam_path}")
    print(f"  Sample : up to {sample_size:,} reads")
    print("=" * 66)

    print("[1/4] Header scan ...")
    h_alerts = scan_headers(bam_path)

    print("[2/4] Generic aux tag scan ...")
    aux = scan_aux_tags(bam_path, sample_size)

    print("[3/4] Adaptive LSB chi-squared ...")
    lsb = scan_lsb_distribution(bam_path, sample_size)

    print("[4/4] Sample Pairs Analysis (RS-SPA) ...")
    spa = scan_spa(bam_path, sample_size)

    level, alerts = _score_threat(h_alerts, aux, lsb, spa)
    color = COLOR.get(level, "")

    print()
    print("--- METRICS ---")
    print(f"  Reads scanned      : {aux['reads_scanned']:,}")
    if lsb:
        print(f"  Bases analyzed     : {lsb['total_bases']:,}")
        print(f"  Tail LSB=0 rate    : {lsb['tail_lsb0_rate']:.5f}  (instrument baseline)")
        print(f"  Head LSB=0 rate    : {lsb['head_lsb0_rate']:.5f}  (payload zone)")
        print(f"  Chi-squared        : {lsb['chi2_statistic']:.4f}  "
              f"(p<0.05={CHI2_P005}, p<0.001={CHI2_P001})")
        print(f"  Bucket entropy     : {lsb['mean_bucket_entropy']:.5f}  "
              f"over {lsb['buckets_checked']} score buckets")
    if spa:
        print(f"  SPA R-statistic    : {spa['r_statistic']:.5f}  "
              f"(clean<0.85, suspicious>0.85, high>1.0)")
        print(f"  Parity asymmetry   : {spa['asymmetry_coeff']:.5f}  "
              f"(near-zero = suspicious)")
        print(f"  Pair counts        : P1={spa['p1']:,}  P2={spa['p2']:,}")
    covert = aux["total_covert_hits"]
    print(f"  Covert AUX tags    : {covert} hex-pattern hits "
          f"across {len(aux['suspicious_tags'])} tag type(s)")
    print()
    print(f"  THREAT LEVEL: {color}{level}{RESET}")
    print()

    if alerts:
        print("--- FORENSIC FLAGS ---")
        for a in alerts:
            print(f"  [!] {a}")
    else:
        print("  No anomalies detected. File appears clean.")
    print("=" * 66)
    return level, alerts, {"lsb": lsb, "spa": spa, "aux": aux}


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: bamsteg_detect.py <file.bam> [--sample N]")
        sys.exit(1)
    path = sys.argv[1]
    n = SAMPLE_SIZE
    for i, arg in enumerate(sys.argv):
        if arg == "--sample" and i + 1 < len(sys.argv):
            n = int(sys.argv[i + 1])
    run_forensics(path, sample_size=n)
