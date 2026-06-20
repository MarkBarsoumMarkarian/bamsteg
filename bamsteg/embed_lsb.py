"""
LSB-mode embedder/extractor for bamsteg.
Hides payload bits in the LSB of Phred quality scores.

Architecture fix: we use a fixed 4-read PRELUDE to store the blob byte count,
then the full PRNG-selected set for the actual blob. This avoids the n-dependent
PRNG selection problem where sample(n=4) and sample(n=34) from same seed differ.

Prelude reads are selected with a separate SEED_KEY so they don't overlap
with payload-carrier reads.
"""

import struct
import pysam
from .crypto import encrypt, decrypt
from .ecc import encode as ecc_encode, decode as ecc_decode
from .selector import select_reads, bytes_to_bits, bits_to_bytes

MAGIC = b"\xBA\xA5\x1E\x61"
VERSION = 1
BITS_PER_BASE = 1

# The prelude stores: MAGIC(4) + VERSION(1) + RESERVED(2) + BLOB_LEN(4) = 11 bytes = 88 bits
# At 32 bits per read (32bp reads), we need 3 reads for the prelude.
# We use PRELUDE_READS=4 for safety margin.
PRELUDE_BYTES = len(MAGIC) + 1 + 2 + 4  # = 11
PRELUDE_READS = 4  # fixed, always 4 reads, enough for 128 bits at any read length >= 11bp

SEED_PRELUDE = "lsb_prelude"
SEED_PAYLOAD = "lsb_payload"

def _pack_prelude(blob_len: int) -> bytes:
    """11-byte prelude: magic + version + reserved + blob_len."""
    return MAGIC + struct.pack("<BHI", VERSION, 0, blob_len)

def _unpack_prelude(data: bytes):
    """Parse prelude, return (version, blob_len)."""
    if not data[:4] == MAGIC:
        raise ValueError(f"LSB magic mismatch. Got {data[:4].hex()}, expected {MAGIC.hex()}")
    version, _reserved, blob_len = struct.unpack("<BHI", data[4:11])
    return version, blob_len

def embed_lsb(
    input_bam: str,
    output_bam: str,
    payload: bytes,
    passphrase: str,
    bits_per_base: int = BITS_PER_BASE,
) -> dict:
    encrypted = encrypt(payload, passphrase)
    blob = ecc_encode(encrypted)
    blob_bits = bytes_to_bits(blob)
    n_blob_bits = len(blob_bits)

    prelude_bytes = _pack_prelude(len(blob))
    prelude_bits = bytes_to_bits(prelude_bytes)

    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        hdr = bam_in.header.to_dict()
        hdr.setdefault("CO", [])
        hdr["CO"].append("bamsteg:lsb:v1")
        all_reads = list(bam_in.fetch(until_eof=True))

    total_reads = len(all_reads)
    sample = all_reads[:1000]
    avg_len = int(sum(r.query_length or 0 for r in sample) / max(len(sample), 1))
    bits_per_read = avg_len * bits_per_base

    if bits_per_read == 0:
        raise ValueError("Zero-length reads.")

    # prelude: fixed PRELUDE_READS reads
    prelude_indices = select_reads(total_reads, PRELUDE_READS, passphrase, SEED_PRELUDE)

    # payload: exactly as many reads as needed
    n_payload_reads = (n_blob_bits + bits_per_read - 1) // bits_per_read
    payload_indices = select_reads(total_reads, n_payload_reads, passphrase, SEED_PAYLOAD)

    # build per-read bit assignments
    prelude_map: dict = {}
    bit_cursor = 0
    for idx in prelude_indices:
        rlen = (all_reads[idx].query_length or 0) * bits_per_base
        chunk = prelude_bits[bit_cursor: bit_cursor + rlen]
        prelude_map[idx] = chunk
        bit_cursor += len(chunk)
        if bit_cursor >= len(prelude_bits):
            break

    payload_map: dict = {}
    bit_cursor = 0
    for idx in payload_indices:
        rlen = (all_reads[idx].query_length or 0) * bits_per_base
        chunk = blob_bits[bit_cursor: bit_cursor + rlen]
        payload_map[idx] = chunk
        bit_cursor += len(chunk)
        if bit_cursor >= n_blob_bits:
            break

    with pysam.AlignmentFile(
        output_bam, "wb", header=pysam.AlignmentHeader.from_dict(hdr)
    ) as bam_out:
        for i, read in enumerate(all_reads):
            modifications = []
            if i in prelude_map:
                modifications = prelude_map[i]
            if i in payload_map:
                modifications = payload_map[i]  # payload wins if overlap (shouldn't happen)

            if modifications and read.query_qualities is not None:
                quals = list(read.query_qualities)
                for j, bit in enumerate(modifications):
                    if j < len(quals):
                        quals[j] = (quals[j] & 0xFE) | (bit & 1)
                read.query_qualities = quals

            bam_out.write(read)

    return {
        "mode": "lsb",
        "payload_bytes": len(payload),
        "blob_bytes": len(blob),
        "prelude_reads": PRELUDE_READS,
        "payload_reads": n_payload_reads,
        "total_reads": total_reads,
        "reads_used_pct": round(100 * (PRELUDE_READS + n_payload_reads) / total_reads, 3),
    }

def extract_lsb(
    input_bam: str,
    passphrase: str,
    bits_per_base: int = BITS_PER_BASE,
) -> bytes:
    """
    Extract LSB payload. No expected_payload_bytes needed anymore --
    blob_len is encoded in the prelude, and prelude always uses PRELUDE_READS reads.
    """
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        all_reads = list(bam_in.fetch(until_eof=True))

    total_reads = len(all_reads)
    sample = all_reads[:1000]
    avg_len = int(sum(r.query_length or 0 for r in sample) / max(len(sample), 1))
    bits_per_read = avg_len * bits_per_base

    # Step 1: read prelude (fixed PRELUDE_READS reads)
    prelude_indices = select_reads(total_reads, PRELUDE_READS, passphrase, SEED_PRELUDE)
    prelude_bits = []
    for idx in prelude_indices:
        r = all_reads[idx]
        if r.query_qualities is None:
            continue
        for q in r.query_qualities:
            prelude_bits.append(q & 1)
            if len(prelude_bits) >= PRELUDE_BYTES * 8:
                break
        if len(prelude_bits) >= PRELUDE_BYTES * 8:
            break

    prelude_bytes = bits_to_bytes(prelude_bits[: PRELUDE_BYTES * 8])
    _version, blob_len = _unpack_prelude(prelude_bytes)

    # Step 2: select payload reads using exact blob_len
    blob_bits_needed = blob_len * 8
    n_payload_reads = (blob_bits_needed + bits_per_read - 1) // bits_per_read
    payload_indices = select_reads(total_reads, n_payload_reads, passphrase, SEED_PAYLOAD)

    collected_bits = []
    for idx in payload_indices:
        r = all_reads[idx]
        if r.query_qualities is None:
            continue
        for q in r.query_qualities:
            collected_bits.append(q & 1)
            if len(collected_bits) >= blob_bits_needed:
                break
        if len(collected_bits) >= blob_bits_needed:
            break

    blob = bits_to_bytes(collected_bits[:blob_bits_needed])
    encrypted = ecc_decode(blob)
    return decrypt(encrypted, passphrase)
