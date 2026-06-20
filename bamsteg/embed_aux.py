import os, struct, pysam
from .crypto import encrypt, decrypt
from .ecc import encode as ecc_encode, decode as ecc_decode
from .selector import select_reads, bytes_to_bits, bits_to_bytes

AUX_TAG = "ZS"
MAGIC = b"\xBA\xA5\x1E\x60"
VERSION = 1

def _build_blob(payload, passphrase):
    encrypted = encrypt(payload, passphrase)
    encoded = ecc_encode(encrypted)
    return MAGIC + struct.pack("<BH", VERSION, 0) + encoded

def _unpack_blob(blob, passphrase):
    if not blob.startswith(MAGIC):
        raise ValueError("No bamsteg magic. Wrong read or not embedded.")
    encrypted = ecc_decode(blob[len(MAGIC)+3:])
    return decrypt(encrypted, passphrase)

def embed_aux(input_bam, output_bam, payload, passphrase):
    blob = _build_blob(payload, passphrase)
    bam_basename = os.path.basename(input_bam)
    with pysam.AlignmentFile(input_bam, "rb") as f:
        hdr = f.header.to_dict()
        hdr.setdefault("CO", []).append("bamsteg:aux:v1")
        all_reads = list(f.fetch(until_eof=True))
    total = len(all_reads)
    selected = select_reads(total, len(blob), passphrase, bam_basename)
    idx_to_pos = {idx: pos for pos, idx in enumerate(selected)}
    with pysam.AlignmentFile(output_bam, "wb",
            header=pysam.AlignmentHeader.from_dict(hdr)) as out:
        for i, read in enumerate(all_reads):
            if i in idx_to_pos:
                pos = idx_to_pos[i]
                read.set_tag(AUX_TAG, f"{pos:06x}:{blob[pos]:02x}", value_type="Z")
            out.write(read)
    return {"mode": "aux", "payload_bytes": len(payload), "blob_bytes": len(blob),
            "reads_used": len(blob), "total_reads": total,
            "reads_used_pct": round(100*len(blob)/total, 3)}

def extract_aux(input_bam, passphrase):
    with pysam.AlignmentFile(input_bam, "rb") as f:
        all_reads = list(f.fetch(until_eof=True))
    collected = {}
    for read in all_reads:
        if read.has_tag(AUX_TAG):
            try:
                pos_hex, byte_hex = read.get_tag(AUX_TAG).split(":")
                collected[int(pos_hex, 16)] = int(byte_hex, 16)
            except (ValueError, AttributeError):
                continue
    if not collected:
        raise ValueError("No bamsteg ZS tags found.")
    blob = bytes([collected.get(i, 0) for i in range(max(collected)+1)])
    return _unpack_blob(blob, passphrase)
