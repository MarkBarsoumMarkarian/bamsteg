import hashlib, random

def derive_seed(passphrase, key):
    h = hashlib.sha256((passphrase + key).encode()).digest()
    return int.from_bytes(h[:8], "little")

def select_reads(total_reads, n_needed, passphrase, key):
    if n_needed > total_reads:
        raise ValueError(f"Need {n_needed} reads but BAM only has {total_reads}.")
    seed = derive_seed(passphrase, key)
    return sorted(random.Random(seed).sample(range(total_reads), n_needed))

def bytes_to_bits(data):
    bits = []
    for byte in data:
        for j in range(7, -1, -1):
            bits.append((byte >> j) & 1)
    return bits

def bits_to_bytes(bits):
    while len(bits) % 8 != 0:
        bits.append(0)
    out = bytearray()
    for i in range(0, len(bits), 8):
        byte = 0
        for j in range(8):
            byte = (byte << 1) | (bits[i+j] & 1)
        out.append(byte)
    return bytes(out)
