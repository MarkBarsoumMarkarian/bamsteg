import struct
from reedsolo import RSCodec, ReedSolomonError

ECC_SYMBOLS = 32
CHUNK_SIZE = 200

def encode(data, ecc_symbols=ECC_SYMBOLS):
    rsc = RSCodec(ecc_symbols)
    chunks = [data[i:i+CHUNK_SIZE] for i in range(0, len(data), CHUNK_SIZE)]
    out = []
    for chunk in chunks:
        enc = bytes(rsc.encode(chunk))
        out.append(struct.pack("<H", len(enc)) + enc)
    return struct.pack("<II", len(chunks), ecc_symbols) + b"".join(out)

def decode(data):
    num_chunks, ecc_symbols = struct.unpack("<II", data[:8])
    rsc = RSCodec(ecc_symbols)
    offset = 8
    result = bytearray()
    for _ in range(num_chunks):
        clen = struct.unpack("<H", data[offset:offset+2])[0]; offset += 2
        chunk = data[offset:offset+clen]; offset += clen
        decoded, _, _ = rsc.decode(chunk)
        result.extend(bytes(decoded))
    return bytes(result)
