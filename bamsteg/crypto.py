import os, struct
from Crypto.Cipher import AES
from Crypto.Protocol.KDF import PBKDF2
from Crypto.Hash import SHA256, HMAC

SALT_LEN = 16; NONCE_LEN = 12; TAG_LEN = 16; PBKDF2_ITER = 200_000

def derive_key(passphrase, salt):
    return PBKDF2(passphrase.encode(), salt, dkLen=32, count=PBKDF2_ITER,
                  prf=lambda p, s: HMAC.new(p, s, SHA256).digest())

def encrypt(payload, passphrase):
    salt = os.urandom(SALT_LEN)
    key = derive_key(passphrase, salt)
    nonce = os.urandom(NONCE_LEN)
    cipher = AES.new(key, AES.MODE_GCM, nonce=nonce)
    ct, tag = cipher.encrypt_and_digest(payload)
    return salt + nonce + tag + struct.pack("<I", len(ct)) + ct

def decrypt(blob, passphrase):
    o = 0
    salt = blob[o:o+SALT_LEN]; o += SALT_LEN
    nonce = blob[o:o+NONCE_LEN]; o += NONCE_LEN
    tag = blob[o:o+TAG_LEN]; o += TAG_LEN
    ct_len = struct.unpack("<I", blob[o:o+4])[0]; o += 4
    ct = blob[o:o+ct_len]
    cipher = AES.new(derive_key(passphrase, salt), AES.MODE_GCM, nonce=nonce)
    try:
        return cipher.decrypt_and_verify(ct, tag)
    except ValueError:
        raise ValueError("Decryption failed: wrong passphrase or corrupted payload.")
