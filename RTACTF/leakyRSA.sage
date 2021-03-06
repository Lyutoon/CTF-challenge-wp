from pwn import remote
import requests
from pwnlib.tubes.tube import *
from hashlib import *
from Crypto.Util.number import *
from sage.combinat.symmetric_group_algebra import epsilon
from tqdm import tqdm
import gmpy2
import random
import math
from Crypto.Hash import SHA256
from Crypto.Cipher import AES
from factordb.factordb import FactorDB
from sage.modules.free_module_integer import IntegerLattice
import itertools
from random import shuffle

# r = remote('hiyoko.quals.seccon.jp', '10042')
# context(log_level='debug')
# ALPHABET = string.ascii_letters + string.digits

# rec = r.recvline().decode()
# print(rec)
# suffix = rec[rec.find('+')+2:rec.find(')')]
# digest = rec[rec.find('==')+3:-1]
# print(f"suffix: {suffix} \ndigest: {digest}")

# for i in itertools.product(ALPHABET, repeat=4):
#     prefix = ''.join(i)
#     guess = prefix + suffix
#     if sha256(guess.encode()).hexdigest() == digest:
#         # log.info(f"Find XXXX: {prefix}")
#         print((f"Find XXXX: {prefix}"))
#         break
# r.sendline(prefix.encode())

# r.interactive()

n = 0x2ac1fcbbf63ffeade11cd2c57c37db18d96d52e433bd9034d4eac2c269ea49a81e5ac41fb631523bb5983adc6fc939c073c13d8a3a42a06accf5a9c304fc444508a8b5833b5431e9af7007bb216c510c62a97eb1fe380bf155b3e497c7d70c2bb921f97eec61e9e9ac7b5d71e47876d20cbfb1a0732e29ec6872041eb67e0ccd39d7b6429bda1581537dda95e79d3aad4df072beada1c72a4ffd86db91918ec9db44ab9c4bebf387ccc1ce7b2540b0d595a4c11823cdbcd8850bd3b666b4a08bd69de515afecc75b283ae47fbf3af6f034f3b0f7848dec935ba8b97e36d2d0a9208df63610cf8825fb729aacaa4c119d0b4c5e230e080d7633f145d22eb06b917fe632c01a373b1c4c8a741bea1d5dd98003e9
e = 0x10001
c = 0xe50a858f715238a9ab44dfd691f6e5ced84e74115e003e31a98b324cf9e8bc9cfe08065f2538cff519e035566b4080742139062e672a0ad3196275cb121ea837de2808f99958bcfe58d1c8996f291412220d01fe65fbf18611b348407b2e2db45b2adcc341926c6d76a9d08fc77db0fedef78cec9e4b881812e60c015c1005dfd0b9408cb3c6f9f98332f165acc3ae98ef97f2a1d98524fe240d3351676ed84ddb73283a6d3efc40bbd466fe3532e579eb9adf07ebbc49af71fb22934a75a69a538eca0fd4e2a5b617abb361a64c553985950dd5201ac7c631580c8bb27d795a196d584ae7c7478bdc1b5ff531ff88e984bceb1e26cf9793f99a11287555d5d2d2a13e1171f77bf8491d8dfa297e9cd6d4b7d
s = 0x14c0af71a961be72e2d4e5ee06337cde2034db1d920f4476e3a3371c8a35e7ba8efabf5c8e8ff86e7297156c4fde5bdc7aabe1516a46c554236104022eb4544f1d7fcb80279595dfe0527bcc373909ce7cc0965ece5ff76b7ee9a5cc31a1b567ed3ddd2364bb596e3c41e4fffb5974f71e788da5c21598e9c6dc32bca162026ba3c410bb1c5c9d5bed4c3b97e3cacbd7b6693f29c74b0756381b658efaa757d448f62a48fbdb06604525222aa51797a1a1e43af4b0c221deef47f84bb5bfa1480cd31242c3d7fba21bdf487709853879dcea284e44cb5ee1a02c558a29740a44e39c7ee3a97ab4805d21cb90b596bd86c51f4e0f783701da73c66f5a4c67d989bb2dba2b8a55a697eb187cc181fc8ce54de370

# P.<x> = PolynomialRing(Zmod(n))

# f = 20211219 + x*s
# f = f.monic()

# print(f.small_roots(X=2**500, beta=0.4, epsilon=5/32))

q = 1922170375642253101720795198135909686494155409464023438804428444232467343105956174957087101662812915933127607975743244342906562612795777257605243515321
mod = q^2

d = inverse(e, q*(q-1))
print(long_to_bytes(int(pow(c, d, mod))))