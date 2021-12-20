from pwn import *
import string
import base64
import math
from libnum import *
import gmpy2
import os
# import random
# from libnum import crt
from tqdm import tqdm
from hashlib import sha256, md5
from Crypto.Hash import SHA256
from Crypto.PublicKey import DSA
from Crypto.Cipher import AES
from itertools import product
from sage.all import *
from Crypto.Util.number import *
import randcrack
import random
from sm4 import SM4Key

def fermat(n):
    a = isqrt(n)
    b2 = a * a - n
    b = isqrt(n)
    count = 0
    while b * b != b2:
        a = a + 1
        b2 = a * a - n
        b = isqrt(b2)
        count += 1
    p = a + b
    q = a - b
    assert n == p * q
    return p, q

n = 0xe72988e811f04091c3291ac28f1e8332193187f3dc5af01579c36badb06671aa9a9543aa07eba8cdab36d787f1ff98a06db995c43cd5c63581ce050e0b9ba856634dabfaf8c7f271fbd026edd6ea1257b16013a526e0581a688cc6a335e7ee4c1b0633f0532d3d0824824195b6b249c70cf0e458609efc01a6575f084e6de53b
e = 0x10001
c = 0x6fadd5d7095bd6f45de69bb4e76080e0ea5f8c5a159de10663133e585b71ae580b99b3e0a8e047a9c51c8091a6b33b01c9ab95668794c3acfb084e939a04cb151757c3b2522da99e03f83e205c7c701066d69b120ca17fcf59061c078d9099e5f4bf6dd6dab206418527035f2c1096861c2896327977ac88c2728faa7504d879

p, q = fermat(n)

d = inverse(e, (p-1)*(q-1))

print(long_to_bytes(int(pow(c, d, n))))