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

e = 0x10001
n1 = 0xa8ed020c3dd125d503bf124052d643ba1405f2c349244122140e79e7d2244304a1590762c61ac83900c2aced76007b2e3f320464fd51fcfad167ebdc87e69329230869e0a3e153b44ed3b04bfe94174bc8b5ee1a3fa8036b6b9e834666aa07229a431b477e589d94f9a4cfed25b195215b0c694b86e874413b8a00cb064809c8e3677632cde9b43b87a0b812c2024b0c821b5c10764fd4de2d18af55d897d94aeded80b71e36fd73014f75641a8c5b38b36faa020e7cf1327a707bb7d42503bcc28768ef184d66b9ba16efd019b68268885a2da302cd326e78b1d473bcf7cd62442ccd25dc85d23aeb5408922b6b00f13584bea394f1bca4cc431f3c29c5d98ec1683453cc0c526abe4aec08781c7a53f50f2047b4995b9bea6a7a9f6b5425b29be6e867764efaa050799f716af78273041372cfe4f3c88a62329f6f1feff99
n2 = 0x16ea1bde86ea11cdec196a9173258efca235da66f8e3d5437e39e1b2e2574dd3f93d65104ca0225d6119519ae9ea9c035e0f85f02212c0992d0705723fa8b97ed6cff860c4d8fb65f0214a0047feca64e662dcbf025fff47590305e90e5d070d39871880828f5e960ab2ef330129ed5752c3b4debe827376a632b06487740fff4b622a88de23649e3e6993cd332b0284b84eb8765d58527209cc202c89d479421131a2f64ae517ee1e62e6c0f329c306569e427113ec6a8b9d96d73e95580d3a33f6add681f9a9156f0681eb1804183dfa8cebbe921d2fb1d43b256f727d46c5859cc5229f7e555ad25397e5cd14620ebbaefa0a0a520bada3ef8b115481734242af6befbd9b069d4a03281094c0f4aca4e6fdcbe2558b104fc2b383e1c70f0e5a07d1a623f9fc2309ca1d09b69aa1869e280fbc50de2adbada7ea545743b12b
c1 = 0x690e49037fee7649033ffaaa71e4730d2d7143fab97beb22e2afdf6eca449cad3f95b60295f592e7e84833e08b3468d61a34c1d1123f4c683c79d68bbe27dd0af203fc50ef7ebe98b1bc1221918470f058a8fb7645eacb569931835bd7f80494dbb67fbaa592ec19d9b4930c787a2ce1267f8088229b5031e710d6cd5720756923ccb64444939a0f09a51c87488650d4d02551fd4ed7a2fd248825ec34c5df8b6077a6d0d75c5832f9140420c92d3d00cf51e3b0665f5a6d031cb369ddebbb5ce77f2176cd12bb0add5aeda6ae88c4ceade0c1fd0ff3960d3ee36a0c6455ae3027f33e660663d0e2298654e19e8c8a06b4de991fac3b4c1673825b3d9f8f5c675f920a7d137f85ba723bf741321904e0c3c601f5c18d02e1e5b7b118e62e91a7926a9b1eda3cc53e2a6cbc95553e1990ec3f6cceddf283410d6e6849a26f89b
c2 = 0x5c4c7dce82753a68dcdbcdce9af52c9b7af2f561c08b8e23b27c6145d4c3df29d498303bee1bd29829a2e0ae9faaf243b387c39d69daccba07dace7bb420115ffaa69f89a3ea4e1ef0e08eb19043e012a090b79e51d6ae8446ca76e88abe5adbdbe25a731d7ee9aa333a84447edafbc360b505ff293c751571c6bf29dee99fdc443b756f182eb588b4a03de3d35dc4f23736d7239cfbd0ca13fa7b234bc4064a2053ab0045f4833250c8c9de91798502b09d4312ee52f3dc5229dfcb73b42f7c3440932839e6e790bb0db1788fbd7c60365121bbe3858ecedd3d48261d081c380e7ddf6ca570c13cc89c0af2011b4978b22d5456d1122dabd7b2068ab30e301a674809732daede77a27ae13e1bc4779e15d51210f6c10be159907ec1a59bfaf8db6cf290a348f734fd88e3c2b7df6bda84665b810cfe55bc3645d8d118c9172

cf = continued_fraction(n1/n2)
fracs = cf.convergents()
for xx in tqdm(fracs):
    q = xx.numerator()
    r = xx.denominator()
    try:
        if n1 % q == 0 and q != 1:
            print('find')
            p = n1 // q
            pp = n2 // r
            d = inverse(e, (p-1)*(q-1))
            print(long_to_bytes(int(pow(c1, d, n1))))
            
            d = inverse(e, (pp-1)*(r-1))
            print(long_to_bytes(int(pow(c2, d, n2))))

            exit(0)
    except:
        pass