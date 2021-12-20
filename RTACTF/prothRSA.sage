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

n = 0xa19028b5c0e77e19fc167374358aa346776e6c20c27499505be59c83ea02014e97af631ba0ccbab881313818fd323c15c82dad8793220ba6679ec4b38787e04d0c1fff0880e04423ea288e443660c63a1607532e47dbaad421723d0546c208447f701cd7e9ee1bb43774d132abbb2e91bf50b67be40ed854dbe6c3071ca3ae3307ac03abd76f74e506594106a22795d4b7938611301248a9957e1a637538a9169cf38daf5d60ffc05ae32ea7e638e16d790ffeebfff655a645c99a513616d3ce00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
e = 0x10001
s = 0x28640a2d7039df867f059cdd0d62a8d19ddb9b08309d265416f96720fa808053a5ebd8c6e8332eae204c4e063f4c8f05720b6b61e4c882e999e7b12ce1e1f812c11cfed72a5c33cfb8f3d34f650e4c19579cf34745f2588aa2fd08a8746257cb789f23ca232346fcf72468a2b160934911902de3f90620aba5874a2d79a33699
c = 0x4595c3c923bd191ba07456611f80e656a197ff528a031e2952adedda532b1fa2caef719c929132a3cdf06d0e55e6a00f7eb1f189a614b26759916ec42f83579a75ab5948186769a1a936b019466f918f29e32852675c464b7f0797c6fdc55efcd54fbe2083761b1df3dde0b9a9a35b96e3b216c54770b444b1f02525f0268c44483c6e84a781fe9111e6912130d69f462c519873043d44e4a3f1f938491feeb591b5831d0abe7399bc87244576decaf2925f287d3c2bb4061d560c919d820e364744f2322c7efd37d42563842bcf9b1d6b46218694dcd49758d311c6896e38cf2b55c7114d78cfdfaeba74720ecf30d9133034799b9735e26ec913cc9f26bb0a

P.<k1, k2> = PolynomialRing(Zmod(n))

p = (2*k1 + 1) * (1<<512) + 1
q = (2*k2 + 1) * (1<<512) + 1

I = [p*q-n, k1*k2-s]

I = Ideal(I)
B = I.groebner_basis()

# print(B)

P.<x> = PolynomialRing(Zmod(n))

f = x^2 + 20395454563409867656617278217828189042208130201653144079515462893108235291959386336828853848840919483776716211380827788569533851096098894111797056850473973215977506955765638320965700609202656074252440878871479437642421736402628440279746509595855847453146967895481905399926423917985459031590852011496691628662604259167416232330264857609102516325571135277439495982072415418716300644648555729595994325721230901044096167777386912138442924390203810728808854693589734150101148074301478029960268708923366959989928005874227107020841823527632715276090921124020041240950439333518284946100275969081529147302101259207747197769243*x + 28363370488384189140808943644607274199262748227466385949335542116022077684563044493960859825144671982446528018743520526298520183802609154816402679682040187831399735385622546475759074209529903290929679123586974003866637756048478255324974862941367199328160197766840085731250929590153259657118084274765647001241
A = 20395454563409867656617278217828189042208130201653144079515462893108235291959386336828853848840919483776716211380827788569533851096098894111797056850473973215977506955765638320965700609202656074252440878871479437642421736402628440279746509595855847453146967895481905399926423917985459031590852011496691628662604259167416232330264857609102516325571135277439495982072415418716300644648555729595994325721230901044096167777386912138442924390203810728808854693589734150101148074301478029960268708923366959989928005874227107020841823527632715276090921124020041240950439333518284946100275969081529147302101259207747197769243

roots = f.small_roots(X=2**512, beta=0.4, epsilon=5/32)
for k2 in roots:
    q = int((2*k2 + 1) * (1<<512) + 1)
    d = int(inverse(e, q-1))
    print(long_to_bytes(int(pow(c, d, q))))
