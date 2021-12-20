# AGCD: `xi = p*qi+ri` for 1<=i<=t, where ri is small.
#        Given some xi, solve for common divisor p.

from Crypto.Util.number import *


def SDA_solver(xs, rho):
    """
    Basic idea: xi / x0 = qi / q0
    Construct lattice to attack.
    B = [2^(rho+1)  x1   x2  ....  xt]
        [           -x0              ]
        [                -x0         ]
        [                   .        ]
        [                     .      ]
        [                       .    ]
        [                         -x0]
    v = (q0, q1, ..., qt)B
      = (q0*2^(rho+1), q0r1-q1r0, ..., q0rt-qtr0)
    B.LLL() to solve for q0 so we can solve for p
    """

    # 1. Construct lattice
    t = len(xs) - 1
    B = Matrix(ZZ, t+1, t+1)
    for i in range(t+1):
        B[i, i] = -xs[0]
        if i == 0:
            B[0, i] = 2^(rho+1)
        else:
            B[0, i] = xs[i]
    # 2. LLL and find p
    v = B.LLL()[0]
    q0 = v[0] // 2^(rho+1)
    p = xs[0] // q0
    return p

def MP_solver(xs, rho):
    """
    Multiivariate polynomial approach (MP)
    similar to Multiivariate-coppersmith construction
    Using polynomials to construct lattice and solve for p
    """

    X = 1<<rho
    m = len(xs)
    PR = PolynomialRing(ZZ, names=[str('x%d' % i) for i in range(1, m+1)])
    h = 3
    u = 1
    variables = PR.gens()
    gg = []
    monomials = [variables[0] ** 0]
    for i in range(m):
        gg.append(xs[i] - variables[i])
        monomials.append(variables[i])
    print(len(monomials), len(gg))
    print('monomials:', monomials)

    B = Matrix(ZZ, len(gg), len(monomials))
    for ii in range(len(gg)):
        for jj in range(len(monomials)):
            if monomials[jj] in gg[ii].monomials():
                B[ii, jj] = gg[ii].monomial_coefficient(monomials[jj]) * monomials[jj]([X] * m)
    B = B.LLL()

    new_poly = []
    for i in range(len(gg)):
        tmp_poly = 0
        for j in range(len(monomials)):
            tmp_poly += monomials[j](variables) * B[i, j] / monomials[j]([X] * m)
        new_poly.append(tmp_poly)
    
    if len(new_poly) > 0:
        Ideal = ideal(new_poly[:m-1])
        GB = Ideal.groebner_basis()
        function_variables = var([str('y%d' % i) for i in range(1, 1 + m)])
        res = solve([pol(function_variables) for pol in GB], function_variables)
        print('got %d basis' % len(GB))
        print('solved result:')
        print(res)
        for tmp_res in res:
            PRRR.<x, y> = PolynomialRing(QQ)
            q = abs(PRRR(res[0][0](x, y)).coefficients()[0].denominator())
            p = xs[-1] // q
            # print(p)
            return p

def generate_testcase(rbits, pbits):
    rbit = int(rbits)
    pbit = int(pbits)
    p = getPrime(pbit)
    ns = [p*getPrime(pbit) + getPrime(rbit) for i in range(5)]
    return p, ns

if __name__ == '__main__':
    # p, ns = generate_testcase(300, 512)
    # recover_p_SDA = SDA_solver(ns, 300)
    # recover_p_MP = MP_solver(xs = ns, rho = 300)
    # if p == recover_p_SDA:
    #     print('Pass SDA method!!')
    # if p == recover_p_SDA:
    #     print('Pass SDA method!!')
    e = 0x10001
    n1 = 0xa8ed020c3dd125d503bf124052d643ba1405f2c349244122140e79e7d2244304a1590762c61ac83900c2aced76007b2e3f320464fd51fcfad167ebdc87e69329230869e0a3e153b44ed3b04bfe94174bc8b5ee1a3fa8036b6b9e834666aa07229a431b477e589d94f9a4cfed25b195215b0c694b86e874413b8a00cb064809c8e3677632cde9b43b87a0b812c2024b0c821b5c10764fd4de2d18af55d897d94aeded80b71e36fd73014f75641a8c5b38b36faa020e7cf1327a707bb7d42503bcc28768ef184d66b9ba16efd019b68268885a2da302cd326e78b1d473bcf7cd62442ccd25dc85d23aeb5408922b6b00f13584bea394f1bca4cc431f3c29c5d98ec1683453cc0c526abe4aec08781c7a53f50f2047b4995b9bea6a7a9f6b5425b29be6e867764efaa050799f716af78273041372cfe4f3c88a62329f6f1feff99
    n2 = 0x16ea1bde86ea11cdec196a9173258efca235da66f8e3d5437e39e1b2e2574dd3f93d65104ca0225d6119519ae9ea9c035e0f85f02212c0992d0705723fa8b97ed6cff860c4d8fb65f0214a0047feca64e662dcbf025fff47590305e90e5d070d39871880828f5e960ab2ef330129ed5752c3b4debe827376a632b06487740fff4b622a88de23649e3e6993cd332b0284b84eb8765d58527209cc202c89d479421131a2f64ae517ee1e62e6c0f329c306569e427113ec6a8b9d96d73e95580d3a33f6add681f9a9156f0681eb1804183dfa8cebbe921d2fb1d43b256f727d46c5859cc5229f7e555ad25397e5cd14620ebbaefa0a0a520bada3ef8b115481734242af6befbd9b069d4a03281094c0f4aca4e6fdcbe2558b104fc2b383e1c70f0e5a07d1a623f9fc2309ca1d09b69aa1869e280fbc50de2adbada7ea545743b12b
    c1 = 0x690e49037fee7649033ffaaa71e4730d2d7143fab97beb22e2afdf6eca449cad3f95b60295f592e7e84833e08b3468d61a34c1d1123f4c683c79d68bbe27dd0af203fc50ef7ebe98b1bc1221918470f058a8fb7645eacb569931835bd7f80494dbb67fbaa592ec19d9b4930c787a2ce1267f8088229b5031e710d6cd5720756923ccb64444939a0f09a51c87488650d4d02551fd4ed7a2fd248825ec34c5df8b6077a6d0d75c5832f9140420c92d3d00cf51e3b0665f5a6d031cb369ddebbb5ce77f2176cd12bb0add5aeda6ae88c4ceade0c1fd0ff3960d3ee36a0c6455ae3027f33e660663d0e2298654e19e8c8a06b4de991fac3b4c1673825b3d9f8f5c675f920a7d137f85ba723bf741321904e0c3c601f5c18d02e1e5b7b118e62e91a7926a9b1eda3cc53e2a6cbc95553e1990ec3f6cceddf283410d6e6849a26f89b
    c2 = 0x5c4c7dce82753a68dcdbcdce9af52c9b7af2f561c08b8e23b27c6145d4c3df29d498303bee1bd29829a2e0ae9faaf243b387c39d69daccba07dace7bb420115ffaa69f89a3ea4e1ef0e08eb19043e012a090b79e51d6ae8446ca76e88abe5adbdbe25a731d7ee9aa333a84447edafbc360b505ff293c751571c6bf29dee99fdc443b756f182eb588b4a03de3d35dc4f23736d7239cfbd0ca13fa7b234bc4064a2053ab0045f4833250c8c9de91798502b09d4312ee52f3dc5229dfcb73b42f7c3440932839e6e790bb0db1788fbd7c60365121bbe3858ecedd3d48261d081c380e7ddf6ca570c13cc89c0af2011b4978b22d5456d1122dabd7b2068ab30e301a674809732daede77a27ae13e1bc4779e15d51210f6c10be159907ec1a59bfaf8db6cf290a348f734fd88e3c2b7df6bda84665b810cfe55bc3645d8d118c9172

    p = SDA_solver([n1, n2], 512)
    assert n1 % p == 0

    pp = next_prime(p)
    q = n1 // p

    r = n2 // pp

    d = inverse(e, (p-1)*(q-1))
    m1 = int(pow(c1, d, n1))

    d = inverse(e, (pp-1)*(r-1))
    m2 = int(pow(c2, d, n2))
    print(long_to_bytes(m1) + long_to_bytes(m2))

    