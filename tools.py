from cmath import phase
from math import log, acos

def norm(a):
    return a.real**2 + a.imag**2

def getMin(a, b):
    prod = a * b.conjugate()
    if prod.real >= norm(a): return abs(a)
    if prod.real >= norm(b): return abs(b)
    return abs(prod.imag / abs(b-a))

def f(x, n):
    ans = 1
    for i in range(2,n+1):
        ans += i ** (-x)
    return ans

def df(x, n):
    ans = 0
    for i in range(2,n+1):
        ans += -x * i ** (-x) * log(i)
    return ans

def M(s1, s2, n):
    ans = 0
    for i in range(2,n+1):
        ans += log(i) * log(i) * i ** (-min(s1.real,s2.real))
    return ans

def isGreater(s1, s2, localM, n):
    L = getMin(f(s1, n),f(s2, n))
    R = localM * norm(s1 - s2) / 8
    return L > R + 1e-20

def varArg(s1, s2, n):
    localM = M(s1,s2,n)
    t0 = 0
    t = 1
    L = s1
    d = s2 - s1
    ans = 0
    while t0 + 1e-20 < 1:
        R = L + t * d
        print(t)
        while not isGreater(L,R,localM,n) and t > 1e-15:
            t /= 2
            R = L + t * d
        if t <= 1e-15:
            print("Failed")
            return -1e30
        ans += phase(f(R, n) / f(L, n))
        t0 += t
        L += d * t
        t = min(2 * t, 1 + 1e-20 - t0)
    return ans

def getRectangle(LD, RU, n):
    LU = complex(LD.real,RU.imag)
    RD = complex(RU.real,LD.imag)
    return (varArg(LD,RD,n) + varArg(RD,RU,n) + varArg(RU,LU,n) + varArg(LU,LD,n)) / (2 * acos(-1))

def newton(x0, R, n):
    xi = x0
    while abs(x0 - xi) <= R and abs(f(xi,n)) > 1e-15:
        xi = xi - f(xi, n) / df(xi, n)
    if abs(x0 - xi) > R: return complex(10,0)
    return xi
