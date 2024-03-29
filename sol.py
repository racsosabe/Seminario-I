import threading
import logging
import time
import multiprocessing
import concurrent.futures
from cmath import phase, exp
from math import log, acos, fsum
from math import exp as expR

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
        ans += exp(-log(i) * x)
    return ans

def df(x, n):
    ans = 0
    for i in range(2,n+1):
        ans += - log(i) * exp(- log(i) * x)
    return ans

def M(s1, s2, n):
    ans = 0
    for i in range(2,n+1):
        ans += log(i) * log(i) * expR(log(i) * (-min(s1.real,s2.real)))
    return ans

def isGreater(s1, s2, localM, n):
    L = getMin(f(s1, n),f(s2, n))
    R = localM * norm(s1 - s2) / 8
    return L > R + 1e-20

def enrect(z1, z2, x):
    return min(z1.real, z2.real) <= x.real <= max(z1.real, z2.real) and min(z1.imag, z2.imag) <= x.imag <= max(z1.imag, z2.imag)

def varArg(s1, s2, n):
    localM = M(s1,s2,n)
    t0 = 0
    t = 1
    L = s1
    d = s2 - s1
    ans = 0
    while t0 + 1e-20 < 1:
        R = L + t * d
        while (not isGreater(L,R,localM,n)) and t > 1e-15:
            t /= 2
            R = L + t * d
        if t <= 1e-15:
            print("Failed")
            if s1.real == s2.real:
                d = 1e-10 if s1.imag < s2.imag else -1e-10
                return varArg(s1 + complex(d,0), s2 + complex(d,0), n)
            else:
                d = 1e-10 if s1.real > s2.real else -1e-10
                return varArg(s1 + complex(0,d), s2 + complex(0,d), n)
        ans += phase(f(R, n) / f(L, n))
        t0 += t
        L += d * t
        t = min(2 * t, 1 + 1e-20 - t0)
    return ans

def getRectangle(LD, RU, n):
    LU = complex(LD.real,RU.imag)
    RD = complex(RU.real,LD.imag)
    return round((varArg(LD,RD,n) + varArg(RD,RU,n) + varArg(RU,LU,n) + varArg(LU,LD,n)) / (2 * acos(-1)))

def norm1(x):
    return max(abs(x.real), abs(x.imag))

def newton(LD, RU, n):
    it = 10000
    x0 = (LD + RU) * 0.5
    xi = x0 + complex(1,0)
    tol = 10**(-20)
    while norm1(x0 - xi) > tol and enrect(LD, RU, x0) and it > 0:
        xi = x0
        x0 = xi - f(xi, n) / df(xi, n)
        it -= 1 
    if min(norm1(x0 - xi), norm1(f(x0, n))) <= 10**(-12) and enrect(LD,RU,x0):
        return x0
    return complex(3,0)

class MutexStack:
    def __init__(self, LD, RU, _n):
        self.stack = []
        self.n = _n
        self._lock = threading.Lock()
        LU = complex(LD.real,RU.imag)
        RD = complex(RU.real,LD.imag)
        self.stack = [LD,RU,varArg(LD,RD,_n),varArg(RD,RU,_n),varArg(RU,LU,_n),varArg(LU,LD,_n)]

    def work(self, name):
        n = self.n
        with self._lock:
            if self.stack:
                LD = self.stack[-6]
                RU = self.stack[-5]
                S = fsum(self.stack[-4:])
                zeros = round(S / (2 * acos(-1)))
                V0 = self.stack[-4]
                V1 = self.stack[-3]
                V2 = self.stack[-2]
                V3 = self.stack[-1]
                for _ in range(6): self.stack.pop()
                if zeros == 0: return [2, complex(3,0)]
                if zeros == 1:
                    O = (LD + RU) / 2
                    res = newton(LD, RU, n)
                    if enrect(LD, RU, res):
                        return [1, res]
                if abs(LD.real - RU.real) > abs(LD.imag - RU.imag):
                    M1 = complex((LD.real + RU.real) / 2, LD.imag)
                    M2 = complex((LD.real + RU.real) / 2, RU.imag)
                    FD = varArg(LD,M1,n)
                    FU = varArg(RU,M2,n)
                    FM = varArg(M1,M2,n)
                    for x in [LD,M2,FD,FM,V2 - FU,V3]:
                        self.stack.append(x)
                    for x in [M1,RU,V0 - FD,V1,FU,-FM]:
                        self.stack.append(x)
                else:
                    M1 = complex(RU.real, (LD.imag + RU.imag) / 2)
                    M2 = complex(LD.real, (LD.imag + RU.imag) / 2)
                    FD = varArg(complex(RU.real,LD.imag),M1,n)
                    FU = varArg(complex(LD.real,RU.imag),M2,n)
                    FM = varArg(M1,M2,n)
                    for x in [LD,M1,V0,FD,FM,V3-FU]:
                        self.stack.append(x)
                    for x in [M2,RU,-FM,V1-FD,V2,FU]:
                        self.stack.append(x)
                return [2, complex(3,0)]
            else:
                return [0, complex(3,0)]


class MutexAnswer:
    def __init__(self):
        self.v = []
        self.len = 0
        self._lock = threading.Lock()

    def push(self, x):
        with self._lock:
            self.v.append(x)
            self.len += 1

    def pop(self):
        with self._lock:
            if self.v:
                res = self.v.pop();
                self.len -= 1
                return res
            else:
                return complex(3,0)

    def size(self):
        return self.len

def solve(S, R, id):
    val = S.work(id)
    while val[0] >= 1:
        if val[0] == 1:
            print("Found answer from {0:d}: ".format(id),val[1])
            R.push(val[1])
        val = S.work(id)

threads = multiprocessing.cpu_count()
S = MutexStack(complex(-1, 0), complex(1, 50000), 3)
R = MutexAnswer()
with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
    for index in range(threads):
        executor.submit(solve, S, R, index)
print(getRectangle(complex(-1,0),complex(1,50000),3))
print("Answers:")
print(R.size())
res = R.pop()
while res.real <= 1:
    print(res.real,res.imag)
    res = R.pop()
