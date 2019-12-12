from mpmath import *
import threading
import logging
import time
import sys
import multiprocessing
import concurrent.futures

def norm2(a):
    return a.real**2 + a.imag**2

def getMin(a, b):
    prod = a * b.conjugate()
    if prod.real >= norm2(a): return abs(a)
    if prod.real >= norm2(b): return abs(b)
    return abs(prod.imag / abs(b-a))

def f(x, n):
    ans = 1
    for i in range(2,n+1):
        ans += exp(-x * log(i))
    return ans

def df(x, n):
    ans = 0
    for i in range(2,n+1):
        ans += - exp(-x * log(i)) * log(i)
    return ans

def M(s1, s2, n):
    ans = 0
    for i in range(2,n+1):
        ans += exp((-min(s1.real,s2.real)) * log(i)) * log(i) * log(i)
    return ans

def isGreater(s1, s2, localM, n):
    L = getMin(f(s1, n),f(s2, n))
    R = localM * norm2(s1 - s2) / 8
    EPS = mpf('10') ** (-20)
    return L > R + EPS

def enrect(z1, z2, x):
    return min(z1.real, z2.real) <= x.real <= max(z1.real, z2.real) and min(z1.imag, z2.imag) <= x.imag <= max(z1.imag, z2.imag)

def varArg(s1, s2, n):
    localM = M(s1,s2,n)
    t0 = mpf(0)
    t = mpf(1)
    L = s1
    d = s2 - s1
    ans = mpf(0)
    EPS1 = mpf('10') ** (-25)
    EPS2 = mpf('10') ** (-80)
    EPS3 = mpf('10') ** (-40)
    while t0 + EPS1 < 1:
        R = L + t * d
        while (not isGreater(L,R,localM,n)) and t > EPS2:
            t /= 2
            R = L + t * d
        if t <= EPS2:
            if s1.real == s2.real:
                d = EPS3 if s1.imag < s2.imag else -EPS3
                return varArg(s1 + mpc(d,0), s2 + mpc(d,0), n)
            else:
                d = EPS3 if s1.real > s2.real else -EPS3
                return varArg(s1 + mpc(0,d), s2 + mpc(0,d), n)
        ans += arg(f(R, n) / f(L, n))
        t0 += t
        L += d * t
        t = min(2 * t, 1 + EPS1 - t0)
    return ans

def getRectangle(LD, RU, n):
    LU = mpc(LD.real,RU.imag)
    RD = mpc(RU.real,LD.imag)
    return round((varArg(LD,RD,n) + varArg(RD,RU,n) + varArg(RU,LU,n) + varArg(LU,LD,n)) / (2 * acos(-1)))

def newton(LD, RU, n):
    it = 10000
    x0 = (LD + RU) * mpf('0.5')
    xi = x0 + mpc(1,0)
    tol = mpf('10') ** (-30)
    while abs(x0 - xi) > tol and enrect(LD, RU, x0) and it > 0:
        xi = x0
        x0 = xi - f(xi, n) / df(xi, n)
        it -= 1 
    if abs(x0 - xi) <= tol and enrect(LD,RU,x0):
        return x0
    return mpc(3,0)

class MutexStack:
    def __init__(self, LD, RU, _n):
        self.stack = []
        self.n = _n
        self._lock = threading.Lock()
        LU = mpc(LD.real,RU.imag)
        RD = mpc(RU.real,LD.imag)
        print(LD, RU, _n)
        self.stack = [LD,RU,varArg(LD,RD,_n),varArg(RD,RU,_n),varArg(RU,LU,_n),varArg(LU,LD,_n)]
        print(self.stack)

    def pop(self):
        n = self.n
        with self._lock:
            if self.stack:
                res = self.stack[-6:]
                for _ in range(6): self.stack.pop()
                return res
            else:
                return []

    def push(self, v):
        with self._lock:
            for x in v:
                self.stack.append(x)

    def size(self):
        return len(self.stack)


class MutexAnswer:
    def __init__(self, threads):
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
                return mpc('3','0')

    def size(self):
        return self.len

def compute(a, n):
    LD = a[0]
    RU = a[1]
    V0 = a[2]
    V1 = a[3]
    V2 = a[4]
    V3 = a[5]
    zeros = round(sum(a[2:]) / (2 * acos(mpf('-1'))))
    if zeros == 0: return [0, mpc(3, 0)]
    if zeros == 1:
        O = (LD + RU) / 2
        res = newton(LD, RU, n)
        if enrect(LD, RU, res):
            return [1, res]
    if abs(LD.real - RU.real) > abs(LD.imag - RU.imag):
        M1 = mpc((LD.real + RU.real) / 2, LD.imag)
        M2 = mpc((LD.real + RU.real) / 2, RU.imag)
        FD = varArg(LD,M1,n)
        FU = varArg(RU,M2,n)
        FM = varArg(M1,M2,n)
        v1 = FD + FM + V2 - FU + V3
        v2 = V0 - FD + V1 + FU - FM
        if round((v1 + v2) / (2 * acos(-1))) != zeros:
            print("Fatal Error. Sum of partitions isn't equal to total")
        return [2, LD, M2, FD, FM, V2 - FU, V3, M1, RU, V0 - FD, V1, FU, -FM]
    else:
        M1 = mpc(RU.real, (LD.imag + RU.imag) / 2)
        M2 = mpc(LD.real, (LD.imag + RU.imag) / 2)
        FD = varArg(mpc(RU.real,LD.imag),M1,n)
        FU = varArg(mpc(LD.real,RU.imag),M2,n)
        FM = varArg(M1,M2,n)
        v1 = V0 + FD + FM + (V3 - FU)
        v2 = -FM + (V1 - FD) + V2 + FU
        if round((v1 + v2) / (2 * acos(-1))) != zeros:
            print("Fatal Error. Sum of partitions isn't equal to total")
        return [2, LD, M1, V0, FD, FM, V3 - FU, M2, RU, -FM, V1 - FD, V2, FU]

def initialize(S, R, n, nodes):
    while S.size() < 6 * nodes:
        cur = S.pop()
        val = compute(cur, n)
        if val[0] == 1:
            R.push(val[1])
        elif val[0] == 2:
            S.push(val[1:])


def solve(S, R, id, n):
    failed = 0
    while failed < 12:
        cur = S.pop();
        if len(cur) == 6:
            val = compute(cur, n)
            if val[0] == 1:
                R.push(val[1])
            elif val[0] == 2:
                S.push(val[1:])
            failed = 0
        else:
            failed += 1

mp.dps = 100
n = int(sys.argv[1])
threads = multiprocessing.cpu_count()
LD = mpc(1 - n,'0')
RU = mpc('1.74','10000')
start = time.time()
S = MutexStack(LD,RU,n)
R = MutexAnswer(threads)
initialize(S, R, n, threads)
with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
    for index in range(threads):
        executor.submit(solve, S, R, index, n)
end = time.time()
res = R.pop()
total = R.size()
while res.real <= 1:
    print(res.real, res.imag)
    res = R.pop()
