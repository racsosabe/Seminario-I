from mpmath import *
import threading
import logging
import time
import sys
import multiprocessing
import concurrent.futures

def norm2(a):
    return a.real**2 + a.imag**2

def getMin(a, b): # Get min(|a + t(b - a)|) with t in [0, 1]
    prod = a * b.conjugate()
    if prod.real >= norm2(a): return abs(a)
    if prod.real >= norm2(b): return abs(b)
    return abs(prod.imag / abs(b-a))

def f(x, n): # f_{n}(x)
    global L
    ans = 1
    for i in range(2,n+1):
        ans += exp(-x * log(i))
    return ans

def df(x, n): # f_{n}'(x)
    global L
    ans = 0
    for i in range(2,n+1):
        ans += - exp(-x * log(i)) * log(i)
    return ans

def M(s1, s2, n): # Upper bound
    global L
    ans = 0
    for i in range(2,n+1):
        ans += exp((-min(s1.real,s2.real)) * log(i)) * log(i) * log(i)
    return ans

def isGreater(s1, s2, n): # Ying et Katz condition for Argument Variation
    L = getMin(f(s1, n),f(s2, n))
    R = M(s1, s2, n) * norm2(s1 - s2) / 8
    EPS = mpf('10') ** (-20)
    return L > R + EPS

def enrect(z1, z2, x): # Check if x belongs to rectangle
    return min(z1.real, z2.real) <= x.real <= max(z1.real, z2.real) and min(z1.imag, z2.imag) <= x.imag <= max(z1.imag, z2.imag)

def varArg(s1, s2, n): # Argument Variation function
    t0 = mpf(0)
    t = mpf(1)
    L = s1
    d = s2 - s1
    ans = mpf(0)
    EPS1 = mpf('10') ** (-25)
    EPS2 = mpf('10') ** (-50)
    EPS3 = mpf('10') ** (-35)
    while t0 + EPS1 < 1:
        t = mpf(abs((s2 - L) / (s2 - s1))) # Start with maximum t possible
        R = L + t * d
        while (not isGreater(L,R,n)) and t > EPS2: # If it doesn't work, divide it by 2
            t /= 2
            R = L + t * d
        if t <= EPS2: # Might be a zero on the segment
            if s1.real == s2.real: # Detour by Re
                d = EPS3 if s1.imag < s2.imag else -EPS3
                return varArg(s1 + mpc(d,0), s2 + mpc(d,0), n)
            else: # Detour by Im
                d = EPS3 if s1.real > s2.real else -EPS3
                return varArg(s1 + mpc(0,d), s2 + mpc(0,d), n)
        ans += arg(f(R, n) / f(L, n)) # No zeros -> Add current variation
        t0 += t
        L += d * t
    return ans

def getRectangle(LD, RU, n):
    LU = mpc(LD.real,RU.imag)
    RD = mpc(RU.real,LD.imag)
    return round((varArg(LD,RD,n) + varArg(RD,RU,n) + varArg(RU,LU,n) + varArg(LU,LD,n)) / (2 * acos(-1)))

def newton(LD, RU, n): # Newton's method
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
        self.stack = [LD,RU,varArg(LD,RD,_n),varArg(RD,RU,_n),varArg(RU,LU,_n),varArg(LU,LD,_n)]

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

def compute(a, n):
    LD = a[0] # Left-down limit
    RU = a[1] # Right-up limit
    V0 = a[2] # Argument variation from bottom side
    V1 = a[3] # Argument variation from right side
    V2 = a[4] # Argument variation from top side
    V3 = a[5] # Argument variation from left side
    zeros = round(sum(a[2:]) / (2 * acos(mpf('-1'))))
    if zeros == 0: return [0, mpc(3, 0)] # No zeros flag
    if zeros == 1:
        O = (LD + RU) / 2
        res = newton(LD, RU, n)
        if enrect(LD, RU, res):
            return [1, res] # Found one zero
    if abs(LD.real - RU.real) > abs(LD.imag - RU.imag): # Largest dimension is Re(s), so split vertically
        M1 = mpc((LD.real + RU.real) / 2, LD.imag) # Bottom side mid point
        M2 = mpc((LD.real + RU.real) / 2, RU.imag) # Top side mid point
        FD = varArg(LD,M1,n) # Just compute argument variation from LD to M1 (M1 to right-down corner can be computed by substracting it from V0)
        FU = varArg(RU,M2,n) # Just compute argument variation from RU to M2 (M2 to left-up corner can be computed by substracting it from V2)
        FM = varArg(M1,M2,n) # Compute argument variation between mid points
        v1 = FD + FM + V2 - FU + V3 # Argument variation from left rectangle
        v2 = V0 - FD + V1 + FU - FM # Argument variation from right rectangle
        if round((v1 + v2) / (2 * acos(-1))) != zeros: # Check if it's consistent
            print("Fatal Error. Sum of partitions isn't equal to total")
        return [2, LD, M2, FD, FM, V2 - FU, V3, M1, RU, V0 - FD, V1, FU, -FM] # Return the rectangles.
    else: # Largest dimension is Im(s), so split horizontally
        M1 = mpc(RU.real, (LD.imag + RU.imag) / 2) # Right side mid point
        M2 = mpc(LD.real, (LD.imag + RU.imag) / 2) # Left side mid point
        FD = varArg(mpc(RU.real,LD.imag),M1,n) # Just compute argument variation from right-down corner to M1 (M1 to right-up corner can be computed by substracting it from V1)
        FU = varArg(mpc(LD.real,RU.imag),M2,n) # Just compute argument variation from left-up corner to M2 (M2 to right-up corner can be computed by substracting it from V3)
        FM = varArg(M1,M2,n) # Compute argument variation between mid points
        v1 = V0 + FD + FM + (V3 - FU) # Argument variation from bottom rectangle
        v2 = -FM + (V1 - FD) + V2 + FU # Argument variation from top rectangle
        if round((v1 + v2) / (2 * acos(-1))) != zeros:
            print("Fatal Error. Sum of partitions isn't equal to total")
        return [2, LD, M1, V0, FD, FM, V3 - FU, M2, RU, -FM, V1 - FD, V2, FU]

def initialize(S, n, nodes):
    # Fill the stack until all nodes can grab a rectangle
    while S.size() < 6 * nodes:
        cur = S.pop()
        if len(cur) == 0: break # Break when no zeros left
        val = compute(cur, n)
        if val[0] == 1: # Found a zero
            with open(outputname, 'a') as out:
                out.write(str(val[1].real))
                out.write(' ')
                out.write(str(val[1].imag))
                out.write(' ')
                out.write(str(abs(f(val[1], n))))
                out.write('\n')
        elif val[0] == 2: # Have two rectangles
            S.push(val[1:])
        # This ignores val[0] == 3 since there are no zeros


def solve(S, id, n):
    global outputname
    failed = 0
    while failed < 12: # Stop after 12 fails
        cur = S.pop()
        if len(cur) == 6:
            val = compute(cur, n)
            if val[0] == 1: # Found a zero
                with open(outputname, 'a') as out:
                    out.write(str(val[1].real))
                    out.write(' ')
                    out.write(str(val[1].imag))
                    out.write(' ')
                    out.write(str(abs(f(val[1], n))))
                    out.write('\n')
            elif val[0] == 2: # Have two rectangles
                S.push(val[1:])
            failed = 0 # Restart failed count since there was progress
        else: # Add one to failed count 
            failed += 1

rate = 1000
limit = 100000
mp.dps = 100
n = int(sys.argv[1])
outputname = 'zeros{:d}.txt'.format(n)
print(outputname)
threads = multiprocessing.cpu_count()

with open('checkpoint', 'r') as Stdin:
    Limag = int(Stdin.read())
    Rimag = rate + Limag

if Limag == 0:
    with open(outputname, 'w') as out:
        out.write('')

L = [0 if i == 0 else log(i) for i in range(n + 1)]

while Limag < limit:
    LD = mpc(1 - n, str(Limag))
    RU = mpc('1.74', str(Rimag))
    print("Starting range ", Limag, Rimag)
    print(LD, RU)
    start = time.time()
    S = MutexStack(LD, RU, n)
    initialize(S, n, threads)
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        for index in range(threads):
            executor.submit(solve, S, index, n)
    end = time.time()
    print("Ended processing range [{:d}, {:d}] in {:f} seconds".format(Limag, Rimag, end - start))
    Limag += rate
    Rimag += rate
    with open('checkpoint', 'w') as check:
        check.write(str(Limag))

print("Ended processing the full range [0, {:d}] :D".format(limit))