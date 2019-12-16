from mpmath import *
import sys

mp.dps = 100
b = mpf(2) * acos(mpf(-1)) / log(mpf(5) / mpf(3))

for line in sys.stdin:
    x, y = [x for x in line.split()]
    Re = mpf(x)
    Im = mpf(y)
    print(Re, fmod(Im, b))

