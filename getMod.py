from mpmath import *
import sys

mp.dps = 100
b = mpf(2) * acos(mpf(-1)) / log(mpf(11) / mpf(2))

for line in sys.stdin:
    x, y = [x for x in line.split()]
    Re = mpf(x)
    Im = mpf(y)
    print(Re, fmod(Im, b))

