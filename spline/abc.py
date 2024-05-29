e = 2.78

def sin(x):
    return (e ** (1j * x) - e ** (-1j * x)) / (2j)

from math import sin as sin_math


alpha = 0.44


print(f"sin({alpha}) = {sin(alpha)}")
print(f"math.sin({alpha}) = {sin_math(alpha)}")