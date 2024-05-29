import numpy as np
from specs import a, b

def create_nodes(n, chebyshev=False):
    if chebyshev:
        return [(b - a) / 2 * np.cos((2 * k - 1) * np.pi / (2 * n)) + (a + b) / 2 for k in range(1, n+1)]
    else:
        return np.linspace(a, b, n)
