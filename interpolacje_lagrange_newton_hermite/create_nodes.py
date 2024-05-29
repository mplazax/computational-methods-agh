import numpy as np


def create_nodes(a, b, n, chebyshev=False):
    if chebyshev:
        return [(b - a) / 2 * np.cos((2 * k - 1) * np.pi / (2 * n)) + (a + b) / 2 for k in range(1, n+1)]
    else:
        return np.linspace(a, b, n)
