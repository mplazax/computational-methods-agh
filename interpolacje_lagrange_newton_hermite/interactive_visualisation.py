import numpy as np
from visualize_interpolation import visualize_interpolation

if __name__ == "__main__":
    f = lambda x: np.exp(-np.sin(x)) + np.cos(x)
    df = lambda x: -np.cos(x) * np.exp(-np.sin(x)) - np.sin(x)
    a, b = -3 * np.pi, 5 * np.pi

    while True:
        metoda = str(input("Podaj metode (l/n/h): "))
        metoda = 'lagrange' if metoda == 'l' else 'newton' if metoda == 'n' else 'hermite'
        n = int(input("Podaj n (int>0): "))
        chebyshev = str(input("Czy chcesz uzyc wezlow Czebyszewa? (t/n): "))
        chebyshev = True if chebyshev == 't' else False
        visualize_interpolation(f, a, b, n, metoda, chebyshev=chebyshev, df=df)
