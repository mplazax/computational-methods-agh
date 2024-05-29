import numpy as np
from metoda_newtona import metoda_newtona
from metoda_siecznych import metoda_siecznych
import specs

def znajdz_miejsca_zerowe(f=specs.f, df=specs.df, a=specs.a, b=specs.b, n=specs.n, metoda='newton', max_iter=specs.max_iter, tol=specs.tol):
    przedzialy = np.linspace(a, b, n+1)
    miejsca_zerowe = []
    for i in range(n):
        x0, x1 = przedzialy[i], przedzialy[i+1]
        if f(x0) * f(x1) <= 0:
            if metoda == 'newton':
                root = metoda_newtona(f, df, (x0 + x1) / 2)
            else:
                root = metoda_siecznych(f, x0, x1)
            if root is not None and a <= root <= b:
                miejsca_zerowe.append(root)
    # Usunięcie duplikatów i sortowanie wyników
    miejsca_zerowe = sorted(set(miejsca_zerowe))
    return miejsca_zerowe

if __name__ == "__main__":
    # Przykładowe funkcje
    def f(x):
        return x**3 - x**2 - x + 1

    def df(x):
        return 3*x**2 - 2*x - 1

    # Wyszukiwanie wszystkich miejsc zerowych na przedziale [-2, 2]
    roots = znajdz_miejsca_zerowe(f, df, -2, 2, 100)
    print("Miejsca zerowe funkcji:", roots)


