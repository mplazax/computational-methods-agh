import specs


def metoda_siecznych(f, x0, x1, tol=specs.tol, max_iter=specs.max_iter):
    for i in range(max_iter):
        if abs(f(x1)) < tol:
            return x1
        # Wzór na kolejny przybliżony pierwiastek
        try:
            x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        except ZeroDivisionError:
            print("Dzielenie przez zero. Zmiana punktów startowych może być konieczna.")
            return None
        x0, x1 = x1, x2
    print("Przekroczono maksymalną liczbę iteracji. Rozwiązanie może nie być dokładne.")
    return x1

if __name__ == "__main__":
    # Przykładowa funkcja: f(x) = x^2 - 4
    def f(x):
        return x**2 - 4

    # Test metody siecznych
    root = metoda_siecznych(f, 1, 2)
    print("Miejsce zerowe funkcji:", root)
