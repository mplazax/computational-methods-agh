import specs


def metoda_newtona(f, df, x0, tol=specs.tol, max_iter=specs.max_iter):
    for i in range(max_iter):
        f_x0 = f(x0)
        if abs(f_x0) < tol:
            return x0
        df_x0 = df(x0)
        if df_x0 == 0:
            print("Pochodna równa zero, nie można kontynuować.")
            return None
        x0 = x0 - f_x0 / df_x0
    print("Przekroczono maksymalną liczbę iteracji. Rozwiązanie może nie być dokładne.")
    return x0


if __name__ == "__main__":
    # Przykładowa funkcja i jej pochodna
    def f(x):
        return x**2 - 4

    def df(x):
        return 2*x

    # Test metody Newtona
    root = metoda_newtona(f, df, 2)
    print("Miejsce zerowe funkcji:", root)
