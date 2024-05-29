import numpy as np


def hermite_coefficients(x, f, df):
    n = len(x)
    q = np.zeros((2 * n, 2 * n))
    # Wstawienie wartości funkcji i pierwszych pochodnych
    for i in range(n):
        q[2 * i][0] = f[i]
        q[2 * i + 1][0] = f[i]
        q[2 * i + 1][1] = df[i]
        if i > 0:
            q[2 * i][1] = (q[2 * i][0] - q[2 * i - 2][0]) / (x[i] - x[i - 1])

    # Wypełnianie tabeli różnic podzielonych
    for i in range(2, 2 * n):
        for j in range(2, i + 1):
            q[i][j] = (q[i][j - 1] - q[i - 1][j - 1]) / (x[i // 2] - x[(i - j) // 2])

    # Wyciągnięcie współczynników
    coefficients = [q[i][i] for i in range(2 * n)]
    return coefficients


def hermite_interpolation(x, coefficients, z):
    n = len(x)
    # Wielomian wynikowy zaczynamy od ostatniego współczynnika
    p = coefficients[-1]
    # Wielomian tworzymy od końca do początku
    for i in range(2 * n - 2, -1, -1):
        p = coefficients[i] + (z - x[i // 2]) * p
    return p


# Testowanie funkcji
if __name__ == "__main__":
    x = np.array([0, 1, 2], dtype=float)  # węzły interpolacyjne
    f = np.sin(x)  # wartości funkcji w węzłach interpolacyjnych
    df = np.cos(x)  # wartości pochodnych w węzłach interpolacyjnych

    coefficients = hermite_coefficients(x, f, df)
    z = 0.5
    H_z = hermite_interpolation(x, coefficients, z)

    print(f"Wartość interpolacji Hermite'a w punkcie z={z}: {H_z}")
