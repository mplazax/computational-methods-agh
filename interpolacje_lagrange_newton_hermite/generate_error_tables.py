import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from accuracy import calculate_accuracy, calculate_variance
from create_nodes import create_nodes
from get_interp_fun import get_interp_fun
from hermite import hermite_interpolation, hermite_coefficients
from lagrange import lagrange_interpolation
from newton import newton_interpolation


def generate_error_tables(f, a, b, n_max, method, chebyshev=False, df=None):
    data = []

    for n in range(1, n_max + 1):
        x_points = create_nodes(a, b, n, chebyshev)
        y_points = [f(x) for x in x_points]
        interp_fun, whos = get_interp_fun(x_points, y_points, method, df)
        x_range = np.linspace(a, b, 1000)
        accuracy = calculate_accuracy(x_range, interp_fun, f)
        variance = calculate_variance(x_range, interp_fun, f)
        data.append([n, accuracy, variance])

    # Create a DataFrame
    df = pd.DataFrame(data, columns=["n", "Dokładność", "Średni błąd kwadratowy"])

    # Display the DataFrame
    print(df)


# def generate_error_tables(f, a, b, n_max, method, chebyshev=False, df=None):
#     accuracies = []
#     variances = []
#
#     for n in range(1, n_max + 1):
#         x_points = create_nodes(a, b, n, chebyshev)
#         y_points = [f(x) for x in x_points]
#         interp_fun, whos = get_interp_fun(x_points, y_points, method, df)
#         x_range = np.linspace(a, b, 1000)
#         accuracies.append(calculate_accuracy(x_range, interp_fun, f))
#         variances.append(calculate_variance(x_range, interp_fun, f))
#
#     label_accuracy = "Bledy maksymalne interpolacji " + whos + f" (Czebyszew: {chebyshev})"
#     label_variance = "Wariancje interpolacji " + whos + f" (Czebyszew: {chebyshev})"
#
#     plt.plot(range(n_max), accuracies, label=label_accuracy)
#     plt.plot(range(n_max), variances, label=label_variance)
#     plt.xlabel("Liczba wezlow n")
#     plt.ylabel("Wartosc bledu / wariancji")
#     plt.legend()
#     plt.show()


if __name__ == "__main__":
    f = lambda x: np.exp(-np.sin(x)) + np.cos(x)
    df = lambda x: -np.cos(x) * np.exp(-np.sin(x)) - np.sin(x)
    a, b = -3 * np.pi, 5 * np.pi

    while True:
        metoda = str(input("Podaj metode (l/n/h): "))
        metoda = 'lagrange' if metoda == 'l' else 'newton' if metoda == 'n' else 'hermite'
        n = int(input("Podaj n_max (int>0): "))
        chebyshev = str(input("Czy chcesz uzyc wezlow Czebyszewa? (t/n): "))
        chebyshev = True if chebyshev == 't' else False
        generate_error_tables(f, a, b, n, metoda, chebyshev=chebyshev, df=df)
