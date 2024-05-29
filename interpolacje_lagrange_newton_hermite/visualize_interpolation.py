import numpy as np
from matplotlib import pyplot as plt

from accuracy import calculate_accuracy
from create_nodes import create_nodes
from get_interp_fun import get_interp_fun
from hermite import hermite_interpolation, hermite_coefficients
from lagrange import lagrange_interpolation
from newton import newton_interpolation


def visualize_interpolation(f, a, b, n, method, chebyshev=False, df=None):
    x_points = create_nodes(a, b, n, chebyshev)
    y_points = [f(x) for x in x_points]
    interp_fun, whos = get_interp_fun(x_points, y_points, method, df)
    x_range = np.linspace(a, b, 1000)
    label_interpolation = "Interpolacja " + whos + f" (Czebyszew: {chebyshev})"
    y_interp = [interp_fun(x) for x in x_range]

    plt.plot(x_range, y_interp, label=label_interpolation)
    plt.scatter(x_points, y_points, label=f"Wezly ({n})")
    plt.plot(x_range, [f(x) for x in x_range], label="Oryginalna funkcja", linestyle='--')
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.legend()

    error = calculate_accuracy(x_range, interp_fun, f)

    plt.annotate(f'Blad: {error:.2f}', xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, color='red')
    plt.show()

