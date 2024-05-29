import numpy as np
import pandas as pd
from sympy import symbols, lambdify, diff, exp
import seaborn as sns
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

expected_x = 0.51574856472875214594

x = symbols('x')  # SymPy symbol
a = -0.7
b = 1.1

# Function that computes f(x) using floats
def f_numeric(x_val):
    return (x_val - 1) * np.exp(-14 * x_val) + x_val**12

# Calculation for the secant method
def calc_xi2(f, xi0, xi1):
    return xi1 - (xi1 - xi0) / (f(xi1) - f(xi0)) * f(xi1)

# Implementing the secant method
def secant_method(f, x0, x1, stop_criterion):
    xi0 = 0.0
    xi1 = x0
    xi2 = x1
    iters = 0

    while not stop_criterion(f, xi1, xi2):
        xi2, xi1, xi0 = calc_xi2(f, xi2, xi1), xi2, xi1
        iters += 1

    return xi2, iters

# Function to iterate through initial guesses and tolerance values
def calculate(a, b, x1, stop_criterion_init, 𝜌_list, step=0.1):
    if not a <= x1 <= b:
        raise Exception(f'x1={x1} is not between {a} and {b}')

    if x1 - a < b - x1:
        a = x1
    else:
        b = x1

    n = int((b - a) / step + 0.5)
    x0_list = [a + step * i for i in range(n)]
    if a == x1:
        x0_list.append(b)

    stop_criterions = [stop_criterion_init(𝜌) for 𝜌 in 𝜌_list]

    df = pd.DataFrame(
        columns=𝜌_list,
        index=[(round(min(x0, x1), 6), round(max(x0, x1), 6)) for x0 in x0_list]
    )

    for i, x0 in enumerate(x0_list):
        for j, stop_criterion in enumerate(stop_criterions):
            df.iloc[i, j] = secant_method(f_numeric, x0, x1, stop_criterion)

    return df

# Visualization function
def show_heatmap(df, annot=True, norm=None, xlabel='x', ylabel='y', title='', **kwargs):
    fig, ax = plt.subplots(figsize=(15, 10))
    sns.heatmap(df, cmap="YlGnBu", annot=annot, norm=norm, **kwargs)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_title(title, fontsize=20)
    plt.show()

# Example usage
stop_criterion_init1 = lambda 𝜌: lambda f, x_prev, x_curr: abs(x_curr - x_prev) < 𝜌
df1 = calculate(a, b, a, stop_criterion_init1, [0.01, 0.001, 0.0001, 1e-5, 1e-7, 1e-10, 1e-15], step=0.1)
show_heatmap(df1.map(lambda cell: cell[1]), xlabel='Wartość rho', ylabel='Punkty startowe x0, x1', title='Liczby iteracji')

# Error heatmap
err_df = df1.map(lambda cell: abs(cell[0] - expected_x))
show_heatmap(err_df, xlabel='Wartość rho', ylabel='Punkty startowe x0, x1', title='Błędy', norm=LogNorm(vmin=1e-15, vmax=1))
