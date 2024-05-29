import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sympy import symbols, lambdify, diff, exp
import seaborn as sns
from matplotlib.colors import LogNorm

import decimal
from decimal import Decimal, getcontext

from approx_algebr.main import plot_fn

# Set the precision you need
getcontext().prec = 50  # Change this to the desired level of decimal places



# Define the variable and function symbolically
x = symbols('x')
f_sym = (x - 1) * exp(-14 * x) + x**12  # Symbolic definition using SymPy's exp

a =-0.7
b =1.1

def f_numeric(x_val):
    # Make sure to convert all numeric literals to Decimal
    return (x_val - 1) * ((-14) * x_val).exp() + x_val**12

plt.figure(figsize=(15, 5))
plot_fn(f_numeric, a, b, step=.01, color='#0070c0')
def newton_raphson(f, variable, x0, stop_criterion):
    x_prev = Decimal('inf')
    x_curr = Decimal(x0)
    iters = 0

    while not stop_criterion(f, x_prev, x_curr):
        f_curr = f(x_curr)
        f_deriv = (f(x_curr + Decimal('1e-20')) - f_curr) / Decimal('1e-20')  # Numerical derivative
        x_prev = x_curr
        x_curr = x_prev - f_curr / f_deriv
        iters += 1

    return x_curr, iters


def calculate_newton(f, variable, a, b, stop_criterion_init, 𝜌_list, step='0.1'):
    a = Decimal(a)
    b = Decimal(b)
    step = Decimal(step)
    x0_list = [a + step * i for i in range(int((b - a) / step) + 1)]
    stop_criterions = [stop_criterion_init(Decimal(𝜌)) for 𝜌 in 𝜌_list]

    df = pd.DataFrame(columns=𝜌_list, index=[str(x0) for x0 in x0_list])
    for x0 in x0_list:
        for 𝜌, stop_criterion in zip(𝜌_list, stop_criterions):
            result, iterations = newton_raphson(f, variable, x0, stop_criterion)
            df.at[str(x0), 𝜌] = (float(result), iterations)

    return df


def show_heatmap(df, annot=True, norm=None, xlabel='x', ylabel='y', title='', **kwargs):
    ax = plt.figure(figsize=(15, 10))

    # Check for NaNs or infinities and handle them
    if df.isnull().values.any() or np.isinf(df.values).any():
        print("Warning: The data contains NaN or infinite values. Adjusting for visualization.")
        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.fillna(df.min().min())  # Replace NaNs with the smallest non-NaN value

    # Adjust vmin and vmax for LogNorm if very small or very large values are present
    vmin = df.min().min() if df.min().min() > 0 else 1e-10
    vmax = df.max().max()

    norm = LogNorm(vmin=vmin, vmax=vmax) if not norm else norm

    s = sns.heatmap(df, cmap="YlGnBu", annot=annot, norm=norm, mask=df.isnull(), **kwargs)
    s.set_xlabel(xlabel, fontsize=16)
    s.set_ylabel(ylabel, fontsize=16)
    plt.title(title, fontsize=20, y=1.025)
    plt.show()

stop_criterion_init2 = lambda 𝜌: lambda f, _, x_curr: abs(f(x_curr)) < 𝜌
stop_criterion_init1 = lambda 𝜌: lambda f, x_prev, x_curr: abs(x_curr - x_prev) < 𝜌

df1 = calculate_newton(f_numeric, x, '-0.7', '1.1', stop_criterion_init1, ['1e-5', '1e-7', '1e-10', '1e-15'], step='0.1')

show_heatmap(df1.map(lambda cell: cell[1]), xlabel='Wartość rho', ylabel='Punkt startowy x0', title='Liczby iteracji')
expected_x = 0.51574856472875214594
# Apply this function to your DataFrame
err_df = abs(df1.map(lambda cell: cell[0]) - expected_x)
show_heatmap(err_df, xlabel='Wartość rho', ylabel='Punkt startowy x0', title='Błędy',
             norm=LogNorm(vmin=10e-15, vmax=10e-1))  # Adjust vmin and vmax according to your dataset


df2 = calculate_newton(f_numeric, x, '-0.7', '1.1', stop_criterion_init2, ['1e-5', '1e-7', '1e-10', '1e-15'], step='0.1')

show_heatmap(df2.map(lambda cell: cell[1]), xlabel='Wartość rho', ylabel='Punkt startowy x0', title='Liczby iteracji')
expected_x = 0.51574856472875214594
# Apply this function to your DataFrame
err_df = abs(df2.map(lambda cell: cell[0]) - expected_x)
show_heatmap(err_df, xlabel='Wartość rho', ylabel='Punkt startowy x0', title='Błędy',
             norm=LogNorm(vmin=10e-15, vmax=10e-1))  # Adjust vmin and vmax according to your dataset
