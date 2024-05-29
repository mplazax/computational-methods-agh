import math
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
from tabulate import tabulate

def abs_diff(F, f, xs):
    return [abs(F(x) - f(x)) for x in xs]


def max_diff(F, f, xs):
    return max(abs_diff(F, f, xs))

def sum_sq_diff(F, f, xs):
    return sum(d ** 2 for d in abs_diff(F, f, xs))

def calc_error(F, f, a, b, N=1000):
    xs = np.linspace(a, b, N)
    diffs = abs_diff(F, f, xs)
    return {
        'max': max(diffs),
        'sq': sum(x ** 2 for x in diffs)
    }

def show_error(F, fs, a, b, N, *, prec=4):
    for f, label in fs:
        err = calc_error(F, f, a, b, N)
        print(tabulate([
                ('Największy bezwzględny błąd', err['max']),
                ('Suma kwadratów różnic      ', err['sq'])
            ], [
                label
            ], tablefmt='fancy_grid', floatfmt=f'.{prec}e')
        )


def plot_fn(fn, a, b, *, label='', title='Wykres', color='b', step=.1, ax=plt):
    n = int((b - a) / step) + 1
    xs = np.linspace(a, b, n)
    ys = np.vectorize(fn)(xs)
    ax.plot(xs, ys, color, label=label)
    if label: ax.legend(loc='best')

    if ax is plt:
        ax.title(title)
        ax.xlabel('x')
        ax.ylabel('y')
    else:
        ax.title.set_text(title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    ax.grid()

    sns.despine()


class PlotFn:
    def __init__(self, f, color='b', label=''):
        self.f = f
        self.color = color
        self.label = label

def plot_fns(fns: list[PlotFn], a, b, *, title='Wykres', step=.1, ax=plt):
    for fn_obj in fns:
        plot_fn(fn_obj.f, a, b, title=title, step=step, ax=ax, color=fn_obj.color, label=fn_obj.label)


class PlotApprox:
    def __init__(self, approx_method, color='b', label='', args=(), kwargs={}):
        self.im = approx_method
        self.color = color
        self.label = label
        self.args = args
        self.kwargs = kwargs


def rich_plot(fn_obj: 'Funkcja wyjściowa',
              im_objs: list[PlotApprox],
              a, b, n, *,
              step=.01, N=1000,
              nodes_calc_method=np.linspace,
              nodes_color='#000',
              title_comp='Porównanie z wyjściową funkcją',
              title_err='Błędy aproksymacji',
              suptitle='Wykresy',
              show_errors_details=False,
              error_prec=4):
    xs = nodes_calc_method(a, b, n)
    ys = np.vectorize(fn_obj.f)(xs)
    W_objs = [PlotFn(obj.im(xs, ys, *obj.args, **obj.kwargs), obj.color, obj.label) for obj in im_objs]

    fig, ax = plt.subplots(1, 2, figsize=(15, 5))
    fig.suptitle(suptitle, fontsize=20)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    # Compare approximation to the original function
    fns = [fn_obj] + W_objs
    plot_fns(fns, a, b, title=title_comp, step=step, ax=ax[0])
    ax[0].scatter(xs, ys, c=nodes_color)
    ax[0].grid(visible=True)

    # Create errors plot
    xs_ = np.linspace(a, b, N)

    for i, W_obj in enumerate(W_objs):
        diffs = abs_diff(fn_obj.f, W_obj.f, xs_)
        ax[1].scatter(xs_, diffs, s=2, c=im_objs[i].color)
    ax[1].set_title(title_err)
    ax[1].grid(visible=True)
    ax[1].set_xlabel('x')
    ax[1].set_ylabel('y')

    plt.show()

    if show_errors_details:
        show_error(fn_obj.f, [(W.f, W.label) for W in W_objs], a, b, N, prec=error_prec)


def get_suptitle(n, m):
    nodes_text = 'węzłów'

    match n:
        case 1:
            nodes_text = 'węzeł'
        case 2 | 3 | 4:
            nodes_text = 'węzły'

    return f'Wielomian {m}. stopnia, {n} {nodes_text}'

g = lambda x: math.e ** (-math.sin(x)) + math.cos(x)

a = -3 * math.pi
b = 5 * math.pi
x = [a, b]


def gauss_jordan(A, B):
    n = len(B)

    for k in range(n):
        # Partial pivoting
        if np.fabs(A[k, k]) < 1e-12:
            for i in range(k + 1, n):
                if np.fabs(A[i, k]) > np.fabs(A[k, k]):
                    # Swap rows
                    A[[k, i]] = A[[i, k]]
                    B[[k, i]] = B[[i, k]]
                    break
        # Division of the pivot row
        pivot = A[k, k]
        A[k] /= pivot
        B[k] /= pivot
        # Elimination loop
        for i in range(n):
            if i == k or A[i, k] == 0: continue
            factor = A[i, k]
            A[i] -= factor * A[k]
            B[i] -= factor * B[k]

    return B


def least_square_approximation(xs, ys, ws: 'weights', m: 'degree of a plynomial'):
    if len(xs) != len(ys):
        raise ValueError('List of x values and list of y values must have the same length')

    n = len(xs)
    G = np.zeros((m + 1, m + 1), float)
    B = np.zeros(m + 1, float)

    sums = [sum(ws[i] * xs[i] ** k for i in range(n)) for k in range(2 * m + 1)]

    for j in range(m + 1):
        for k in range(m + 1):
            G[j, k] = sums[j + k]

        B[j] = sum(ws[i] * ys[i] * xs[i] ** j for i in range(n))

    A = gauss_jordan(G, B)
    return lambda x: sum(A[i] * x ** i for i in range(m + 1))


def plot(n, m, ws=None):
    if ws is None:
        ws = [1] * n

    rich_plot(
        PlotFn(g, "#777", "Wyjściowa funkcja"),
        [
            PlotApprox(least_square_approximation, '#0070c0', 'Aproksymacja średniokwadratowa', args=(ws, m)),
        ],
        a, b, n,
        show_errors_details=True,
        nodes_color='#073763',
        suptitle=get_suptitle(n, m)
    )


def find_best_approximation(F, a, b, ns, ms, ws, N=1000, max_err=(float('inf'), float('inf'))):
    max_matrix = np.zeros((len(ns), len(ms)))
    sq_matrix = np.zeros((len(ns), len(ms)))

    best = [(0, 0)] * 2
    ns = sorted(ns)
    ms = sorted(ms)

    i = -1
    for n in ns:
        print(f'Calculating for {n} nodes...')
        i += 1
        j = -1
        xs = np.linspace(a, b, n)
        ys = np.vectorize(F)(xs)

        for m in ms:
            if m > n:
                break
            j += 1
            f = least_square_approximation(xs, ys, ws[i][j], m)
            err = calc_error(g, f, a, b, N)

            for k, (matrix, err) in enumerate(zip((max_matrix, sq_matrix), (err['max'], err['sq']))):
                if err > max_err[k]:
                    break
                matrix[i][j] = err
                ii, jj = best[k]
                if err < matrix[ii][jj]:
                    best[k] = i, j

    return {
        'matrix': {
            'max': pd.DataFrame(max_matrix, columns=ms, index=ns),
            'sq': pd.DataFrame(sq_matrix, columns=ms, index=ns)
        },
        'best': {
            'max': (ns[best[0][0]], ms[best[0][1]]),
            'sq': (ns[best[1][0]], ms[best[1][1]])
        }
    }


def show_err_heatmap(df, annot=True):
    plt.figure(figsize=(15, 10))
    s = sns.heatmap(df, cmap="YlGnBu", norm=LogNorm(), annot=annot)
    s.set_xlabel('Stopień wielomianu', fontsize=16)
    s.set_ylabel('Liczba węzłów', fontsize=16)


n_from = 1
n_to = 100

ns = range(n_from, n_to + 1)  # Number of nodes
ms = range(n_from, (n_to - 1) // 2 + 1)  # Degree of the polynomial

res = find_best_approximation(g, a, b, ns, ms)

show_err_heatmap(res['matrix']['max'], annot=False)

show_err_heatmap(res['matrix']['max'], annot=False, norm=LogNorm())

show_err_heatmap(res['matrix']['max'].iloc[0:21, 0:10], annot=False)

show_err_heatmap(res['matrix']['max'].iloc[0:21, 0:10], annot=False, norm=LogNorm())

show_err_heatmap(res['matrix']['sq'], annot=False)

show_err_heatmap(res['matrix']['sq'], annot=False, norm=LogNorm())

show_err_heatmap(res['matrix']['sq'].iloc[0:21, 0:10], annot=False)

show_err_heatmap(res['matrix']['sq'].iloc[0:21, 0:10], annot=False, norm=LogNorm())

n = 7
for m in range(2, 4):
    plot(n, m)

n = 9
for m in range(2, 5):
    plot(n, m)

n = 15
for m in range(2, 8):
    plot(n, m)

n = 22
for m in [*range(3, 11, 3), 10]:
    plot(n, m)

n = 100
for m in (3, 7, 18, 33, 49):
    plot(n, m)

m = 3
for n in range(2 * m + 1, 20, 2):
    plot(n, m)

m = 3
for n in range(103, 106):
    plot(n, m)

m = 7
for n in range(2 * m + 1, 25, 3):
    plot(n, m)

m = 15
for n in [31, *range(2 * m + 20, 151, 20)]:
    plot(n, m)

sq_err_df  = res['matrix']['sq']
max_err_df = res['matrix']['max']

m = 3, 4, 5, 7, 10, 15, 25, 35, 49
n = 7, 10, 15, 20, 25, 35, 40, 50, 60, 75, 85, 100
max_err_df.loc[n, m]

m = 3, 4, 5, 7, 10, 15, 25, 35, 49
n = 7, 10, 15, 20, 25, 35, 40, 50, 60, 75, 85, 100
sq_err_df.loc[n, m]

show_err_heatmap(res['matrix']['max'], annot=False, norm=LogNorm(vmin=10**-1, vmax=10**3))

show_err_heatmap(res['matrix']['sq'], annot=False, norm=LogNorm(vmin=10**-1, vmax=10**7))

plt.figure(figsize=(15, 10))
sns.heatmap(res['matrix']['max'], annot=False, cmap='YlGnBu', vmin=10**-1, vmax=10)
plt.show()

plt.figure(figsize=(15, 10))
sns.heatmap(res['matrix']['sq'], annot=False, cmap='YlGnBu', vmin=10**-1, vmax=2000)
plt.show()





