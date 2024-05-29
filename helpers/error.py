import numpy as np


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