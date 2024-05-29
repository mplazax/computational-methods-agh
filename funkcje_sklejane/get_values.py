from eval_spline import eval_spline
from cubic_spline_coeffs import cubic_spline_coeffs
from quad_spline_coeffs import quad_spline_coeffs
from specs import f, a, b
from numpy import linspace
from create_nodes import create_nodes


def get_values(n):
    x_range = linspace(a, b, 1000)
    x_uni = create_nodes(n, False)
    x_cheb = create_nodes(n, True)
    y_uni = [f(x) for x in x_uni]
    y_cheb = [f(x) for x in x_cheb]

    vals = {
        "f": [f(x) for x in x_range],
        "3": eval_spline(x_range, x_uni, cubic_spline_coeffs(x_uni, y_uni, False)),
        "2": eval_spline(x_range, x_uni, quad_spline_coeffs(x_uni, y_uni, False)),
        "3h": eval_spline(x_range, x_cheb, cubic_spline_coeffs(x_cheb, y_cheb, False)),
        "2h": eval_spline(x_range, x_cheb, quad_spline_coeffs(x_cheb, y_cheb, False)),
        "3l": eval_spline(x_range, x_uni, cubic_spline_coeffs(x_uni, y_uni, True)),
        "2l": eval_spline(x_range, x_uni, quad_spline_coeffs(x_uni, y_uni, True)),
        "3hl": eval_spline(x_range, x_cheb, cubic_spline_coeffs(x_cheb, y_cheb, True)),
        "2hl": eval_spline(x_range, x_cheb, quad_spline_coeffs(x_cheb, y_cheb, True)),
        }

    return vals
