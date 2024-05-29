import numpy as np

from create_nodes import create_nodes
from hermite import hermite_coefficients, hermite_interpolation
from lagrange import lagrange_interpolation
from newton import newton_interpolation


def get_interp_fun(x_points, y_points, method=None, df=None):
    if method == 'lagrange':
        interp_fun = lambda x: lagrange_interpolation(x_points, y_points, x)
        whos = "Lagrange'a"
    elif method == 'newton':
        interp_fun = lambda x: newton_interpolation(x_points, y_points, x)
        whos = "Newtona"
    elif method == 'hermite':
        df = [df(x) for x in x_points]
        coefficients = hermite_coefficients(x_points, y_points, df)
        interp_fun = lambda x: hermite_interpolation(x_points, coefficients, x)
        whos = "Hermita"
    else:
        interp_fun = lambda x: x
        whos = ""

    return interp_fun, whos
