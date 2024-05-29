import numpy as np
from specs import df


def cubic_spline_coeffs(x_points, y_points, clamped):
    n = len(x_points) - 1
    num_coeffs = 4 * n
    A = np.zeros((num_coeffs, num_coeffs))
    b = np.zeros(num_coeffs)

    # Równość funkcji w punktach
    for i in range(n):
        xi, xi1 = x_points[i], x_points[i+1]
        yi, yi1 = y_points[i], y_points[i+1]
        # Set up the polynomial at knots
        # S_i(xi) = y_i
        A[2*i][4*i] = xi**3
        A[2*i][4*i+1] = xi**2
        A[2*i][4*i+2] = xi
        A[2*i][4*i+3] = 1
        b[2*i] = yi
        # S_i(xi+1) = y_i+1
        A[2*i+1][4*i] = xi1**3
        A[2*i+1][4*i+1] = xi1**2
        A[2*i+1][4*i+2] = xi1
        A[2*i+1][4*i+3] = 1
        b[2*i+1] = yi1
    # od 0 do 2 * n - 1
    
    # Równość pochodnych w punktach
    for i in range(1, n):
        xi = x_points[i]

        A[3*n+i-1][4*(i-1)] = 3*xi*xi
        A[3*n+i-1][4*(i-1)+1] = 2*xi
        A[3*n+i-1][4*(i-1)+2] = 1

        A[3*n+i-1][4*(i)] = -3*xi*xi
        A[3*n+i-1][4*(i)+1] = -2*xi
        A[3*n+i-1][4*(i)+2] = -1

        b[3*n+i-1] = 0
    # od 3 * n do 4 * n - 2

    # Równość drugich pochodnych w punktach
    for i in range(1, n):
        xi= x_points[i]

        A[3*n+i-2][4*(i-1)] = 6*xi
        A[3*n+i-2][4*(i-1)+1] = 2

        A[3*n+i-2][4*(i)] = -6*xi
        A[3*n+i-2][4*(i)+1] = -2

        b[3*n+i-2] = 0
    # zaczyna się od 3n - 1
    # kończy się na 4n - 3

    if clamped:
        x1 = x_points[0]
        A[4*n-2][0] = 3*x1*x1
        A[4*n-2][1] = 2*x1
        A[4*n-2][2] = 1
        b[4*n-2] = df(x1)

        xn = x_points[-1]
        A[4*n-1][-4] = 3*xn*xn
        A[4*n-1][-3] = 2*xn
        A[4*n-1][-2] = 1
        b[4*n-1] = df(xn)

    else:
        # free boundary
        A[4*n-2][0] = 6*x_points[0]
        A[4*n-2][1] = 2
        b[4*n-2] = 0

        A[4*n-1][-4] = 6*x_points[-1]
        A[4*n-1][-3] = 2
        b[4*n-1] = 0

    coeffs = np.linalg.solve(A, b)
    return [tuple(coeffs[i:i+4]) for i in range(0, len(coeffs), 4)]
