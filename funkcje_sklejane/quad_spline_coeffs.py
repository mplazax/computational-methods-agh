import numpy as np
from specs import df


def quad_spline_coeffs(x_points, y_points, clamped):
    n = len(x_points) - 1
    num_coeffs = 3 * n
    A = np.zeros((num_coeffs, num_coeffs))
    b = np.zeros(num_coeffs)

    # Równość funkcji w punktach
    for i in range(n):
        xi, xi1 = x_points[i], x_points[i+1]
        yi, yi1 = y_points[i], y_points[i+1]
        # Set up the polynomial at knots
        # S_i(xi) = y_i
        A[2*i][3*i] = xi**2
        A[2*i][3*i+1] = xi
        A[2*i][3*i+2] = 1
        b[2*i] = yi
        # S_i(xi+1) = y_i+1
        A[2*i+1][3*i] = xi1**2
        A[2*i+1][3*i+1] = xi1
        A[2*i+1][3*i+2] = 1
        b[2*i+1] = yi1

    # Równość pochodnych w punktach
    for i in range(1, n):
        xi = x_points[i]

        A[2*n+i-1][3*(i-1)] = 2*xi
        A[2*n+i-1][3*(i-1)+1] = 1

        A[2*n+i-1][3*(i)] = -2*xi
        A[2*n+i-1][3*(i)+1] = -1

        b[2*n+i-1] = 0
    
    # clamped boudary condition
    

    if clamped:
        A[3*n-2][0] = 2*x_points[0]
        A[3*n-2][1] = 1
        b[3*n-2] = df(x_points[0])

        A[3*n-1][3*n-3] = 2*x_points[n-1]
        A[3*n-1][3*n-2] = 1
        b[3*n-1] = df(x_points[n-1])

    else:
        A[3*n-2][0] = 2
        b[3*n-2] = 0

        A[3*n-1][3*n-3] = 2
        b[3*n-1] = 0


    # Solve the linear system
    coeffs = np.linalg.solve(A, b)
    return [tuple(coeffs[i:i+3]) for i in range(0, len(coeffs), 3)]

