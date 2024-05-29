import numpy as np


def newton_interpolation(x_points, y_points, x):
    n = len(x_points)
    divided_diff = np.copy(y_points)
    for j in range(1, n):
        for i in range(n-1, j-1, -1):
            divided_diff[i] = (divided_diff[i] - divided_diff[i-1]) / (x_points[i] - x_points[i-j])
    result = divided_diff[-1]
    for i in range(n-2, -1, -1):
        result = result * (x - x_points[i]) + divided_diff[i]
    return result
