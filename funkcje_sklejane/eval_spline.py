def eval_spline(x_range, x_points, coefficients):
    
    y_points = [0 for _ in x_range]
    j = 0
    for i in range(len(x_points) - 1):
        while x_points[i] <= x_range[j] < x_points[i+1]:
            x = x_range[j]
            coefs = coefficients[i]
            for k in range(len(coefs)):
                y_points[j] += coefs[k] * (x ** (len(coefs) - 1 - k))
            j += 1
    return y_points
