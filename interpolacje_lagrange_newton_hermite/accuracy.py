def calculate_accuracy(x_range, interpolation, f):
    return max(abs(interpolation(x) - f(x)) for x in x_range)


def calculate_variance(x_range, interpolation, f):
    return sum((interpolation(x) - f(x)) ** 2 for x in x_range) / len(x_range)
