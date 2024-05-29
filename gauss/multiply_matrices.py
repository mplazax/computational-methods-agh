def multiply_matrices(a, b):
    """Multiply two matrices."""
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]