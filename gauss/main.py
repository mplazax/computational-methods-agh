import random
from gauss.multiply_matrices import multiply_matrices
from gauss.solve_gauss import solve_gauss
from gauss.solve_thomas import solve_thomas


class Solver:


    def __init__(self, n):
        self.A = [[0 for _ in range(n)] for _ in range(n)]
        self.n = n
        self.m = 3
        self.k = 5

    def get_ab(self, A):
        x = random.choices([-1, 1], k=self.n)
        x = [[x[i]] for i in range(self.n)]

        b = multiply_matrices(A, x)

        AB = [[A[i][j] for j in range(self.n)] + [b[i]] for i in range(self.n)]
        return AB


    def solve_gauss(self, A):
        return(solve_gauss(self.get_ab(A)))

    def solve_thomas(self, A):
        return solve_thomas(self.get_ab(A))

    def set_n(self, new_n):
        self.n = new_n
        self.A = [[0 for _ in range(self.n)] for _ in range(self.n)]

    def first(self, n: int):
        self.set_n(n)
        self.A = [[1 for _ in range(n)]] + [[1 / (i + 1 + j + 1 - 1) for j in range(n)] for i in range(1, n)]

    def second(self, n):
        self.set_n(n)
        for i in range(self.n):
            self.n = n
            for j in range(i, self.n):
                self.A[i][j] = 2 * (i + 1) / (j + 1)

            for j in range(i):
                self.A[i][j] = self.A[j][i]

    def third_a(self, n):
        self.set_n(n)
        for i in range(n):
            self.A[i][i] = self.k
            if i < n - 1:
                self.A[i][i+1] = 1 / (i + 1 + self.m)

            if i > 0:
                self.A[i][i-1] = self.k / (i + 1 + self.m + 1)

            for j in range(i - 1):
                self.A[i][j] = 0

            for j in range(i + 2, n):
                self.A[i][j] = 0

    def third_b(self, n):
        self.set_n(n)
        for i in range(n):
            self.A[i][i] = -self.m * (i + 1) * self.k

            if i < n - 1:
                self.A[i][i+1] = i + 1

            if 0 < i:
                self.A[i][i-1] = self.m / (i + 1)

            for j in range(i - 1):
                self.A[i][j] = 0

            for j in range(i + 2, n):
                self.A[i][j] = 0


# TODO dodac mierzenie czasu obliczen

