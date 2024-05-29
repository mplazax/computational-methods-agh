import numpy as np


def thomas(a,b,c,d):
    """ A is the tridiagnonal coefficient matrix and d is the RHS matrix"""
    N = len(a)
    cp = np.zeros(N,dtype='float64')
    dp = np.zeros(N,dtype='float64')
    X = np.zeros(N,dtype='float64')
    
    cp[0] = c[0]/b[0]  
    dp[0] = d[0]/b[0]

    for i in np.arange(1,(N),1):
        dnum = b[i] - a[i]*cp[i-1]
        cp[i] = c[i]/dnum
        dp[i] = (d[i]-a[i]*dp[i-1])/dnum
    
    X[(N-1)] = dp[N-1]

    for i in np.arange((N-2),-1,-1):
        X[i] = (dp[i]) - (cp[i])*(X[i+1])
    
    return(X)

def is_tridiagonal(A):
    n = A.shape[0]
    for i in range(n):
        for j in range(n):
            if abs(i-j) > 1 and A[i, j] != 0:
                return False
    return True

def solve_thomas(AB):
    d = [row[-1] for row in AB]
    A = [row[:-1] for row in AB]
    if not is_tridiagonal(A):
        raise Exception("solve_thomas: not a tridiagonal matrix")
    
    a = [A[i][i-1] for i in range(1, len(A))]
    b = [A[i][i] for i in range(len(A))]
    c = [A[i][i+1] for i in range(len(A) - 1)]

    return solve_thomas(a, b, c, d)

if __name__ == "__main__":
    # TODO
    pass