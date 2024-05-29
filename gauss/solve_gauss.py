# Python3 program to demonstrate working of
# Gaussian Elimination method

# wykazać że zasadnicza różnica występuje w czasie pomiędzy metodą gaussa i thomasa

# uwarunkowanie układu równań

# function to get matrix content
def solve_gauss(AB):
    N = len(AB)
    singular_flag = convert_to_row_echelon(AB)

    if (singular_flag != -1):

        print("Singular Matrix.")

        if (AB[singular_flag][N]):
            print("Inconsistent System.")
        else:
            print("May have infinitely many solutions.")

        return

    return back_substitute(AB)

def swap_row(AB, i, j):
    N = len(AB)
    for k in range(N + 1):

        temp = AB[i][k]
        AB[i][k] = AB[j][k]
        AB[j][k] = temp

def convert_to_row_echelon(AB):
	N = len(AB)
	for k in range(N):
	
		i_max = k
		v_max = AB[i_max][k]

		for i in range(k + 1, N):
			if (abs(AB[i][k]) > v_max):
				v_max = AB[i][k]
				i_max = i
                        
		if AB[k][i_max] == 0:
			return k

		if (i_max != k):
			swap_row(AB, k, i_max)

		for i in range(k + 1, N):

			coeff = AB[i][k]/AB[k][k]

			for j in range(k + 1, N + 1):
				AB[i][j] -= AB[k][j]*coeff

			AB[i][k] = 0

		print(AB)

	print(AB)
	return -1

def back_substitute(AB):
    N = len(AB)
    x = [None for _ in range(N)]

    for i in range(N-1, -1, -1):
        x[i] = AB[i][N]

        for j in range(i + 1, N):
            x[i] -= AB[i][j]*x[j]

        x[i] = (x[i]/AB[i][i])

    print("\nSolution for the system:")
    for i in range(N):
        print("{:.8f}".format(x[i])) 
    return x


if __name__ == "__main__":
	AB = [
		[3.0, 2.0, -4.0, 3.0],
		[2.0, 3.0, 3.0, 15.0], 
		[5.0, -3, 1.0, 14.0]
		]
	solve_gauss(AB)

