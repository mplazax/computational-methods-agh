from numpy import e, sin, cos, pi

f = lambda x: (x-1)*e**(-16*x)+x**12
df = lambda x: e**(-16*x)-16*e**(-16*x)+12*x**11
a = -0.7
b = 1.1
tol=1e-6
max_iter=100
n=100