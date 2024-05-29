from numpy import e, sin, cos, pi

f = lambda x: e**(-sin(x)) + cos(x)
df = lambda x: -cos(x)*e**(-sin(x)) - sin(x)
a = -3 * pi
b = 5 * pi
