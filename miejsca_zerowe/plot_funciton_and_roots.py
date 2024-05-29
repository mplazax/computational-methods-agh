import matplotlib.pyplot as plt
import numpy as np
# kryterium stopu gdzie jesli zejdzie nachylenie ponizej danego epsilon to zatrzymujemy
# metoda newtona: parametry: kryterium stopu, epsilon kryterium stopu, punkt startowy
# metoda siecznych: ustawienie punktów startowych, ile wykonanych iteracji, kryterium stopu 

def plot_function_and_roots(f, a, b, roots):
    # Define x values for plotting
    x_values = np.linspace(a, b, 400)
    y_values = f(x_values)

    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, label='e^(-sin(x)) + cos(x)')
    plt.scatter(roots, [f(root) for root in roots], color='red', zorder=5, label='Roots')
    plt.title('Function Plot and its Roots')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.axhline(0, color='black',linewidth=0.5)
    plt.axvline(0, color='black',linewidth=0.5)
    plt.grid(True)
    plt.legend()
    plt.show()