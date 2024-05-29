from get_values import get_values
from create_nodes import create_nodes
import matplotlib.pyplot as plt
import numpy as np
from specs import a, b, f

def make_plot(title: str) -> None:
    plt.ylabel('Wartości y')
    plt.xlabel('Wartości x')
    plt.title(title)
    plt.legend()
    plt.show()


def main(n: int):
    values = get_values(n)
    
    original = values['f']
    quad_even_ntrl = values["2"]
    cubic_even_ntrl = values["3"]
    quad_cheb_ntrl = values["2h"]
    cubic_cheb_ntrl = values["3h"]
    quad_even_clpd = values["2l"]
    cubic_even_clpd = values["3l"]
    quad_cheb_clpd = values["2hl"]
    cubic_cheb_clpd = values["3hl"]

    cheb_x = create_nodes(n, True)
    even_x = create_nodes(n, False)

    cheb_y = [f(x) for x in cheb_x]
    even_y = [f(x) for x in even_x]

    x_range = np.linspace(a, b, 1000)

    add_to_plot = lambda y, label, color: plt.plot(x_range, y, label=label, color=color)

    plot_quad = lambda y: add_to_plot(y, 'Spline 2 stopnia', 'blue')
    plot_cubic = lambda y: add_to_plot(y, 'Spline 3 stopnia', 'red')
    plot_original = lambda : add_to_plot(original, 'Funkcja oryginalna', 'orange')

    plot_knots = lambda case: plt.plot(cheb_x, cheb_y, color='green') if case == 'cheb' else plt.plot(even_x, even_y, color='green')

    # wykres dla rownomiernych wezlow warunek normal
    plot_quad(quad_even_ntrl)
    plot_cubic(cubic_even_ntrl)
    plot_original()
    plot_knots('even')
    make_plot('Interpolacja dla równomiernych węzłów z warunkiem natural')

    # wykres dla czebyszewa natural
    plot_quad(quad_cheb_ntrl)
    plot_cubic(cubic_cheb_ntrl)
    plot_original()
    plot_knots('cheb')
    make_plot('Interpolacja dla węzłów czebyszewa z warunkiem natural')

    # rownomierne clamped
    plot_quad(quad_even_clpd)
    plot_cubic(cubic_even_clpd)
    plot_original()
    plot_knots('even')
    make_plot('Interpolacja dla węzłów równomiernych z warunkiem clamped')

    # chebyshev clamped
    plot_quad(quad_cheb_clpd)
    plot_cubic(cubic_cheb_clpd)
    plot_original()
    plot_knots('cheb')
    make_plot('Interpolacja dla węzłów czebyszewa z warunkiem clamped')



if __name__ == "__main__":
    n = int(input("Podaj liczbe wezlow: "))
    main(n)

