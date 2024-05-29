from znajdz_miejsca_zerowe import znajdz_miejsca_zerowe


def main():
    tryb_pracy = int(input("Podaj tryb pracy (0 - wartości domyślne): "))
    if tryb_pracy == 0:
        print("Miejsca zerowe dla metody newtona: ", znajdz_miejsca_zerowe(metoda="newton"))
        print("Miejsca zerowe dla metody siecznych: ", znajdz_miejsca_zerowe(metoda="sieczne"))

    while True:
        tol = float(input("Podaj dokładność: "))
        max_iter = int(input("Podaj maks. liczbę iteracji: "))
        n = int(input("Podaj liczbę przedziałów: "))
        newton = znajdz_miejsca_zerowe(n=n, max_iter=max_iter, tol=tol, metoda="newton")
        print("Miejsca zerowe dla metody newtona: ", newton)
        print("Miejsca zerowe dla metody siecznych: ", znajdz_miejsca_zerowe(n=n, max_iter=max_iter, tol=tol, metoda="sieczne"))
    


if __name__ == "__main__":
    main()
