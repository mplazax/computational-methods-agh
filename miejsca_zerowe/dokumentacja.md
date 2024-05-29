## Dokumentacja - Znajdowanie wszystkich miejsc zerowych funkcji w zadanym przedziale

### Opis

Skrypt służy do znajdowania wszystkich miejsc zerowych funkcji na określonym przedziale. Implementuje to poprzez podział przedziału na mniejsze podprzedziały, sprawdzenie zmiany znaku funkcji na końcach tych podprzedziałów i użycie metody siecznych lub metody Newtona dla tych, które potencjalnie zawierają pierwiastki.

### Funkcje

#### `metoda_siecznych(f, x0, x1, tol, max_iter)`

Znajduje przybliżenie miejsca zerowego funkcji metodą siecznych.

##### Parametry

- `f`: funkcja, dla której szukamy miejsca zerowego.
- `x0`, `x1`: dwa początkowe przybliżenia pierwiastka.
- `tol` (opcjonalnie): tolerancja błędu, domyślnie `1e-6`.
- `max_iter` (opcjonalnie): maksymalna liczba iteracji, domyślnie `100`.

##### Zwracana wartość

Zwraca przybliżenie miejsca zerowego funkcji lub `None`, jeśli metoda nie zdoła znaleźć pierwiastka.

#### `metoda_newtona(f, df, x0, tol, max_iter)`

Znajduje przybliżenie miejsca zerowego funkcji metodą Newtona.

##### Parametry

- `f`: funkcja, dla której szukamy miejsca zerowego.
- `df`: pochodna funkcji `f`.
- `x0`: początkowe przybliżenie pierwiastka.
- `tol` (opcjonalnie): tolerancja błędu, domyślnie `1e-6`.
- `max_iter` (opcjonalnie): maksymalna liczba iteracji, domyślnie `100`.

##### Zwracana wartość

Zwraca przybliżenie miejsca zerowego funkcji lub `None`, jeśli metoda nie zdoła znaleźć pierwiastka.

#### `znajdz_miejsca_zerowe(f, df, a, b, n, metoda)`

Wyszukuje wszystkie miejsca zerowe funkcji na zadanym przedziale [a, b].

##### Parametry

- `f`: funkcja, dla której szukamy miejsc zerowych.
- `df`: pochodna funkcji `f`.
- `a`, `b`: przedział, w którym szukamy miejsc zerowych.
- `n`: liczba podprzedziałów na które dzielony jest główny przedział.
- `metoda`: wybór metody rozwiązywania ('newton' lub inne dla metody siecznych).

##### Zwracana wartość

Lista miejsc zerowych funkcji w zadanym przedziale.

### Przykład użycia

```python
def f(x):
    return x**3 - x**2 - x + 1

def df(x):
    return 3*x**2 - 2*x - 1

roots = znajdz_miejsca_zerowe(f, df, -2, 2, 100)
print("Miejsca zerowe funkcji:", roots)
```

### Uwagi

- Należy zachować ostrożność przy wyborze punktów startowych oraz liczby podprzedziałów, aby uniknąć pominięcia pierwiastków.
- Funkcja nie obsługuje automatycznie specjalnych przypadków, takich jak wielokrotne pierwiastki czy funkcje, których pochodna jest równa zero w punkcie zerowym.
