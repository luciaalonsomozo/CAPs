import numpy as np
import matplotlib.pyplot as plt

# Parámetros fijos
b = -3
k = -10

# Número total de iteraciones
iterations = 20

# Valor inicial
x0 = 1e-40

# Crear lista para almacenar la órbita
x_values = []

# Iterar la función
x = x0
for _ in range(iterations):
    x = b + x - k / (1 + np.exp(x))
    x_values.append(x)

# Crear figura
plt.figure(figsize=(12, 6))
plt.plot(range(iterations), x_values, '.-')
plt.title(f"Evolución de xₙ para sistema dinámico con b={b} y k={k}")
plt.xlabel("n (iteración)")
plt.ylabel("xₙ")
plt.grid(True)
plt.show()
