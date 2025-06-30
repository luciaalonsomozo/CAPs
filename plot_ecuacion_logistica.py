import numpy as np
import matplotlib.pyplot as plt

# --- AÑADE ESTAS LÍNEAS ---
plt.rcParams['text.usetex'] = False  # Asegurarse de que no use LaTeX externo
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}' # Aunque no use LaTeX externo, esto a veces ayuda con mathtext
plt.rcParams['mathtext.fontset'] = 'cm' # Usa la fuente Computer Modern para mathtext
# --------------------------

# Definir el rango de valores para r
r_values = np.arange(3.56995, 4.0, 0.0001)

# Número de iteraciones y las que se van a graficar
iterations = 2000
last = 1000

# Inicializar el array para los valores de x
x = 1e-5 * np.ones(r_values.shape)

# Crear la figura y el eje
plt.figure(figsize=(10, 6))

# Iterar sobre los valores de r
for i in range(iterations):
    x = r_values * x * (1 - x)
    # En las últimas iteraciones, graficar los puntos
    if i >= (iterations - last):
        plt.plot(r_values, x, ',k', alpha=0.75, markersize=2)

# Configurar los ejes
plt.xlim(3.56995, 4)
# plt.title("Diagrama de bifurcación del mapa logístico")
plt.xlabel(r"$\mu$")
plt.ylabel("x")
plt.show()
