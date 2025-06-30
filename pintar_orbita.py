import numpy as np
import matplotlib.pyplot as plt

# Función logística
def f(x):
    return 3.45 * x * (1 - x)

# Función para graficar el diagrama de telaraña
def cobweb(f, x0, n):
    x = np.linspace(0, 1, 500)
    y = f(x)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(x, y, label=r"$f(x)=3x(1-x)$", color="blue")
    ax.plot(x, x, label="y=x", linestyle="--", color="black")  # Línea identidad

    # Iteraciones
    xn = x0
    for i in range(n):
        
        yn = f(xn)
        if i > 150:
            ax.plot([xn, xn], [xn, yn], color="red")  # Línea vertical
            ax.plot([xn, yn], [yn, yn], color="red")  # Línea horizontal
            # print(n)
        xn = yn

    ax.set_xlabel("$x$")
    ax.set_ylabel("$f(x)$")
    ax.legend()
    ax.set_title(f"Cobweb Plot para $x_0 = {x0}$")
    plt.show()

# Punto inicial arbitrario
cobweb(f, x0=0.97, n=200)
