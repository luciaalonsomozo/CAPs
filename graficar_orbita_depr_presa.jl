module Graficar

export graficar

using Plots
gr()  # Establecer GR como el backend

    function graficar(f, x0, b, k, n_max=50)
        # Inicialización del vector de valores
        x = zeros(n_max)
        x[1] = x0

        # Iterar para calcular x_n+1
        for n in 2:n_max
            x[n] = f(x[n-1], b, k)
        end

        # Graficar los resultados
        # Graficar los resultados con diferentes colores para puntos y líneas
        p = plot(1:n_max, x, xlabel="n", ylabel="x_n", title="Dos depredadores y una presa", linewidth=2, label="", color=:blue)
        scatter!(1:n_max, x, color=:red, label="", marker=:circle, markersize=5)  # Cambia el color de los puntos y el tamañ
        display(p)

        sleep(100)
    end
end