import pandas as pd
import matplotlib.pyplot as plt

# Leer CSV
df = pd.read_csv("datos_logistica_2.csv")

# Limpiar y normalizar columna 'estabilidad'
df['estabilidad'] = df['estabilidad'].astype(str).str.strip().str.upper()

# Convertir 'periodo' y 'nu' a numéricos
df['periodo'] = pd.to_numeric(df['periodo'], errors='coerce')
df['nu'] = pd.to_numeric(df['nu'], errors='coerce')

# Eliminar filas con datos faltantes
df = df.dropna(subset=['periodo', 'nu'])

print(f"Datos originales: {len(df)} filas")

# Filtrar filas únicas para (nu, periodo, estabilidad)
df_filtrado = df.drop_duplicates(subset=['nu', 'periodo', 'estabilidad'])
print(f"Datos tras filtrar duplicados: {len(df_filtrado)} filas")

# Crear figura con 1 fila y 2 columnas para las dos gráficas
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

estabilidades = ['E', 'I']
colores = {'E': 'green', 'I': 'red'}

for ax, est in zip(axes, estabilidades):
    datos = df_filtrado[df_filtrado['estabilidad'] == est]
    if datos.empty:
        ax.set_title(f"No hay datos con estabilidad '{est}'")
        ax.axis('off')
        continue

    # Scatter con puntos en (nu, periodo)
    ax.scatter(datos['nu'], datos['periodo'], color=colores[est], s=5, alpha=0.7)
    ax.set_title(f"Órbitas con estabilidad '{est}'")
    ax.set_xlabel('nu')
    ax.set_ylabel('periodo')
    ax.grid(True)

plt.tight_layout()
plt.show()
