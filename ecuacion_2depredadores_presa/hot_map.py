import pandas as pd
import matplotlib.pyplot as plt

# Leer CSV
df = pd.read_csv("datos_dp_11_06_aumentando_periodo.csv")

# Limpiar y normalizar columna 'estabilidad'
df['estabilidad'] = df['estabilidad'].astype(str).str.strip().str.upper()

# Convertir período a número
df['periodo'] = pd.to_numeric(df['periodo'], errors='coerce')
df = df.dropna(subset=['periodo'])

print(f"Datos originales: {len(df)} filas")

# Estadísticas generales antes del filtrado
total_estables = len(df[df['estabilidad'] == 'E'])
total_inestables = len(df[df['estabilidad'] == 'I'])

# Índices de los que se mantienen después del filtrado
indices_filtrados = df.sort_values('periodo', ascending=False).drop_duplicates(
    subset=['b', 'k', 'estabilidad']
).index

# Filas eliminadas (las que no están en el filtrado)
df_eliminado = df.loc[~df.index.isin(indices_filtrados)]

eliminadas_estables = len(df_eliminado[df_eliminado['estabilidad'] == 'E'])
eliminadas_inestables = len(df_eliminado[df_eliminado['estabilidad'] == 'I'])

# Filtrar para conservar mayor período por (b, k, estabilidad)
df_filtrado = df.loc[indices_filtrados]

print(f"Datos tras filtrado por (b, k, estabilidad): {len(df_filtrado)} filas")
print("Valores únicos de estabilidad:", df_filtrado['estabilidad'].unique())

print(f"\nEstadísticas:")
print(f" - Total de órbitas estables: {total_estables}")
print(f" - Total de órbitas inestables: {total_inestables}")
print(f" - Órbitas eliminadas con estabilidad 'E': {eliminadas_estables}")
print(f" - Órbitas eliminadas con estabilidad 'I': {eliminadas_inestables}")

# Gráfico: número de órbitas distintas por (b,k)
conteos = df.groupby(['b', 'k']).size().reset_index(name='cantidad')

plt.figure(figsize=(8, 6))
sc = plt.scatter(conteos['b'], conteos['k'], c=conteos['cantidad'], cmap='Reds', s=5)
plt.title("Número de órbitas periódicas distintas por (b, k)")
plt.xlabel("b")
plt.ylabel("k")
plt.grid(True)
plt.colorbar(sc, label='Número de órbitas')
plt.tight_layout()
plt.show()
