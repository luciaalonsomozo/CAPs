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

# Paleta de colores por período
colores_azules = ['lightblue', 'cornflowerblue', 'royalblue', 'navy']
colores_rosas = ['lightpink', 'hotpink', 'deeppink', 'mediumvioletred', 'darkmagenta']

# Asignar colores a períodos
colores_por_periodo = {}
for p in range(2, 6):  # 2 a 5
    colores_por_periodo[p] = colores_azules[p - 2]
for p in range(6, 11):  # 6 a 10
    colores_por_periodo[p] = colores_rosas[p - 6]

# --- MODIFICACIONES PRINCIPALES AQUÍ ---

# Crear figura con 2x1 subplots (2 filas, 1 columna)
# Ajustamos el figsize para una mejor visualización vertical
fig, axes = plt.subplots(2, 1, figsize=(7, 10)) # Ancho 7, Alto 10 (ejemplo, ajusta según necesites)

estabilidades = ["E", "I"]  # Solo estas dos

# Definir los límites explícitos para el eje k
k_min_manual = -45
k_max_manual = -5

for i, est in enumerate(estabilidades):
    ax = axes[i] # axes[0] para el primer plot, axes[1] para el segundo
    datos = df_filtrado[df_filtrado['estabilidad'] == est]
    if datos.empty:
        ax.set_title(f"No hay datos con estabilidad '{est}'")
        ax.axis('off')
        continue

    for p in sorted(colores_por_periodo):
        sub = datos[datos['periodo'] == p]
        if not sub.empty:
            ax.scatter(sub['b'], sub['k'], color=colores_por_periodo[p], label=f"Período {p}", s=3)

    if i == 0:
        ax.set_title(f"Órbitas estables")
    else:
        ax.set_title(f"Órbitas inestables")

    ax.set_xlabel("b")
    # Configurar el label 'k' en vertical
    # Usamos rotation=0 y ha='right' para que 'k' se lea horizontalmente al lado del eje
    # Si quieres que 'k' se lea de abajo a arriba, usa rotation=90
    ax.set_ylabel("k", rotation=0, ha='right')

    ax.legend(fontsize='small')
    ax.grid(True)
    
    # Aplicar los límites de k definidos manualmente
    ax.set_ylim(k_min_manual, k_max_manual)

plt.tight_layout()
plt.show()