import pandas as pd
import numpy as np

# Leer CSV
df = pd.read_csv("datos_estabilidad.csv", header=None)

# Asegurarse de que trabajamos con las dos columnas
df_subset = df.iloc[:, :2].copy()

# Tomar solo la parte entera de cada número
for col in df_subset.columns:
    df_subset[col] = df_subset[col].apply(np.floor).astype(int)

# Eliminar duplicados
distintas = df_subset.drop_duplicates()

# Mostrar resultados
print(f"Número total de filas: {len(df_subset)}")
print(f"Número de filas distintas (por parte entera): {len(distintas)}")
