import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("predator_prey_model/predator_prey_data.csv")

# Clean data
df['stability'] = df['stability'].astype(str).str.strip().str.upper()
df['period'] = pd.to_numeric(df['period'], errors='coerce')
df = df.dropna(subset=['period'])

# Filter to keep only one orbit per (b, k, period) – ignoring stability
indices_filtrados = df.sort_values('period', ascending=False).drop_duplicates(
    subset=['b', 'k', 'period']
).index
df_filtrado = df.loc[indices_filtrados]

# Count how many different periods appear for each (b, k)
counts = df_filtrado.groupby(['b', 'k']).size().reset_index(name='count')

# Plot
plt.figure(figsize=(8, 6))
sc = plt.scatter(counts['b'], counts['k'], c=counts['count'], cmap='Reds', s=5)
plt.title("Number of distinct periodic orbits per (b, k)")
plt.xlabel("b")
plt.ylabel("k")
plt.grid(True)
plt.colorbar(sc, label='Number of orbits')

plt.tight_layout()
plt.show()
