import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv("predator_prey_model/predator_prey_data.csv")

# Clean and normalize the 'stability' column
df['stability'] = df['stability'].astype(str).str.strip().str.upper()

# Convert 'period' to numeric
df['period'] = pd.to_numeric(df['period'], errors='coerce')
df = df.dropna(subset=['period'])

print(f"Original data: {len(df)} rows")

# General statistics before filtering
total_stable = len(df[df['stability'] == 'S'])
total_unstable = len(df[df['stability'] == 'U'])

# Indices of retained rows after filtering
filtered_indices = df.sort_values('period', ascending=False).drop_duplicates(
    subset=['b', 'k', 'stability']
).index

# Rows that were removed (not in filtered)
df_removed = df.loc[~df.index.isin(filtered_indices)]

removed_stable = len(df_removed[df_removed['stability'] == 'S'])
removed_unstable = len(df_removed[df_removed['stability'] == 'U'])

# Filter to keep the row with the highest period per (b, k, stability)
df_filtered = df.loc[filtered_indices]

print(f"Data after filtering by (b, k, stability): {len(df_filtered)} rows")
print("Unique values of 'stability':", df_filtered['stability'].unique())

print(f"\nStatistics:")
print(f" - Total stable orbits: {total_stable}")
print(f" - Total unstable orbits: {total_unstable}")
print(f" - Removed stable orbits: {removed_stable}")
print(f" - Removed unstable orbits: {removed_unstable}")

# Color palettes for periods
blue_colors = ['lightblue', 'cornflowerblue', 'royalblue', 'navy']
pink_colors = ['lightpink', 'hotpink', 'deeppink', 'mediumvioletred', 'darkmagenta']

# Assign colors by period
colors_by_period = {}
for p in range(2, 6):  # Periods 2 to 5
    colors_by_period[p] = blue_colors[p - 2]
for p in range(6, 11):  # Periods 6 to 10
    colors_by_period[p] = pink_colors[p - 6]

# Create figure with 1 row and 2 columns (subplots)
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

stabilities = ["S", "U"]  # Only these two types

# Set manual limits for the k axis
k_min_manual = -45
k_max_manual = -5

for i, st in enumerate(stabilities):
    ax = axes[i]
    data = df_filtered[df_filtered['stability'] == st]
    if data.empty:
        ax.set_title(f"No data with stability '{st}'")
        ax.axis('off')
        continue

    for p in sorted(colors_by_period):
        subset = data[data['period'] == p]
        if not subset.empty:
            ax.scatter(subset['b'], subset['k'], color=colors_by_period[p], label=f"Period {p}", s=1)

    ax.set_title("Stable orbits" if st == "S" else "Unstable orbits")
    ax.set_xlabel("b")
    ax.set_ylabel("k", rotation=0, ha='right')  # Label 'k' horizontally beside the axis
    ax.legend(fontsize='small')
    ax.grid(True)
    ax.set_ylim(k_min_manual, k_max_manual)

plt.tight_layout()
plt.show()
