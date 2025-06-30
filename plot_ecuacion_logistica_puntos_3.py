import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm # Import colormap module

# --- AÑADE ESTAS LÍNEAS ---
plt.rcParams['text.usetex'] = False  # Asegurarse de que no use LaTeX externo
# plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}' # Aunque no use LaTeX externo, esto a veces ayuda con mathtext
plt.rcParams['mathtext.fontset'] = 'cm' # Usa la fuente Computer Modern para mathtext
# --------------------------

# Definir el rango de valores para r
# Adjusted start of r_values to include mu=0.5
r_values = np.arange(3.56995, 4.0, 0.0001)

# Número de iteraciones y las que se van a graficar
iterations = 2000
last = 1000

# Inicializar el array para los valores de x
x = 1e-5 * np.ones(r_values.shape)

# Crear la figura y el eje
plt.figure(figsize=(14, 8)) # Increased figure size again for more data

# Iterar sobre los valores de r para el diagrama de bifurcación
for i in range(iterations):
    x = r_values * x * (1 - x)
    # En las últimas iteraciones, graficar los puntos del diagrama
    if i >= (iterations - last):
        plt.plot(r_values, x, ',k', alpha=0.3, markersize=0.5) # Even lighter background for more overlaid points

# --- TODA LA DATA EXTRAÍDA DE LAS TABLAS LATEX PARA ÓRBITAS ESTABLES ---
stable_orbits_data = [
    
    # From first table (interval 3.56995 < mu < 4)
    (3.6055, 10, [0.9013725844806769, 0.3205291246046368, 0.785242583712226, 0.6080195080502302,
                  0.8593052539924549, 0.4359039925717974, 0.8865625309544296, 0.362603006536723,
                  0.8333106096381467, 0.5008185072158021]),
    (3.6275, 6, [0.7708464170152596, 0.6407696472097291, 0.8349921205547574,
                 0.499797962675559, 0.9068748519292111, 0.30635273103215666]),
    (3.6624, 8, [0.7432272471261291, 0.6989342701099414, 0.7706610919783411, 0.6473020004139705,
                 0.8361336867566139, 0.5018006032834138, 0.9155881258685924, 0.2830541128946623]),
    (3.68172, 11, [0.7250638969792778, 0.7339370471142395, 0.7189421957416284, 0.74394422954138,
                   0.7013353082604845, 0.7711881428644469, 0.6496652335294268, 0.8379606423178432,
                   0.499913529697628, 0.9204299724713629, 0.2696441592094289]),
    (3.6872, 9, [0.2657954282681484, 0.7195505115478071, 0.7440680107198037, 0.7021564684114944,
                 0.771114321076577, 0.6507798302377812, 0.8379731406680769, 0.500626468692557,
                 0.9217985529104166]),
    (3.7020, 7, [0.7038164287080472, 0.7717147122649582, 0.6521854682437301, 0.8397601372116263,
                 0.4981524079974901, 0.9254873628668389, 0.25529178595838636]),
    (3.74, 5, [0.842538311847857, 0.49617646838741175, 0.9349453234664683,
               0.2274763953239883, 0.6572335095050293]),
    (3.7742, 7, [0.2010272647771743, 0.6061942788221429, 0.9009874979539385, 0.33669270775832894,
                 0.8428948347524774, 0.49979138592692884, 0.9435498357474521]),
    (3.82843, 3, [0.5139418448568985, 0.9563633487732167, 0.15976993160815986]),
    (3.8415, 6, [0.1486296807799159, 0.4860991796294752, 0.9596326961720238,
                 0.1488111995779169, 0.4865890772387863, 0.9596840952793992]),
    (3.8442, 6, [0.9610472346930975, 0.14390934682529946, 0.47360331308723286, 0.958371418795577,
                 0.15336682863451812, 0.49915185778191484]),
    (3.9058, 5, [0.8489540057737763, 0.5008450376230533, 0.9764472109128072,
                 0.08982580405252519, 0.3193269943655632])
  
]
# --- FIN DE LA DATA ---

# Choose a colormap for distinct colors.
# 'tab20' has 20 distinct colors. If more than 20 stable orbits are plotted,
# colors will repeat, but it's usually still distinguishable.
colors = [
                    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
                    "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7"]

# Plot each stable orbit
for idx, (mu_val, period_val, x_points) in enumerate(stable_orbits_data):
    # Create an array of mu_val repeated for each x_point
    mu_array = np.full_like(x_points, mu_val, dtype=float) # Ensure float type
    # Plot the points with a distinct color and marker
    plt.plot(mu_array, x_points, 'o', color=colors[idx], markersize=6, # Slightly smaller markers
             label=f'$\mu={mu_val}$ (P={period_val})', alpha=0.9)

# Configurar los ejes
plt.xlim(r_values.min(), r_values.max()) # Set x-limits based on r_values range
plt.xlabel(r"$\mu$")
plt.ylabel("x")
plt.title("Diagrama de bifurcación del mapa logístico con órbitas estables")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize='small') # Smaller font for legend
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()