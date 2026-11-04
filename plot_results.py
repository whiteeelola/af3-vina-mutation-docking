import matplotlib.pyplot as plt
import numpy as np

# Data
labels = ['Wild-Type', 'Mutant']
energies = [-7.018, -6.859]

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 5))

# Define colors: blue for Wild-Type, red for Mutant
colors = ['#1f77b4', '#d62728']

# Draw the bars
bars = ax.bar(labels, energies, color=colors, edgecolor='black', width=0.5)

# Add text labels showing the exact values above the bars
for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, yval - 0.2, f'{yval} kcal/mol', 
            ha='center', va='top', fontweight='bold', fontsize=12)

# Axis settings
ax.set_ylabel('Binding Affinity (kcal/mol)', fontsize=12)
ax.set_title('Pyrimethamine Binding Affinity to PfDHFR', fontsize=14, fontweight='bold')
ax.set_ylim(-8.5, 0) # Fix Y-axis scale to make the comparison clear

# Remove top and right spines for a cleaner look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add horizontal grid lines
ax.yaxis.grid(True, linestyle='--', alpha=0.7)

# Save the figure in high resolution and display it
plt.tight_layout()
plt.savefig('binding_affinity_chart.png', dpi=300) 
print("Chart saved as 'binding_affinity_chart.png'")
plt.show()