import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set Seaborn style
sns.set(style="white", font_scale=1.2)  

# Path to your CSV
csv_file = "file path/ affinity.csv"

# Read CSV
df = pd.read_csv(csv_file)

# Drop rows with missing affinity
df = df.dropna(subset=["Minimized_Affinity"])

# Sort by affinity (most negative = strongest binder)
df = df.sort_values(by="Minimized_Affinity")

# Create plot
plt.figure(figsize=(14,7))
bars = plt.bar(df["Ligand"], df["Minimized_Affinity"], color='skyblue', edgecolor='black')

# Highlight top 5 ligands in orange
for i, rect in enumerate(bars):
    if i < 5:
        rect.set_color('orange')

# Add numerical labels on top of bars
for rect in bars:
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width()/2.0, height + 0.05, f'{height:.2f}',
             ha='center', va='bottom', fontsize=10)

# Rotate x-axis labels for readability
plt.xticks(rotation=45, ha='right')

# Axis labels and title
plt.ylabel("Minimized Affinity (kcal/mol)", fontsize=12)
plt.title("Minimized Affinity of Antiviral Compounds to WNV NS5 Protein",
          fontsize=16, weight='bold')

# Remove grid lines for clean look
plt.grid(False)

# Tight layout
plt.tight_layout()
plt.show()
