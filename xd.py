import numpy as np
import matplotlib.pyplot as plt
from Bio.Align import substitution_matrices
from Bio import pairwise2

# Load the BLOSUM62 substitution matrix
blosum62 = substitution_matrices.load('BLOSUM62')

# Extract the unique amino acids
amino_acids = sorted(list(set([aa for pair in blosum62.keys() for aa in pair])))

# Initialize a matrix for the scores
matrix_size = len(amino_acids)
matrix = np.zeros((matrix_size, matrix_size))

# Fill the matrix with the BLOSUM62 scores
for (aa1, aa2), score in blosum62.items():
    i = amino_acids.index(aa1)
    j = amino_acids.index(aa2)
    matrix[i, j] = score
    matrix[j, i] = score

# Create a heatmap using matplotlib
fig, ax = plt.subplots()
cax = ax.matshow(matrix, cmap='coolwarm')

# Add color bar
fig.colorbar(cax)

# Set the ticks and labels
ax.set_xticks(range(matrix_size))
ax.set_yticks(range(matrix_size))
ax.set_xticklabels(amino_acids)
ax.set_yticklabels(amino_acids)

# Rotate the x labels for better readability
plt.xticks(rotation=90)

# Add title
plt.title("BLOSUM62 Matrix", pad=20)

# Show the heatmap
plt.show()