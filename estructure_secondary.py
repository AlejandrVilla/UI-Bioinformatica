import sys
import networkx as nx
import matplotlib.pyplot as plt


def calculate_energy_matrix(sequence, alpha):
    n = len(sequence)
    E = [[0] * n for _ in range(n)]

    for d in range(1, n):
        for i in range(n - d):
            j = i + d
            min_energy = float('inf')

            # Caso 1: E(i+1, j)
            min_energy = min(min_energy, E[i + 1][j])

            # Caso 2: E(i, j-1)
            min_energy = min(min_energy, E[i][j - 1])

            # Caso 3: E(i+1, j-1) + α(ri, rj)
            if sequence[i] + sequence[j] in ['CG', 'GC', 'AU', 'UA']:
                min_energy = min(
                    min_energy, E[i + 1][j - 1] + alpha[sequence[i] + sequence[j]])

            # Caso 4: min{E(i, k-1) + E(k, j)} para i < k ≤ j
            for k in range(i + 1, j):
                min_energy = min(min_energy, E[i][k - 1] + E[k][j])

            E[i][j] = min_energy

    return E


def predict_secondary_structure(sequence, alpha):
    n = len(sequence)
    E = calculate_energy_matrix(sequence, alpha)

    score = E[0][n - 1]

    return E, score


def traceback(sequence, E, alpha):
    n = len(sequence)
    traceback_pairs = ["."] * n
    paired_positions = []

    def traceback_helper(i, j):
        nonlocal traceback_pairs, paired_positions

        if i >= j:
            return

        if sequence[i] + sequence[j] in ['CG', 'GC', 'AU', 'UA'] and E[i][j] == E[i + 1][j - 1] + alpha[sequence[i] + sequence[j]]:
            paired_positions.append((i + 1, j + 1))
            traceback_helper(i + 1, j - 1)
            return

        # Caso 1: E(i+1, j)
        if E[i][j] == E[i + 1][j]:
            traceback_helper(i + 1, j)
            return

        # Caso 2: E(i, j-1)
        if E[i][j] == E[i][j - 1]:
            traceback_helper(i, j - 1)
            return

        else:
            # Caso 4: min{E(i, k-1) + E(k, j)} para i < k ≤ j
            for k in range(i + 1, j):
                if E[i][j] == E[i][k - 1] + E[k][j]:
                    traceback_helper(i, k - 1)
                    traceback_helper(k, j)
                    break

    traceback_helper(0, n - 1)
    return traceback_pairs, paired_positions


def identify_structures(sequence, pairs):
    structures = []
    visited = [False] * len(sequence)
    pairs_dict = {i: j for i, j in pairs}

    # Identificar Bulbos y Hélices
    for i, j in pairs:
        if not visited[i-1] and not visited[j-1]:
            visited[i-1] = visited[j-1] = True
            if abs(i - j) == 1:
                structures.append((i, j, "Bulbo"))
            else:
                structures.append((i, j, "Hélice"))

    unpaired = set(range(1, len(sequence) + 1)) - \
        {i for pair in pairs for i in pair}

    # Identificar Bucles
    for i in unpaired:
        if (i + 1 in unpaired) or (i - 1 in unpaired):
            structures.append((i, i, "Lazo"))

    # Identificar Bucles Internos
    for i, j in pairs:
        if (i - 1 in pairs_dict) and (pairs_dict[i - 1] == j + 1):
            if abs(i - (i - 1)) > 1 and abs(j - (j + 1)) > 1:
                structures.append((i, j, "Bucle Interno"))

    # Identificar Bucle
    for i, j in pairs:
        if (i + 1 in unpaired) and (j - 1 in unpaired) and ((i + 1, j - 1) not in pairs):
            structures.append((i, j, "Bucle"))

    return structures


def plotData(sequence, matches, output_file):
    plt.figure(figsize=(max(4, len(sequence)/8), max(4, len(sequence)/8)))

    node_colors = {
        "A": "#FF6B6B",
        "C": "#FFD93D",
        "G": "#6BCB77",
        "U": "#4D96FF"
    }

    G = nx.Graph()
    for i, c in enumerate(sequence):
        G.add_node(i+1, base=c, color=node_colors.get(c, "#FFFFFF"))

    for i in range(1, len(sequence)):
        G.add_edge(i, i+1, color="#00092C", style="-")

    for e in matches:
        G.add_edge(e[0], e[1], color="#FF5F00", style="--")

    pos = nx.kamada_kawai_layout(G)

    nx.draw_networkx_nodes(G, pos, node_color=nx.get_node_attributes(
        G, "color").values(), node_size=200, edgecolors="#00092C", linewidths=1.5)
    nx.draw_networkx_labels(G, pos, labels=nx.get_node_attributes(G, "base"))

    edge_colors = nx.get_edge_attributes(G, "color").values()
    edge_styles = nx.get_edge_attributes(G, "style").values()
    edges = G.edges()
    nx.draw_networkx_edges(G, pos, edgelist=edges,
                           edge_color=edge_colors, style=edge_styles, width=3)
    structures = identify_structures(sequence, matches)
    labeled_positions = set()

    for structure in structures:
        if len(structure) == 3:
            i, j, name = structure
            if (i, j) not in labeled_positions:
                x = (pos[i][0] + pos[j][0]) / 2
                y = (pos[i][1] + pos[j][1]) / 2
                plt.text(x, y, name, fontsize=8, ha='center',
                         va='center', bbox=dict(facecolor='white', alpha=0.5))
                labeled_positions.add((i, j))
                labeled_positions.add((j, i))
        else:
            i, name = structure
            if i not in labeled_positions:
                x, y = pos[i]
                plt.text(x, y, name, fontsize=8, ha='center',
                         va='center', bbox=dict(facecolor='white', alpha=0.5))
                labeled_positions.add(i)

    for i in range(0, len(sequence), 10):
        plt.annotate(str(i+1), xy=pos[i+1], xytext=(
            5, 5), textcoords='offset points', arrowprops=dict(arrowstyle="->", color='black'))

    plt.axis('off')
    plt.savefig(output_file)
    plt.close()


if __name__ == "__main__":
    sequence = "ACGTGCTGTAGCGGC"
    alpha = {
        'CG': -1,
        'GC': -1,
        'AU': -1,
        'UA': -1,
    }
    filepath = 'static/second_struct.png'
    E, score = predict_secondary_structure(sequence, alpha)

    traceback_pairs = traceback(sequence, E, alpha)
    plotData(sequence, matches=traceback_pairs[1], output_file=filepath)