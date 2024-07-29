class Node:
    def __init__(self, label, i):
        self.name = label
        self.index = i
        self.children = []
        self.distance = 0

def pre_order_print(node, edges):
    print(node.name)
    print("Hijos de {}, cantidad: {}".format(node.name, len(node.children))) 
    for child in node.children:
        edges.append((node.name, child.name, child.distance))
        pre_order_print(child, edges)


def neighbor_joining(matrix, labels):
    nodes = []
    for i, label in enumerate(labels):
        nodes.append(Node(label, i))

    while len(nodes) > 2:
        n = len(nodes)
        total_d = [0 for _ in range(n)]
        r = [0 for _ in range(n)]
        d_s = []
        print(matrix)
        print(len(nodes))
        for i in range(n):
            for j in range(n):
                total_d[i] += matrix[i][j]

        for i in range(n):
            row = []
            for j in range(n):
                if i != j:
                    row.append((n - 2) * matrix[i][j] - total_d[i] - total_d[j]);
                else:
                    row.append(0);

            d_s.append(row)

        min_val = 100000
        min_indices = [-1,-1]

        for i in range(len(d_s)):
            for j in range(len(d_s[i])):
                if d_s[i][j] < min_val:
                    min_val = d_s[i][j]
                    min_indices = [i, j]
        a = min_indices[0]
        b = min_indices[1]

        delta = (total_d[a] - total_d[b]) / (n - 2);

        dist_a_to_b = matrix[a][b]
        new_dist_a = (dist_a_to_b + delta) / 2
        new_dist_b = dist_a_to_b - new_dist_a

        new_node = Node("Node_{}".format(len(nodes)), len(nodes))
        new_node.children.append(Node(nodes[a].name, nodes[a].index))
        new_node.children[0].children = nodes[a].children
        new_node.children[0].distance = new_dist_a

        new_node.children.append(Node(nodes[b].name, nodes[b].index))
        new_node.children[1].children = nodes[b].children
        new_node.children[1].distance = new_dist_b

        new_matrix = []
        new_row = [0]

        for i in range(n):
            if i != a and i != b:
                new_dist = (matrix[i][a] + matrix[i][b] - dist_a_to_b) / 2;
                new_row.append(new_dist)
        new_matrix.append(new_row)

        for i in range(n):
            if i != a and i != b:
                new_row = []
                for j in range(i, n):
                    if j != a and j != b:
                        new_row.append(matrix[i][j])
                new_matrix.append(new_row)

        size = len(new_matrix)
        complete_matrix = [[0 for _ in range(size)] for _ in range(size)]

        for i in range(size):
            c = i
            for j in range(len(new_matrix[i])):
                complete_matrix[i][c] = new_matrix[i][j]
                c += 1
        
        for i in range(size):
            for j in range(i):
                complete_matrix[i][j] = complete_matrix[j][i]

        new_matrix = complete_matrix

        matrix = new_matrix
        nodes = [node for i, node in enumerate(nodes) if i != a and i != b]
        nodes.insert(0, new_node)

    aux = nodes[1]
    nodes = nodes[0]
    nodes.children.append(Node(aux.name, aux.index))
    nodes.children[len(nodes.children) - 1].children = aux.children
    nodes.children[len(nodes.children) - 1].distance = matrix[0][1]
    print(nodes)


    G = nx.Graph()
    edges = []
    pre_order_print(nodes, edges)


    G.add_weighted_edges_from((u, v, w) for u, v, w in edges)

    # Draw the graph
    pos = nx.spring_layout(G)  # Position the nodes using the spring layout
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=1000, font_size=16, font_weight='bold', edge_color='gray')

    # Draw edge labels
    edge_labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red')

    # Show the plot
    plt.show()
    return nodes

