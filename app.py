from flask import Flask
from flask import render_template
from flask import request
from needleman_wunsch import *
import matplotlib
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.cluster.hierarchy import dendrogram, linkage, average
import subprocess
import numpy as np
import networkx as nx
from estructure_secondary import predict_secondary_structure, traceback, plotData
from Bio.Align import substitution_matrices
from Bio import pairwise2

# app = Flask(__name__)
app = Flask(__name__, static_url_path="/static", static_folder='./static')

DEBUG = True

nw_path = "./nw.exe"
sw_path = "./sw.exe"
sa_path = './sa.exe'
nw_alingments = []
sw_alingments = []
sa_alignments = []

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    # Initialize the matrix
    n = len(seq1)
    m = len(seq2)
    F = [[0] * (m + 1) for _ in range(n + 1)]

    # Fill the first row and column
    for i in range(1, n + 1):
        F[i][0] = F[i - 1][0] + gap
    for j in range(1, m + 1):
        F[0][j] = F[0][j - 1] + gap

    # Fill the rest of the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = F[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete_score = F[i - 1][j] + gap
            insert_score = F[i][j - 1] + gap
            F[i][j] = max(match_score, delete_score, insert_score)

    # The final alignment score is in the bottom-right cell
    return F[n][m]

def string_type(string):
    adn = True
    arn = True
    prot = True

    for a in string:
        if adn and a not in "ACGT":
            adn = False
        if arn and a not in "ACGU":
            arn = False
        if prot and a not in "ACDEFGHIKLMNPQRSTVWY":
            prot = False
        
    if adn: 
        return "ADN"
    elif arn:
        return "ARN"
    elif prot:
        return "PROTEINS"

    return None

def adn_to_arn_transcription(adn_string):
    arn_transcription = ''
    for a in adn_string:
        if a == 'C':
            arn_transcription += 'G'
        elif a == 'G':
            arn_transcription += 'C'
        elif a == 'T':
            arn_transcription += 'A'
        elif a == 'A':
            arn_transcription += 'U'
        
    return arn_transcription

class Chunk:
    def __init__(self, seq, color):
        self.seq = seq
        self.type = color


def get_dot_matrix(seq1, seq2):
    matrix = []
    for i in range(len(seq1)):
        matrix.append([])
        for j in range(len(seq2)):
            matrix[i].append(" ")
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
                matrix[i][j] = "*"
    return matrix


@app.route('/', methods=['POST', 'GET'])
def home():
    ret = False
    if request.method == 'POST':

        return render_template('home.html')

    else:
        return render_template('home.html')

@app.route('/basics', methods=['POST', 'GET'])
def basics():
    ret = False
    if request.method == 'POST':
        seq = request.form['m']
        seq_type = string_type(seq)
        seq_len = len(seq)
        transcription = None
        if seq_type == "ADN":
            transcription = adn_to_arn_transcription(seq)
        ret = True
        return render_template('basics.html', seq=seq, ret=ret, seq_type=seq_type, n_elements=seq_len, transcription=transcription)
    

    else:
        return render_template('basics.html', seq="")


@app.route('/pairwise', methods=['POST', 'GET'])
def pairwise():
    ret = False
    if request.method == 'POST':
        if request.form['algorithm'] == "nw":
            m = request.form['m']
            n = request.form['n']
            match = request.form['match']
            mismatch = request.form['mismatch']
            indel = request.form['indel']
            args = [m, n, match, mismatch, indel]
            subprocess.run([nw_path] + args)
            nw_alingments = []
            with open("nw.txt", "r") as nw_file:
                score = int(nw_file.readline())
                n_alingments = int(nw_file.readline())
                cur = -1
                i = 0
                for line in nw_file.readlines():
                    if line == "#\n":
                        cur += 1;
                        i += 1
                        nw_alingments.append([[],[],[]])
                        continue
                    thing = line.split()
                    if cur >= len(nw_alingments):
                        break
                    if len(thing) < 4:
                        continue
                    nw_alingments[cur][0].append(Chunk(thing[0], thing[3]))
                    nw_alingments[cur][1].append(Chunk(thing[1], thing[3]))
                    nw_alingments[cur][2].append(Chunk(thing[2].replace('-', ' '), thing[3]))
                    i += 1


            ret = True
            print(f"score {score}")
            print(f"Alignments: {n_alingments}")

            return render_template('index.html', seq1=m, seq2=n, match=match, mismatch=mismatch, gap=indel, ret=ret, algorithm="nw", nw_alingments=nw_alingments, score=score, n_alingments=n_alingments)
        elif request.form['algorithm'] == "sw":
            m = request.form['m']
            n = request.form['n']
            match = request.form['match']
            mismatch = request.form['mismatch']
            indel = request.form['indel']
            args = [m, n, match, mismatch, indel]
            subprocess.run([sw_path] + args)
            sw_alingments = []
            with open("sw.txt", "r") as sw_file:
                score = int(sw_file.readline())
                n_alingments = int(sw_file.readline())
                cur = -1
                for line in sw_file.readlines():
                    if line == "#\n":
                        cur += 1;
                        sw_alingments.append([[],[],[]])
                        continue
                    thing = line.split()
                    if cur >= len(sw_alingments):
                        break
                    if len(thing) < 4:
                        continue
                    sw_alingments[cur][0].append(Chunk(thing[0], thing[3]))
                    sw_alingments[cur][1].append(Chunk(thing[1], thing[3]))
                    sw_alingments[cur][2].append(Chunk(thing[2].replace('-', ' '), thing[3]))
            ret = True
            print(f"score {score}")
            print(f"Alignments: {n_alingments}")

            return render_template('index.html', seq1=m, seq2=n, match=match, mismatch=mismatch, gap=indel, ret=ret, algorithm="sw", sw_alingments=sw_alingments, score=score, n_alingments=n_alingments)
        elif request.form['algorithm'] == "dm":
            m = request.form['m']
            n = request.form['n']
            dot_matrix = get_dot_matrix(m, n)
            n_array = []
            for a in n:
                n_array.append(a)
            dot_matrix.insert(0, n_array)

            for i in range(len(dot_matrix)):
                if i == 0:
                    dot_matrix[i].insert(0, " ")
                else:
                    dot_matrix[i].insert(0, m[i - 1])

            return render_template('index.html', seq1=m, seq2=n, match="1", mismatch="-1", gap="-2", ret=ret, algorithm="dm", dot_matrix=dot_matrix)
        else:
            m = request.form['m']
            n = request.form['n']
            match = int(request.form['match'])
            mismatch = int(request.form['mismatch'])
            indel = int(request.form['indel'])

            fasta_alignments = []

            matrix = substitution_matrices.load('blosum62')

            global_alignments = pairwise2.align.globalds(m, n, matrix, indel, indel)
            
            for alignment in global_alignments:
                seq1_aligned, seq2_aligned, score, begin, end = alignment
                fasta_alignments.append([f"\n>secuencia1_alineada\n{seq1_aligned}", f"\n>secuencia2_alineada\n{seq2_aligned}"])

            return render_template('index.html', seq1=m, seq2=n, match=match, mismatch=mismatch, gap=indel, ret=True, algorithm="fasta", fasta_alignments=fasta_alignments, score=score, n_alingments=len(fasta_alignments))
    else:
        
        return render_template('index.html', seq1="", seq2="", algorithm="nw", match="1", mismatch="-1", gap="-2")
    
@app.route('/multiple', methods=['POST', 'GET'])
def multiple():
    if request.method == 'POST':
        text = request.form['sequences']
        sequences = text.strip().split('\n')
        algorithm = request.form['algorithm']
        if algorithm == "sa":
            match = request.form['match']
            mismatch = request.form['mismatch']
            gap = request.form['gap']
            args = [str(len(sequences))]
            args.extend(sequences)
            args.extend([match, mismatch, gap])
            print(args)
            subprocess.run([sa_path] + args)
            star_alignments = []    
            with open("sa.txt", "r") as sa_file:
                star = sa_file.readline()
                global_score = int(sa_file.readline())
                for line in sa_file.readlines():
                    star_alignments.append(line.strip())
            return render_template('multiple.html', sequences=text, algorithm="sa", sa_alignments=star_alignments, star=star, score=global_score)
        else:
            match = int(request.form['match'])
            mismatch = int(request.form['mismatch'])
            gap = int(request.form['gap'])
            dist_matrix = []
            for i in range(len(sequences)):
                dist_matrix.append([])
                for j in range(len(sequences)):
                    if i == j:
                        dist_matrix[i].append(0)
                    else:
                        dist_matrix[i].append(100 - needleman_wunsch(sequences[i], sequences[j], match=match, mismatch=mismatch, gap=gap))
            print(dist_matrix)
            labels = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
            labels = labels[:len(dist_matrix)]

            neighbor_joining(dist_matrix, labels)
            return render_template('multiple.html', sequences=text, algorithm="nj")

    else:
        return render_template('multiple.html', sequences="")
    
        

    
@app.route('/cluster', methods=['POST', 'GET'])
def cluster():
    if request.method == 'POST':
        mode = request.form['mode']
        dist_matrix = request.form['m']

        rows = dist_matrix.strip().split('\n')
        matrix = np.zeros((len(rows), len(rows)))
        i = 0
        for row in rows:
            nums = row.strip().split()
            j = 0
            for num in nums:
                matrix[i][j] = int(num)
                matrix[j][i] = int(num)

                j += 1
            i += 1

        print(matrix)

        upper_triangular_indices = np.triu_indices_from(matrix, k=1)
        upper_triangular_values = matrix[upper_triangular_indices]

        # Step 4: Condense the matrix into a 1D array
        condensed_array = upper_triangular_values

        Z = linkage(condensed_array, method=mode)

        plt.figure(figsize=(8,7))   
        dendrogram(Z)
        plt.title('Dendrograma')
        plt.xlabel('Indice del cluster')
        plt.ylabel('Distancias')
        plt.savefig("static/cluster.png")
        return render_template('clustering.html', mode=mode, matrix=dist_matrix)

    else:
        
        return render_template('clustering.html', mode="single", matrix="")
    


@app.route('/phylogeny', methods=['POST', 'GET'])
def phylogeny():
    if request.method == 'POST':
        dist_matrix = request.form['dm']

        rows = dist_matrix.strip().split('\n')
        matrix = np.zeros((len(rows), len(rows)))
        i = 0
        for row in rows:
            nums = row.strip().split()
            j = 0
            for num in nums:
                matrix[i, j] = int(num)
                matrix[j, i] = int(num)

                j += 1
            i += 1

        print(matrix)
        f = open('graph.dot','w')
        f.write("graph {\n")
        labels = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
        labels = labels[:len(matrix)]
        original_labels = labels[:]

        labels_dict = {}
        for label in original_labels:
            f.write(label + ';\n')
            labels_dict[label] = label

        all_labels = labels[:]
        node_heights = [0 for _ in range(len(labels))]

        distance_matrix = matrix.copy()
        original_matrix = matrix.copy()
        for k in range(len(matrix) - 1):
            min_edge = 100000
            min_row = 0
            min_col = 0
            for i in range(len(distance_matrix)):
                for j in range(i):
                    if min_edge > distance_matrix[i, j]:
                        min_edge = distance_matrix[i, j]
                        min_row = i;
                        min_col = j;
            new_labels = labels[:]
            new_labels.pop(min_row)
            node = new_labels[min_col] + labels[min_row]

            all_labels.append(node)
            node_heights.append(min_edge / 2)

            f.write(node + ';\n')


            f.write(node + " -- " + new_labels[min_col] + ' [label="' + str(min_edge / 2 - node_heights[all_labels.index(new_labels[min_col])]) + '"];\n')
            f.write(node + " -- " + labels[min_row] + ' [label="' + str(min_edge / 2 - node_heights[all_labels.index(labels[min_row])]) + '"];\n')

            new_labels[min_col] = new_labels[min_col] + labels[min_row]

            new_distance_matrix = np.zeros((len(new_labels), len(new_labels)))

            for i in range(len(new_labels)):
                for j in range(i):
                    ii = 0
                    jj = 0
                    if i >= min_row:
                        ii = i + 1
                    else:
                        ii = i
                    if j >= min_row:
                        jj = j + 1
                    else:
                        jj = j

                    if i == min_col:
                        new_distance_matrix[i, j] = (distance_matrix[ii, jj] + distance_matrix[min_row, jj]) / 2
                    elif j == min_col:
                        new_distance_matrix[i, j] = (distance_matrix[ii, jj] + distance_matrix[ii, min_row]) / 2
                    else:
                        new_distance_matrix[i, j] = distance_matrix[ii, jj]
                    new_distance_matrix[j, i] = new_distance_matrix[i, j]

            labels = new_labels[:]
            distance_matrix = np.copy(new_distance_matrix) 

        f.write('}')
        f.close()

        
        args = ["Graphviz/bin/dot", "-Tpng", "graph.dot","-o","static/arbol.png"]

        subprocess.run(args)

        return render_template('filogenia.html', matrix=dist_matrix)

    else:
        return render_template('filogenia.html', matrix="")


@app.route('/secondary', methods=['POST', 'GET'])
def secondary():
    if request.method == 'POST':
        sequence = request.form['sequence']
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

        return render_template('secondary.html', sequence=sequence)
    else:
        return render_template('secondary.html', sequence="")


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
        dist = "{:.2f}".format(child.distance)
        dist = float(dist)
        edges.append((node.name, child.name, dist))

        pre_order_print(child, edges)


def neighbor_joining(matrix, labels):
    nodes = []
    all_labels = {}
    nodes_to_label = []


    for i, label in enumerate(labels):
        nodes.append(Node(label, i))
        all_labels[label] = label
        nodes_to_label.append(label)

    while len(nodes) > 2:
        n = len(nodes)
        total_d = [0 for _ in range(n)]
        r = [0 for _ in range(n)]
        d_s = []
        #print(matrix)
        #print(len(nodes))
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

        all_labels[new_node.name] = new_node.name

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


    G = nx.Graph()
    edges = []
    pre_order_print(nodes, edges)

    # Create a dictionary for labels, excluding certain nodes
    labels = {node: label for node, label in all_labels.items() if node in nodes_to_label}
    G.add_weighted_edges_from((u, v, w) for u, v, w in edges)

    # Draw the graph
    pos = nx.circular_layout(G) # Position the nodes using the spring layout
    nx.draw(G, pos, labels=labels, with_labels=True, node_color='lightblue', node_size=1000, font_size=16, font_weight='bold', edge_color='gray')

    # Draw edge labels
    edge_labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red')

    # Show the plot
    plt.savefig('static/neigh_join.png')
    return nodes



if __name__ == "__main__":
    app.run(host='0.0.0.0', port='5000', debug=True)

