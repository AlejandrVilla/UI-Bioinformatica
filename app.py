from flask import Flask
from flask import render_template
from flask import request
from needleman_wunsch import *
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import subprocess
import numpy as np

# app = Flask(__name__)
app = Flask(__name__, static_url_path="/static", static_folder='./static')

DEBUG = True

nw_path = "./nw.exe"
sw_path = "./sw.exe"
sa_path = './sa.exe'
nw_alingments = []
sw_alingments = []
sa_alignments = []


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
def nw():
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

            return render_template('index.html',ret=ret, algorithm="nw", nw_alingments=nw_alingments, score=score, n_alingments=n_alingments)
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

            return render_template('index.html',ret=ret, algorithm="sw", sw_alingments=sw_alingments, score=score, n_alingments=n_alingments)
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

            return render_template('index.html',ret=ret, algorithm="dm", dot_matrix=dot_matrix)

    else:
        
        return render_template('index.html')
    
@app.route('/multiple', methods=['POST', 'GET'])
def multiple():
    if request.method == 'POST':
        text = request.form['sequences']
        sequences = text.strip().split('\n')
        algorithm = request.form['algorithm']
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
        return render_template('multiple.html', algorithm="sa", sa_alignments=star_alignments, star=star, score=global_score)
    else:
        return render_template('multiple.html')
        

    
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
        return render_template('clustering.html')

    else:
        
        return render_template('clustering.html')

if __name__ == "__main__":
    app.run(host='0.0.0.0', port='5000', debug=True)