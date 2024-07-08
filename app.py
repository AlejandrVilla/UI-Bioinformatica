from flask import Flask
from flask import render_template
from flask import request
from needleman_wunsch import *
import subprocess

# app = Flask(__name__)
app = Flask(__name__, static_url_path="/static", static_folder='./static')

DEBUG = True

nw_path = "./nw.exe"
sw_path = "./sw.exe"

nw_alingments = []
sw_alingments = []


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
                for line in nw_file.readlines():
                    nw_alingments.append(line.split())
            nw_alingments = [[res[0], res[1], res[2].replace('-', ' ')] for res in nw_alingments]
            ret = True
            print(f"score {score}")
            print(f"Alignments: {n_alingments}")

            return render_template('index.html',ret=ret, algorithm="nw", nw_alingments=nw_alingments, score=score, n_alingments=n_alingments)
        else:
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
                for line in sw_file.readlines():
                    sw_alingments.append(line)
            ret = True
            print(f"score {score}")
            print(f"Alignments: {n_alingments}")

            return render_template('index.html',ret=ret, algorithm="sw", sw_alingments=sw_alingments, score=score, n_alingments=n_alingments)

    else:
        return render_template('index.html')

if __name__ == "__main__":
    app.run(host='0.0.0.0', port='5000', debug=True)