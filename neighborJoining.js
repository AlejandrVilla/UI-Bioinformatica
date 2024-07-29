// utils/neighborJoining.js
export function neighborJoining(matrix, labels) {
    let nodes = labels.map((label, i) => ({ name: label, index: i, children: [] }));
    while (nodes.length > 2) {
        const n = nodes.length;
        const totalD = Array(n).fill(0);
        const R = Array(n).fill(0);
        let D_s = [];
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                totalD[i] += matrix[i][j];
            }
        }

        for(let i = 0 ;  i < n; i++){
            let row = [];
            for(let j = 0 ; j < n; j++){
                if(i != j)
                    row[j] = (n - 2) * matrix[i][j] - totalD[i] - totalD[j];
                else
                    row[j] = 0;
            }
            D_s.push(row);
        }

        let minVal = Infinity;
        let minIndices = [-1, -1];

        for (let i = 0; i < D_s.length; i++) {
            for (let j = 0; j < D_s[i].length; j++) {
                if (D_s[i][j] < minVal) {
                    minVal = D_s[i][j];
                    minIndices = [i, j];
                }
            }
        }
    
        let a = minIndices[0];
        let b = minIndices[1];
        let delta = (totalD[a] - totalD[b]) / (n - 2);

        const distAtoB = matrix[a][b];
        const newDistA = (distAtoB + delta) / 2;
        const newDistB = distAtoB - newDistA;

        const newNode = {
            name: `Node_${nodes.length}`,
            children: [
                { ...nodes[a], distance: newDistA },
                { ...nodes[b], distance: newDistB },
            ],
            index: nodes.length,
        };
        let newMatrix = [];
        let newRow = [0];
        for (let i = 0; i < n; i++) {
            if (i !== a && i !== b) {
                const newDist = (matrix[i][a] + matrix[i][b] - distAtoB) / 2;
                newRow.push(newDist);
            }
        }
        newMatrix.push(newRow);
        for (let i = 0; i < n; i++) {
            if (i !== a && i !== b) {
                const newRow = [];
                for (let j = i; j < n; j++) {
                    if (j !== a && j !== b) {
                        newRow.push(matrix[i][j]);
                    }
                }
                newMatrix.push(newRow);
            }
        }
        let size = newMatrix.length;
        let completeMatrix = Array.from({ length: size }, () => Array(size).fill(0));
        for (let i = 0; i < size; i++) {
            let c = i;
            for (let j = 0; j < newMatrix[i].length; j++) {
                completeMatrix[i][c] = newMatrix[i][j];
                c++;
            }
        }
        for (let i = 0; i < size; i++) {
            for (let j = 0; j < i; j++) {
                completeMatrix[i][j] = completeMatrix[j][i];
            }
        }
        newMatrix = completeMatrix;

        matrix = newMatrix;
        nodes = nodes.filter((_, i) => i !== a && i !== b);
        nodes.unshift(newNode);
    }
    const aux = nodes[1];
    nodes = nodes[0];
    nodes.children.push({...aux, distance: matrix[0][1]})
    console.log(nodes)
    return nodes;
  }
  