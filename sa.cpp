#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#define str string
#define vstr vector<str>
#define pstrstr pair<str,str>
#define vi vector<int>
#define vpii vector<pair<int,int>>
#define vvi vector<vi>
#define vvvi vector<vvi>
#define vpstrstr vector<pair<str,str>>
#define vvpstrstr vector<vpstrstr>
#define inf INT32_MAX
using namespace std;

vstr cadenas;
int n_cadenas;
vvi matriz_scores;
int global_max_score=-inf;
int indx_max_score;
vvpstrstr alineamientos;
str star;
vpstrstr alineamiento_multiple;
vstr alineamiento_sol;

void print_cadenas(ofstream& file)
{
    for(auto cadena : cadenas)
        file<<cadena<<"\n";
}

void print_scores(ofstream& file)
{
    for(auto fila : matriz_scores)
    {
        for(auto score : fila)
            file<<score<<" ";
        file<<"\n";
    }
}

void print_alineamientos(ofstream& file)
{
    for(int f=0 ; f<n_cadenas ; ++f)
    {
        for(int c=0 ; c<n_cadenas ; ++c)
        {
            if(f!=c)
            {
                file<<f+1<<" "<<c+1<<"\n";
                file<<alineamientos[f][c].first<<" "<<alineamientos[f][c].second<<"\n";
            }
        }
        file<<"\n";
    }
}

void print_sol(ofstream& file)
{
    for(auto alineamiento : alineamiento_sol)
        file<<alineamiento;
}

// alineamiento global de dos cadenas
void nw(vvi& matriz, vvvi& res, str& cadena1, str& cadena2)
{
    int n = cadena1.size();
    int m = cadena2.size();
    int max_score = 0;
    for(int i=1 ; i<n+1 ; ++i)
    {
        matriz[i][0] = matriz[i-1][0] - 2;
        res[i][0].push_back(1);    // up
    }
    for(int j=1 ; j<m+1 ; ++j)
    {
        matriz[0][j] = matriz[0][j-1] - 2;
        res[0][j].push_back(-1);    // left
    }
    for(int f=1 ; f<n+1 ; ++f)
    {
        for(int c=1 ; c<m+1 ; ++c)
        {
            int value=0;
            if(cadena1[f-1] == cadena2[c-1]) value=1;
            else value=-1;
            int diag, up, left, max_value;
            diag = matriz[f-1][c-1] + value;
            up = matriz[f-1][c] - 2;
            left = matriz[f][c-1] - 2;
            max_value = max(diag, max(up, left));
            // -1 left
            // 1 up
            // 0 diag
            if(diag == max_value) res[f][c].push_back(0);
            if(up == max_value) res[f][c].push_back(1);
            if(left == max_value) res[f][c].push_back(-1);
            matriz[f][c] = max_value;
        }
    }
}

// get path
void get_sol(int f, int c, str cadena1, str cadena2, str sol1, str sol2, pstrstr& alineamiento, vvvi& res, bool ok)
{
    if(f<0 || c<0)
        return;
    if(f==0 && c==0) 
    {
        // push solution
        reverse(sol1.begin(), sol1.end());
        reverse(sol2.begin(), sol2.end());
        alineamiento = make_pair(sol2, sol1);
        ok = true;
        // cout<<solution<<'\n';
        return;
    }
    for(int i=0 ; i<res[f][c].size() && !ok ; ++i)
    {
        // cout<<"size: "<<res[f][c].size()<<" pos: "<<f<<" "<<c<<"\n";
        char caracter1 = ' ';
        char caracter2 = ' ';
        str solution1 = "";
        str solution2 = "";
        if(res[f][c][i] == -1) //left
        {
            caracter1 = cadena2[c-1];
            caracter2 = '_';
            solution1 = sol1 + caracter1;
            solution2 = sol2 + caracter2;
            get_sol(f, c-1, cadena1, cadena2, solution1, solution2, alineamiento, res, ok);
        }
        else if(res[f][c][i] == 1)   // up
        {
            caracter1 = '_';
            caracter2 = cadena1[f-1];
            solution1 = sol1 + caracter1;
            solution2 = sol2 + caracter2;
            get_sol(f-1, c, cadena1, cadena2 , solution1, solution2, alineamiento, res, ok);
        }
        else if(res[f][c][i] == 0)  // diag
        {
            caracter1 = cadena2[c-1];
            caracter2 = cadena1[f-1];
            solution1 = sol1 + caracter1;
            solution2 = sol2 + caracter2;
            get_sol(f-1, c-1, cadena1, cadena2, solution1, solution2, alineamiento, res, ok);
        }
    }
}

// calculo matriz de scores
void process()
{
    for(int i=0 ; i<n_cadenas ; ++i)
    {
        for(int j=i+1 ; j<n_cadenas ; j++)
        {
            str cadena1 = cadenas[i];
            str cadena2 = cadenas[j];
            int n = cadena1.size();
            int m = cadena2.size();
            vvi matriz(n+1, vi(m+1,0));
            vvvi res(n+1, vvi(m+1, vi(0)));
            nw(matriz, res, cadena1, cadena2);
            matriz_scores[i][j] = matriz[n][m];
            matriz_scores[j][i] = matriz[n][m];
            pstrstr alineamiento;
            str sol1="";
            str sol2="";
            get_sol(n, m, cadena1, cadena2, sol1, sol2, alineamiento, res, false);
            // cout<<alineamiento.first<<" "<<alineamiento.second<<"\n";
            alineamientos[i][j] = alineamiento;
            alineamientos[j][i] = make_pair(alineamiento.second, alineamiento.first);
        }
    }
}

// mejor score, indice de la estrella
void max_score()
{
    for(int f=0 ; f<n_cadenas ; ++f)
    {
        int acc = 0;
        for(int c=0 ; c<n_cadenas ; ++c)
            acc += matriz_scores[f][c];
        matriz_scores[f][n_cadenas] = acc;
        if(acc > global_max_score)
        {
            global_max_score = acc;
            indx_max_score = f;
        }
    }
}

// alineamiento multiple
void sol()
{
    // estrella final
    str alineamiento_estrella=star;
    for(auto alineamiento : alineamiento_multiple)
    {
        str tmp2;
        str alg1 = alineamiento.first;
        int i=0, j=0;
        while(i<alineamiento_estrella.size() || j<alg1.size())
        {
            if(alineamiento_estrella[i] == alg1[j] && alg1[j] != '_')
            {
                tmp2.push_back(alineamiento_estrella[i]);
                i++, j++;
            }
            else if(alineamiento_estrella[i] == '_')
            {
                tmp2.push_back('_');
                i++;
            }
            else if(alg1[j] == '_')
            {
                tmp2.push_back('_');
                j++;
            }
        }
        // cout<<tmp2<<"\n";
        alineamiento_estrella = tmp2;
    }
    alineamiento_sol.push_back(alineamiento_estrella);

    str alineamiento_cadena;
    // otros alineamientos
    for(auto alineamiento : alineamiento_multiple)
    {
        str alineamiento_cadena = alineamiento.first;
        str tmp2;
        // str alg1 = alineamiento.second;
        int i=0, j=0;
        while(i<alineamiento_cadena.size() || j<alineamiento_estrella.size())
        {
            if(alineamiento_cadena[i] == alineamiento_estrella[j])
            {
                tmp2.push_back(alineamiento.second[i]);
                i++, j++;
            }
            else if(alineamiento_cadena[i] == '_')
            {
                tmp2.push_back('_');
                i++;
            }
            else if(alineamiento_estrella[j] == '_')
            {
                tmp2.push_back('_');
                j++;
            }
        }
        // cout<<tmp2<<"\n";
        alineamiento_cadena = tmp2;
        alineamiento_sol.push_back(alineamiento_cadena);
    }
}

int main(int argc, char* argv[])
{
    int n,m;
    int match, mismatch, gap;
    int num_seq;

    num_seq = std::stoi(argv[1]);

    for (int i = 0; i < num_seq; i++)
    {
        cadenas.push_back(argv[i + 2]);
    }
    match = std::stoi(argv[num_seq + 2]);
    mismatch = std::stoi(argv[num_seq + 3]);
    gap = std::stoi(argv[num_seq + 4]);


    ofstream file("sa.txt");
    str cadena;

    n_cadenas = cadenas.size();
    matriz_scores.assign(n_cadenas, vi(n_cadenas+1,0));
    alineamientos.assign(n_cadenas, vpstrstr(n_cadenas));
    process();
    max_score();
    //print_scores(file);
    // file<<"\nAligments:\n";
    // print_alineamientos(file);
    star = cadenas[indx_max_score];
    file<< star;
    file<< global_max_score<< "\n";
    for(int c=0 ; c<n_cadenas ; ++c)
    {
        if(c!=indx_max_score)
        {
            alineamiento_multiple.push_back(alineamientos[indx_max_score][c]);
        }
    }
    sol();
    print_sol(file);
    file.close();
}