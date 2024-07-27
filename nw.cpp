// Alineamiento global

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#define str string
#define vi vector<int>
#define vvi vector<vi>
#define vvvi vector<vvi>
#define vpstrstr vector<pair<str,str>>
#define inf INT_MAX
using namespace std;

str String1, String2;
int n,m;
int match, missmatch, gap;

// print score matrix
void print(vvi& matrix, ofstream& file)
{
    for(int i=0 ; i<n+1 ; ++i)
    {
        for(int j=0 ; j<m+1 ; ++j)
            file<<matrix[i][j]<<" ";
        file<<"\n";
    }
    file<<"\n";
}

// print path matrix
void print(vvvi& matrix, ofstream& file)
{
    for(int i=0 ; i<n+1 ; ++i)
    {
        for(int j=0 ; j<m+1 ; ++j)
        {
            for(int k=0 ; k < matrix[i][j].size() ; ++k)
                file<<matrix[i][j][k];
            file<<" ";
        }
        file<<"\n";
    }
    file<<"\n";
}

// initialize arrays
void ini(vvi& matrix, vvvi& res)
{
    for(int i=1 ; i<n+1 ; ++i)
    {
        matrix[i][0] = matrix[i-1][0] + gap;
        res[i][0].push_back(1);    // up
    }
    for(int j=1 ; j<m+1 ; ++j)
    {
        matrix[0][j] = matrix[0][j-1] + gap;
        res[0][j].push_back(-1);    // left
    }
}

// get best score
void process(vvi& matrix, vvvi& res)
{
    // String2 < String1
    for(int f=1 ; f<n+1 ; ++f)  // String1
    {
        for(int c=1 ; c<m+1 ; ++c)  // String2
        {
            int value = 0;
            if(String1[f-1] == String2[c-1]) value = match;
            else value = missmatch;
            int diag, up, left, max_value;
            diag = matrix[f-1][c-1] + value;
            up = matrix[f-1][c] + gap;
            left = matrix[f][c-1] + gap;
            max_value = max(diag, max(up, left));
            // -1 left
            // 1 up
            // 0 diag
            if(diag == max_value) res[f][c].push_back(0);
            if(up == max_value) res[f][c].push_back(1);
            if(left == max_value) res[f][c].push_back(-1);
            matrix[f][c] = max_value;
        }
    }
}

// get path
void get_sol(int f, int c, str sol1, str sol2, vpstrstr& solutions, vvvi& res)
{
    if(f<0 || c<0)
        return;
    if(f==0 && c==0) 
    {
        // push solution
        solutions.push_back(make_pair(sol1, sol2));
        // cout<<solution<<'\n';
        return;
    }
    for(int i=0 ; i<res[f][c].size() ; ++i)
    {
        // cout<<"size: "<<res[f][c].size()<<" pos: "<<f<<" "<<c<<"\n";
        char caracter1 = ' ';
        char caracter2 = ' ';
        str solution1 = "";
        str solution2 = "";
        if(res[f][c][i] == -1) //left
        {
            caracter1 = String2[c-1];
            caracter2 = '-';
            solution1 = sol1 + caracter1;
            solution2 = sol2 + caracter2;
            get_sol(f, c-1, solution1, solution2, solutions, res);
        }
        else if(res[f][c][i] == 1)   // up
        {
            caracter1 = '-';
            caracter2 = String1[f-1];
            solution1 = sol1 + caracter1;
            solution2 = sol2 + caracter2;
            get_sol(f-1, c, solution1, solution2, solutions, res);
        }
        else if(res[f][c][i] == 0)  // diag
        {
            caracter1 = String2[c-1];
            caracter2 = String1[f-1];
            solution1 = sol1 + caracter1;
            solution2 = sol2 + caracter2;
            get_sol(f-1, c-1, solution1, solution2, solutions, res);
        }
    }
}

// verify score
void print_score(str& str1, str& str2, ofstream& file)
{
    if( str1.size() != str2.size())
    {
        file<<"strings must have same length\n";
        return;
    }
    int score = 0;
    for(int i=0 ; i<str1.size() ; i++)
    {
        if(str1[i]=='_' || str2[i]=='_')
            score-=2;
        else if(str1[i] == str2[i])
            score+=1;
        else
            score-=1;
    }
    file<<"score for "<<str1<<" and "<<str2<<": "<<score<<"\n";
}

str separator(str& A, str& B)
{
    str sep;
    for(int i=0 ; i<A.size() ; ++i)
    {
        if(A[i] == B[i])
            sep += '|';
        else
            sep += '-';
    }
    return sep;
}

int main(int argc, char* argv[])
{
    // freopen("input.txt", "r", stdin);
    String1 = argv[1];
    String2 = argv[2];
    match = stoi(argv[3]);
    missmatch = stoi(argv[4]);
    gap = stoi(argv[5]);
    cout<<match<<" "<<missmatch<<" "<<gap<<"\n";
    // cin>>String1>>String2;

    ofstream file("nw.txt");
    // string2 always bigger
    if (String2.size() > String1.size())
        swap(String1, String2);
    // file<<"Strings:\n";
    // file<<String1<<"\n";
    // file<<String2<<"\n\n";
    n = String1.size();
    m = String2.size();
    vvi matrix(n+1, vi(m+1, 0));
    vvvi res(n+1, vvi(m+1, vi(0)));
    ini(matrix, res);
    // print(matrix, file);

    // get score matrix
    process(matrix, res);

    // file<<"Final matrix\n";
    // print(matrix, file);
    // print(res, file);

    vpstrstr solutions;
    str sol1="";
    str sol2="";
    // get solutions
    get_sol(n, m, sol1, sol2, solutions, res);    
    file<<matrix[n][m]<<"\n";
    file<<solutions.size()<<"\n";
    for(auto pair_string : solutions)
    {
        reverse(pair_string.first.begin(), pair_string.first.end());
        reverse(pair_string.second.begin(), pair_string.second.end());
        str sep = separator(pair_string.first, pair_string.second);
        std::string cur;
        int start = 0;
        if (pair_string.first[0] == pair_string.second[0]) cur = "match";
        else if (pair_string.first[0] == '-' || pair_string.second[0] == '-') cur = "gap";
        else cur = "mismatch";
        file << "#\n";
        for (int i = 1; i < pair_string.first.length(); i++)
        {
            std::string type;
            if (pair_string.first[i] == pair_string.second[i]) type = "match";
            else if (pair_string.first[i] == '-' || pair_string.second[i] == '-') type = "gap";
            else type = "mismatch";
            if (type != cur)
            {
                file << pair_string.first.substr(start, i - start) << ' ' << pair_string.second.substr(start, i - start) << ' ' << sep.substr(start, i - start) << " " << cur <<'\n';
                start = i;
                cur = type;
            }
            if (i == pair_string.first.length() - 1) file << pair_string.first.substr(start, i - start) << ' ' << pair_string.second.substr(start, i - start) << ' ' << sep.substr(start, i - start) << " " << cur <<'\n';
        }
        // print_score(String1, String, file);
    }
    file.close();
}