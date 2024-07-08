#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

class SmithWaterman 
{
public:

    std::string s;
    std::string t;

    int match;
    int mismatch;
    int gap;

    std::vector<std::vector<bool>> up_matrix;
    std::vector<std::vector<bool>> left_matrix;
    std::vector<std::vector<bool>> diag_matrix;

    std::vector<std::vector<int>> values_matrix;
public:
    SmithWaterman(std::string s, std::string t, int match_, int mismatch_, int gap_)
    {
        this->s = s;
        this->t = t;

        this->match = match_;
        this->mismatch = mismatch_;
        this->gap = gap_;

        up_matrix = std::vector<std::vector<bool>>(s.size() + 1, std::vector<bool>(t.size() + 1, false));
        left_matrix = std::vector<std::vector<bool>>(s.size() + 1, std::vector<bool>(t.size() + 1, false));
        diag_matrix = std::vector<std::vector<bool>>(s.size() + 1, std::vector<bool>(t.size() + 1, false));
        values_matrix = std::vector<std::vector<int>>(s.size() + 1, std::vector<int>(t.size() + 1, 0));

        for (int i = 0; i < s.size() + 1; i++)
        {
            values_matrix[i][0] = 0;
        }
        for (int i = 0; i < t.size() + 1; i++)
        {
            values_matrix[0][i] = 0;
        }

        int rows = s.size();
        int cols = t.size();

        for (int diagonal = 0; diagonal < rows + cols - 1; diagonal++) {
        
            int start_row = std::max(0, diagonal - cols + 1);
            int start_col = std::min(diagonal, cols - 1);

            int diag_len = std::min(rows - start_row - 1, start_col) + 1;

            for (int k = 0; k < diag_len; k++) {
                int i = start_row + 1 + k, j = start_col + 1 - k;
                int new_val;
                if (s[i - 1] == t[j - 1])
                {
                    new_val = std::max({values_matrix[i - 1][j] + gap, values_matrix[i][j - 1] + gap, values_matrix[i - 1][j - 1] + match, 0});
                }
                else
                {
                    new_val = std::max({values_matrix[i - 1][j] + gap, values_matrix[i][j - 1] + gap, values_matrix[i - 1][j - 1] + mismatch, 0});
                }
                if (new_val == values_matrix[i - 1][j] - 2 && new_val != 0)
                {
                    up_matrix[i][j] = true;
                }
                if (new_val == values_matrix[i][j - 1] - 2 && new_val != 0)
                {
                    left_matrix[i][j] = true;
                }
                if ((new_val == values_matrix[i - 1][j - 1] + match && s[i - 1] == t[j - 1]) || (new_val  == values_matrix[i - 1][j - 1] + mismatch && s[i - 1] != t[j - 1]) && new_val != 0)
                {
                    diag_matrix[i][j] = true;
                }

                values_matrix[i][j] = new_val;
                
            }
        }

    }

    std::vector<std::vector<int>> get_score_matrix()
    {
        return this->values_matrix;
    }

    int get_score()
    {
        return this->values_matrix[s.size()][t.size()];
    }


    unsigned long long optimal_allignments_num()
    {
        std::vector<std::vector<unsigned long long>> dp(s.size() + 1, std::vector<unsigned long long>(t.size() + 1, 0));
        int rows = s.size();
        int cols = t.size();

        for (int diagonal = 0; diagonal < rows + cols - 1; diagonal++)
        {
        
            int start_row = std::max(0, diagonal - cols + 1);
            int start_col = std::min(diagonal, cols - 1);

            int diag_len = std::min(rows - start_row - 1, start_col) + 1;

            for (int i = 0; i < s.size() + 1; i++)
            {
                dp[i][0] = 1;
            }
            for (int i = 0; i < t.size() + 1; i++)
            {
                dp[0][i] = 1;
            }

            for (int k = 0; k < diag_len; k++) 
            {
                int i = start_row + 1 + k, j = start_col + 1 - k;
                if (up_matrix[i][j])
                {
                    dp[i][j] += dp[i - 1][j];
                }
                if (left_matrix[i][j])
                {
                    dp[i][j] += dp[i][j - 1];
                }
                if (diag_matrix[i][j])
                {
                    dp[i][j] += dp[i - 1][j - 1];
                }
            }
        }

        return dp[rows][cols];
    }

    std::pair<std::string, std::string> get_one_alignment()
    {
        int i = this->s.size();
        int j = this->t.size();

        std::string alignment_s = "";
        std::string alignment_t = "";

        while (i > 0 || j > 0)
        {
            if (i == 0)
            {
                while (j > 0)
                {
                    alignment_s.insert(0, 1, '_');
                    alignment_t.insert(0, 1, this->t[j - 1]);
                    j--;
                }
                break;
            }
            else if (j == 0)
            {
                while (i > 0)
                {
                    alignment_s.insert(0, 1, this->s[i - 1]);
                    alignment_t.insert(0, 1, '_');
                    i--;
                }
                break;
            }
            else if (diag_matrix[i][j])
            { 
                alignment_s.insert(0, 1, this->s[i - 1]);
                alignment_t.insert(0, 1, this->t[j - 1]);
                i--;
                j--;
            }
            else if (up_matrix[i][j])
            {
                alignment_s.insert(0, 1, this->s[i - 1]);
                alignment_t.insert(0, 1, '_');
                i--;
            }
            else
            {
                alignment_s.insert(0, 1, '_');
                alignment_t.insert(0, 1, this->t[j - 1]);
                j--;
            }
        }

        return std::make_pair(alignment_s, alignment_t);
    }

    bool is_valid_alignment(int x, int y)
    {
        int xx = x;
        int yy = y;
        for (int i = this->values_matrix[x][y]; i >= 0; i--)
        {
            if (xx < 0 || yy < 0) return false;
            if (i != this->values_matrix[xx][yy])
            {
                return false;
            }
            if (i != 0 && this->s[xx - 1] != this->t[yy - 1])
            {
                return false;
            }
            xx--;
            yy--;
        }
        return true;
    }

    std::vector<std::string> get_allignments()
    {
        int max_val = 0;
        
        std::vector<int> xs;
        std::vector<int> ys;

        for (int i = 0; i < this->values_matrix.size(); i++)
        {
            for (int j = 0; j < this->values_matrix[i].size(); j++)
            {
                if (max_val < this->values_matrix[i][j]) 
                {
                    if (is_valid_alignment(i, j))
                    {
                        max_val = this->values_matrix[i][j];
                        xs = {i};
                        ys = {j};
                    }
                }
                else if (max_val == this->values_matrix[i][j])
                {
                    if (is_valid_alignment(i, j))
                    {
                        xs.push_back(i);
                        ys.push_back(j);
                    }
                }
            }
        }
        std::vector<std::string> alignments;
        for (int i = 0; i < xs.size(); i++)
        {
            alignments.push_back(this->s.substr(xs[i] - max_val, (xs[i] - 1) - (xs[i] - max_val) + 1));
        }
        return alignments;
    }

    int get_max_score()
    {
        int max_val = 0;
        for (int i = 0; i < this->values_matrix.size(); i++)
        {
            for (int j = 0; j < this->values_matrix[i].size(); j++)
            {
                if (max_val < this->values_matrix[i][j])
                {
                    if (is_valid_alignment(i,j)) max_val = this->values_matrix[i][j];
                }
            }
        }
        return max_val;
    }

    std::pair<std::pair<int,int>,std::pair<int,int>> get_positions()
    {
        int max_val = 0;
        int x, y;
        for (int i = 0; i < this->values_matrix.size(); i++)
        {
            for (int j = 0; j < this->values_matrix[i].size(); j++)
            {
                if (max_val < this->values_matrix[i][j]) 
                {
                    max_val = this->values_matrix[i][j];
                    x = i;
                    y = j;
                }    
            }
        }

        return std::make_pair(std::make_pair(x - max_val, x - 1), std::make_pair(y - max_val, y - 1));
    }

    std::string get_alignment()
    {
        int max_val = 0;
        int x, y;
        for (int i = 0; i < this->values_matrix.size(); i++)
        {
            for (int j = 0; j < this->values_matrix[i].size(); j++)
            {
                if (max_val < this->values_matrix[i][j]) 
                {
                    max_val = this->values_matrix[i][j];
                    x = i;
                    y = j;
                }    
            }
        }

        return this->s.substr(x - max_val, (x - 1) - (x - max_val) + 1);
    }

    int get_num_allignments()
    {
        int max_val = 0;
        int count = 0;
        for (int i = 0; i < this->values_matrix.size(); i++)
        {
            for (int j = 0; j < this->values_matrix[i].size(); j++)
            {
                if (max_val < this->values_matrix[i][j] && is_valid_alignment(i, j)) 
                {
                    max_val = this->values_matrix[i][j];
                    count = 1;
                }
                else if (max_val == this->values_matrix[i][j] && is_valid_alignment(i, j)) count++;
            }
        }

        return count;
    }

    



};

std::string get_random_dna_sequence(int length)
{
    std::string alphabet = "ACGT";
    std::string random_seq = "";
    for (int i = 0; i < length; i++)
    {
        int p = rand() % 4;
        random_seq.push_back(alphabet[p]);
    }
    return random_seq;
}

std::pair<std::string, std::string> get_sequences_from_file(std::string filename)
{
    std::ifstream file(filename.c_str());
    std::string s, t;
    file >> s;
    file >> t;

    return std::make_pair(s, t); 
}

int main(int argc, char* argv[])
{
    // freopen("input.txt", "r", stdin);
    std::string String1, String2;
    int n,m;
    int match, missmatch, gap;

    String1 = argv[1];
    String2 = argv[2];
    match = std::stoi(argv[3]);
    missmatch = std::stoi(argv[4]);
    gap = std::stoi(argv[5]);
    std::cout<<match<<" "<<missmatch<<" "<<gap<<"\n";
    // cin>>String1>>String2;

    std::ofstream file("sw.txt");
    // string2 always bigger
    if (String2.size() > String1.size())
        swap(String1, String2);
    // file<<"Strings:\n";
    // file<<String1<<"\n";
    // file<<String2<<"\n\n";
    SmithWaterman alineacion(String1, String2, match, missmatch, gap);

    std::vector<std::vector<int>> score_matrix = alineacion.get_score_matrix();

    file << alineacion.get_max_score() << '\n';
    
    file << alineacion.get_num_allignments() << '\n';

    auto alineamientos = alineacion.get_allignments();

    for (const auto & a: alineamientos)
    {
        file << a << '\n';
    }
    file.close();

    return 0;
}