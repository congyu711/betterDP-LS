#include "main.h"
using namespace std;
bool genTikz(string filename)
{
    int inc=3;
    ofstream outfile(filename);
    outfile << 
R"(\documentclass{minimal}
\usepackage{tikz}
\usetikzlibrary{positioning}
\usepackage[active,tightpage]{preview} 
\PreviewEnvironment{tikzpicture}
\setlength\PreviewBorder{4mm}
\begin{document}
\begin{tikzpicture})" <<endl;
    // \node[shape=circle,draw=black] (s) at (0,0) {s};
    for(int i=0;i<left_num;i++)
        outfile<<R"(\node[shape=circle,draw=black] ()"<<to_string(currentPermutation[0][i])+"l"<<") at (5,"<<inc*(i+1)<<") {"\
        <<to_string(currentPermutation[0][i])+"l"<<"};\n";
    for(int i=0;i<right_num;i++)
        outfile<<R"(\node[shape=circle,draw=black] ()"<<to_string(currentPermutation[1][i])+"r"<<") at (15,"<<inc*(i+1)<<") {"\
        <<to_string(currentPermutation[1][i])+"r"<<"};\n";
    for(int ii=0;ii<left_num;ii++)
    {
        int i=currentPermutation[0][ii];
        for(int j=0;j<l2r[i].size();j++)
        {
            // \path [-] (s) edge[thick,bend left=0] node[] {1} (a);
            outfile<<R"(\path [-] ()"<<to_string(i)+"l) edge[thick] node[] {} ("<<to_string(position[1][l2r[i][j]])+"r);\n";
        }
    }
    
    outfile<<
R"(\end{tikzpicture}
\end{document})"<<endl;
    return 0;
}
// int main()
// {
//     vector<vector<int>> g(4);
//     g[0].push_back(2);
//     g[0].push_back(2);
//     g[1].push_back(1);
//     g[1].push_back(2);
//     g[2].push_back(3);
//     g[3].push_back(2);
//     g[3].push_back(1);
//     genTikz(g, 4, 4, "a.tex");
// }