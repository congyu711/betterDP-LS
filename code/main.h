#pragma once
#include <bits/stdc++.h>
using namespace std;

const int N=3000;

class prob
{
public:
    // graph
    int leftNum, rightNum;
    vector<int> l2r[N], r2l[N];
    // vector<int> lmovable(N),rmovable(N);
    int movable[2][N];
    // the current permutation of two sides. from [line number] to [index]
    int Permutation[2][N];
    // the optimal permutation
    int CPopt[2][N];
    // idx to line number
    int position[2][N];

    // M matrix in the paper
    // m[0][i][j]=k means: in the left part of the graph,
    // the number of crossings of i's and j's edges when i is above j (i,j are line number)
    // O(dlog(d))
    int m[2][N][N];

    // Delta matrix for insert and swap(in the paper)
    int Deltai[2][N][N], Deltas[2][N][N], Delta[2][N][N];

    // guarantee j<i
    // op[i][j] =0 -> do nothing
    //          1 -> insert i before j
    //          2 -> insert j after i
    //          3 -> swap i j
    int op[N][N];

    // current number of crossings in the bipartite graph.
    int currentsol, opt = 0x7fffffff;

    // tabu search table
    int TStable[2][N];

    // aux for random perturbation
    int shuffleaux[N];
    int auxsz;
    // dp
    int dp[N];
    tuple<int, int, int> ops[N][2];
    // functions
    int readGraph(string);
    void maintainPositions();
    int getcurrentsolution();
    void randPerturbation(int,double);
    void computeM(int, int);
    void computeDelta(int, int);
    void localsearch(bool USETS=false);
    void checkcp();
};

int runtimes = 0;
int timelimit = 5;
int tmp[N * N]; // aux array for mergesort
