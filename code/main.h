#include <bits/stdc++.h>
#include "counttime.cpp"
using namespace std;

const int N=1000;
mt19937 gen(time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count());
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
    int cnt;        // TS step counter

    // aux for random perturbation
    int shuffleaux[N];
    int auxsz;
    // dp
    int dp[N];
    tuple<int, int, int> ops[N][2];

    double optTime;
    // functions
    int readGraph(string);
    void maintainPositions();
    int getcurrentsolution();
    void randPerturbation(int,double);
    void computeM(int, int);
    void computeDelta(int, int);
    void localsearch(bool, bool);
    void checkcp();
    void tabu(int, int, int);
    void Delta_dp(int,bool);
};
class gene
{
public:
    int Permutation[2][N];
};

int timelimit = 700;
int tmp[N*N]; // aux array for mergesort

const int groupsize=50;
vector<gene> genes(groupsize);
priority_queue<pair<int,int>> pq;   // opt, index
unordered_set<int> groupidxs;
int opt=INT32_MAX;
double optTime;