// main.cpp
#include "main.h"
#include "counttime.cpp"
// #include "segment_tree.cpp"
// #include "genTikz.cpp"

mt19937 gen(time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count());
int prob::readgraph(string filepath)
{
    int a;
    ifstream fin(filepath);
    fin >> a;
    fin >> left_num >> right_num;
    for (int i = 0; i < left_num; i++)
    {
        int idx, tmp;
        fin >> tmp;
        //////////////////////////////
        // ignore movable bit
        // lmovable[i]=0;
        //////////////////////////////

        fin >> idx;
        movable[0][idx] = tmp;
        movable[0][idx] ^= 1;
        // currentPermutation[0][i] = i; not idx
        currentPermutation[0][i] = i;
        while (fin.get() != '\n')
        {
            fin >> tmp;
            // this means
            // 1) vertex idx has an edge to vertex tmp-left_num; or
            // 2) vertex idx has an edge to that vertex of the tmp-left_num th line in the input file.
            l2r[idx].push_back(tmp - left_num);
        }
    }
    for (int i = 0; i < right_num; i++)
    {
        int idx, tmp;
        fin >> tmp >> idx;
        if (idx >= right_num)
            idx -= left_num;
        movable[1][idx] = tmp;
        movable[1][idx] ^= 1;

        currentPermutation[1][i] = idx;
        //////////////////////////////
        // ignore movable bit
        // rmovable[i]=0;
        //////////////////////////////
    }

    ///////////////////////////////////////
    // here use 2)
    // because NewGraph.cpp use 2)      :(
    for (int i = 0; i < left_num; i++)
    {
        for (auto &e : l2r[i])
            e = currentPermutation[1][e];
    }
    ///////////////////////////////////////

    ///////////////////////////////////////
    for (int i = 0; i < right_num; i++)
        currentPermutation[1][i] = i;
    ///////////////////////////////////////

    // compute r2l graph
    for (int i = 0; i < left_num; i++)
        for (auto &e : l2r[i])
            r2l[e].push_back(i);
    fin.close();
    return 0;
}
// just maintain positions array according to the currentpermutation array
void prob::maintainpositions()
{
    // maintain LRpositions array
    for (int i = 0; i < left_num; i++)
        position[0][currentPermutation[0][i]] = i;
    for (int i = 0; i < right_num; i++)
        position[1][currentPermutation[1][i]] = i;
}
// mergesort: compute the number of inversion
int tmp[N * N];
int mergesort(vector<int> &a, int l, int r)
{
    int res = 0;
    if (r - l <= 1)
        return res;
    int mid = (l + r) >> 1;
    res += mergesort(a, l, mid);
    res += mergesort(a, mid, r);
    int p = l, q = mid, s = l;
    while (s < r)
    {
        if (p >= mid || (q < r && a[p] > a[q]))
        {
            tmp[s++] = a[q++];
            res += (mid - p);
        }
        else
            tmp[s++] = a[p++];
    }
    for (int i = l; i < r; i++)
        a[i] = tmp[i];
    return res;
}
// compute total number of edge crossings
int prob::getcurrentsolution()
{
    maintainpositions();
    vector<int> aux;
    for (int i = 0; i < left_num; i++)
    {
        vector<int> auxx;
        for (auto e : l2r[currentPermutation[0][i]])
            auxx.push_back(position[1][e]);
        sort(auxx.begin(), auxx.end());
        for (auto e : auxx)
            aux.push_back(e);
    }
    return mergesort(aux, 0, aux.size());
}


/**
 * @brief give a random permutation of currentpermutation array.
 * 
 * @param lr 0 for the left side, 1 for the right side.
 * @param p p percent of movable vertices will be rearranged. default is 100%
 */
void prob::random_perturbation(int lr,double p=1.0)
{
    srand(gen());
    ///////////////////////////////////////
    int num = (lr ? right_num : left_num);
    ///////////////////////////////////////
    vector<bool> shuffleflag(num,false);
    auxsz = 0;
    for (int i = 0; i < num; i++)
    {
        // can not move
        if (movable[lr][currentPermutation[lr][i]] == 0 || gen()>p*0xffffffff)
        {
            shuffleaux[auxsz++] = currentPermutation[lr][i];
            shuffleflag[currentPermutation[lr][i]]=true;
        }
    }
    random_shuffle(currentPermutation[lr], currentPermutation[lr] + num);
    auxsz = 0;
    for (int i = 0; i < num; i++)
    {
        if (shuffleflag[currentPermutation[lr][i]])
            currentPermutation[lr][i] = shuffleaux[auxsz++];
    }
}
void prob::computeM(int lr, int num)
{
    ///////////////////////////////////////
    auto *g = (lr ? r2l : l2r);
    ///////////////////////////////////////

    for (int i = 0; i < num; i++)
    {
        for (int j = i + 1; j < num; j++)
        {
            m[lr][i][j]=0;
            vector<int> vi, vj;
            for (auto e : g[i])
                vi.push_back(position[lr ^ 1][e]);
            for (auto e : g[j])
                vj.push_back(position[lr ^ 1][e]);
            sort(vi.begin(), vi.end());
            sort(vj.begin(), vj.end());
            // compute the number of inversion
            auto sz = vi.size() + vj.size();
            vector<int> vij(sz), vji(sz);
            for (int i = 0; i < vi.size(); i++)
            {
                vij[i] = vi[i];
                vji[i + vj.size()] = vi[i];
            }
            for (int i = 0; i < vj.size(); i++)
            {
                vji[i] = vj[i];
                vij[i + vi.size()] = vj[i];
            }
            if(movable[lr][i]||movable[lr][j])
            {
                m[lr][i][j] = mergesort(vij, 0, sz);
                m[lr][j][i] = mergesort(vji, 0, sz);
            }
            else if(position[lr][i]<position[lr][j])
                m[lr][i][j] = mergesort(vij,0,sz);
            else if(position[lr][i]>position[lr][j])
                m[lr][j][i] = mergesort(vji,0,sz);
        }
    }
}
void prob::computeDelta(int lr, int num)
{

    for (int i = num - 1; i >= 0; i--)
    {
        if (movable[lr][currentPermutation[lr][i]] == 0)
        {    
            for(int j=0;j<num;j++)  Deltai[lr][i][j]=0;
            continue;
        }
        for (int j = i - 1; j >= 0; j--)
        {
            int ii=currentPermutation[lr][i],jj=currentPermutation[lr][j];
            if (j == i - 1)
                Deltai[lr][i][j] = m[lr][ii][jj] - m[lr][jj][ii];
            else
                Deltai[lr][i][j] = Deltai[lr][i][j + 1] + m[lr][ii][jj] - m[lr][jj][ii];
        }
        for (int j = i + 1; j < num; j++)
        {
            int ii=currentPermutation[lr][i],jj=currentPermutation[lr][j];
            if (j == i + 1)
                Deltai[lr][i][j] = m[lr][jj][ii] - m[lr][ii][jj];
            else
                Deltai[lr][i][j] = Deltai[lr][i][j - 1] - m[lr][ii][jj] + m[lr][jj][ii];
        }
    }
    // compute Delta_swap matrix
    // memset(Deltas, 0, sizeof(Deltas));
    for (int i = 0; i < num; i++)
    {
        if (movable[lr][currentPermutation[lr][i]] == 0)
        {
            for(int j=0;j<num;j++)  Deltas[lr][i][j]=0;
            continue;
        }
        // j<=i meaningless
        // j=i+1 can be handled by inserting.
        for (int j = i + 1; j < num; j++)
        {
            if (movable[lr][currentPermutation[lr][j]] == 0)
            {    
                Deltas[lr][i][j]=0;
                continue;
            }
            Deltas[lr][i][j] = Deltai[lr][j][i] + Deltai[lr][i][j - 1];
            Deltas[lr][j][i] = Deltai[lr][j][i] + Deltai[lr][i][j - 1];
        }
    }
}
/**
 * @brief dp2 is a perturbation.
 * 
 * @return true if currentpermutation is changed;
 * @return false otherwise.
 */
bool prob::dp2()
{
    int lr = 0;
    bool f = 0;
    int cnt = 100;
    int presol = 0;
    maintainpositions();
    computeM(0,left_num);
    computeM(1,right_num);
    currentsol=getcurrentsolution();
    int n = min(left_num, right_num) * 2;
    for (int k = 4; k < n; k+=7)
    {
        lr = k & 1;
        /////////////////////////////////////
        auto *g = (lr ? r2l : l2r);
        int num = k / 2;
        /////////////////////////////////////
        maintainpositions();
        computeDelta(lr, num);
        // here delta matrix should be non-positive.
        for (int i = 0; i < num; i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (movable[lr][currentPermutation[lr][j]] && movable[lr][currentPermutation[lr][i]])
                {
                    Delta[lr][i][j] = min(min(Deltai[lr][i][j], Deltai[lr][j][i]), Deltas[lr][i][j]);
                    if (Delta[lr][i][j] >= 0)
                        op[i][j] = 0;
                    else if (Delta[lr][i][j] == Deltai[lr][i][j])
                        op[i][j] = 1;
                    else if (Delta[lr][i][j] == Deltai[lr][j][i])
                        op[i][j] = 2;
                    else
                        op[i][j] = 3;
                }
                else if (movable[lr][currentPermutation[lr][i]])
                {
                    Delta[lr][i][j] = Deltai[lr][i][j];
                    if (Delta[lr][i][j] < 0)
                        op[i][j] = 1;
                    else
                        op[i][j] = 0, Delta[lr][i][j] = 0;
                }
                else if (movable[lr][currentPermutation[lr][j]])
                {
                    Delta[lr][i][j] = Deltai[lr][j][i];
                    if (Delta[lr][i][j] < 0)
                        op[i][j] = 2;
                    else
                        op[i][j] = 0, Delta[lr][i][j] = 0;
                }
                else
                {
                    op[i][j] = 0;
                    Delta[lr][i][j] = 0;
                }
            }
        }
        // dp
        for (int i = 0; i < num; i++)
        {
            dp[i] = 0;
            if (i > 1)
            {
                int minn = 0, minidx = -1;
                for (int j = 1; j <= i; j++)
                {
                    if (minn > Delta[lr][i][j] + dp[j - 1])
                    {
                        minn = Delta[lr][i][j] + dp[j - 1];
                        minidx = j - 1;
                    }
                }
                if (minidx == -1)
                {
                    dp[i] = 0;
                    ops[i][0] = make_tuple(0, 0, 0);
                    ops[i][1] = make_tuple(-1, i - 1, 0);
                }
                else
                {
                    dp[i] = minn;
                    int j = minidx;
                    ops[i][0] = make_tuple(op[i][j + 1], i, j + 1);
                    ops[i][1] = make_tuple(-1, j, 0);
                }
            }
            else if (i == 1)
            {
                dp[i] = Delta[lr][0][1];
                if (op[0][1] != 0)
                {
                    ops[i][0] = make_tuple(op[0][1], 1, 0);
                    ops[i][1] = make_tuple(0, 0, 0);
                }
            }
        }

        // adjust Permutations.
        // O(n) because there is no intersaction between operations.
        int pos = num - 1;
        int ppos = -1;
        while (pos != ppos)
        {
            ppos = pos;

            int i, j, l, t;
            tie(l, i, j) = ops[pos][0];

            // int count[(lr?left_num:right_num)]={0};
            vector<int> count((lr?left_num:right_num),0);
            set<int> idxs;
            switch (l)
            {
            case 1:
                // insert i before j
                // maintain M matrix
                f=1;
                for(int ii=j;ii<i;ii++)
                {
                    for(auto ridx:g[currentPermutation[lr][ii]])
                    {
                        count[ridx]++;
                        idxs.insert(ridx);
                    }
                }
                for(auto ridxa:g[currentPermutation[lr][i]])
                {
                    for(auto ridxb:idxs)
                    {
                        m[lr^1][ridxb][ridxa]+=count[ridxb];
                        m[lr^1][ridxa][ridxb]-=count[ridxb];
                    }
                }
                tabutable[lr][currentPermutation[lr][i]] = cnt;
                t = currentPermutation[lr][i];
                for (int ii = i; ii >= j + 1; ii--)
                {
                    currentPermutation[lr][ii] = currentPermutation[lr][ii - 1];
                }
                currentPermutation[lr][j] = t;
                break;

            case 2:
                // 2 -> insert j after i
                // maintain M
                f=1;
                for(int ii=j+1;ii<=i;ii++)
                {
                    for(auto ridx:g[currentPermutation[lr][ii]])
                    {
                        count[ridx]++;
                        idxs.insert(ridx);
                    }
                }
                for(auto ridxa:g[currentPermutation[lr][j]])
                {
                    for(auto ridxb:idxs)
                    {
                        m[lr^1][ridxb][ridxa]-=count[ridxb];
                        m[lr^1][ridxa][ridxb]+=count[ridxb];
                    }
                }
                tabutable[lr][currentPermutation[lr][j]] = cnt;
                t = currentPermutation[lr][j];
                for (int ii = j; ii <= i - 1; ii++)
                {
                    currentPermutation[lr][ii] = currentPermutation[lr][ii + 1];
                }
                currentPermutation[lr][i] = t;
                break;

            case 3:
            {
                // 3 -> swap
                //  maintain M
                // insert i before j
                f=1;
                for(int ii=j;ii<i;ii++)
                {
                    for(auto ridx:g[currentPermutation[lr][ii]])
                    {
                        count[ridx]++;
                        idxs.insert(ridx);
                    }
                }
                for(auto ridxa:g[currentPermutation[lr][i]])
                {
                    for(auto ridxb:idxs)
                    {
                        m[lr^1][ridxb][ridxa]+=count[ridxb];
                        m[lr^1][ridxa][ridxb]-=count[ridxb];
                    }
                }
                vector<int> count1((lr?left_num:right_num),0);
                set<int> idxs1;
                // insert j after i-1
                for(int ii=j+1;ii<=i-1;ii++)
                {
                    for(auto ridx:g[currentPermutation[lr][ii]])
                    {
                        count1[ridx]++;
                        idxs1.insert(ridx);
                    }
                }
                for(auto ridxa:g[currentPermutation[lr][j]])
                {
                    for(auto ridxb:idxs1)
                    {
                        m[lr^1][ridxb][ridxa]-=count1[ridxb];
                        m[lr^1][ridxa][ridxb]+=count1[ridxb];
                    }
                }
                tabutable[lr][currentPermutation[lr][i]] = cnt;
                tabutable[lr][currentPermutation[lr][j]] = cnt;
                swap(currentPermutation[lr][i], currentPermutation[lr][j]);

                break;
            }

            default:
                break;
            }
            tie(l, i, j) = ops[pos][1];

            if (l == (-1))
                pos = i;
        }
        maintainpositions();
        currentsol+=dp[num-1];
        if (lr == 1 && presol == currentsol)
        {
            // new opt
            if(opt>currentsol)
            {
                opt=currentsol;
            }
        }
        ed = system_clock::now();
        if (counttime(st, ed) > timelimit)
            break;
    }
    return f;
}
int dp3aux[2][N];
/**
 * @brief dp3 is a perturbation.
 * 
 * @return true if currentpermutation is changed;
 * @return false otherwise.
 */
bool prob::dp3()
{
    static prob* p=new prob;
    // block size and new `left_num` & `right_num`
    int blocksz[2];
    blocksz[0]=sqrt(left_num);
    blocksz[1]=sqrt(right_num);
    p->left_num=ceil((double)left_num/blocksz[0]);
    p->right_num=ceil((double)right_num/blocksz[1]);
    memcpy(dp3aux[0],currentPermutation[0],left_num*4);
    memcpy(dp3aux[1],currentPermutation[1],right_num*4);
    // movable
    auto getnewidx=[&](int x,int lr)->int{return ceil((double)(x+1)/blocksz[lr])-1;};
    for(int i=0;i<p->left_num;i++)  p->movable[0][i]=1;
    for(int i=0;i<p->right_num;i++)  p->movable[1][i]=1;
    for(int i=0;i<left_num;i++)
    {
        if(movable[0][i]==0)
        {
            p->movable[0][getnewidx(i,0)]=0;
            i=(getnewidx(i,0)+1)*blocksz[0];
        }
    }
    for(int i=0;i<right_num;i++)
    {
        if(movable[1][i]==0)
        {
            p->movable[1][getnewidx(i,1)]=0;
            i=(getnewidx(i,1)+1)*blocksz[1];
        }
    }

    // graph
    // do not consider multi edge
    for(int i=0;i<left_num;i++)
    {
        unordered_set<int> toR;
        for(auto &e:l2r[i]) toR.insert(e);
        for(auto &e:toR)    p->l2r[getnewidx(i,0)].push_back(getnewidx(e,1));
    }
    for (int i = 0; i < p->left_num; i++)
        for (auto &e : p->l2r[i])
            p->r2l[e].push_back(i);
    
    // permutations
    for(int i=0;i<p->left_num;i++)  p->currentPermutation[0][i]=i;
    for(int i=0;i<p->right_num;i++)  p->currentPermutation[1][i]=i;
    p->maintainpositions();

    // local search (avoid recursive function call)
    p->currentsol=p->getcurrentsolution();
    p->localsearch();
    p->maintainpositions();
    // get new permutation.
    int lidx=0,ridx=0;
    for(int i=0;i<p->left_num;i++)
    {
        int c=p->currentPermutation[0][i];
        for(int j=0;j<(c!=p->left_num-1?blocksz[0]:left_num%blocksz[0]);j++)
            currentPermutation[0][lidx++]=dp3aux[0][c*blocksz[0]+j];
    }
    for(int i=0;i<p->right_num;i++)
    {
        int c=p->currentPermutation[1][i];
        for(int j=0;j<(c!=p->right_num-1?blocksz[1]:right_num%blocksz[1]);j++)
            currentPermutation[1][ridx++]=dp3aux[1][c*blocksz[1]+j];
    }
    maintainpositions();
    computeM(0,left_num);
    computeM(1,right_num);
    currentsol=getcurrentsolution();
    
    delete p;
    bool f=0;
    for(int i=0;i<10;i++)
    {
        int LR=gen()%2;
        int k=gen()%min(left_num,right_num);
        if(currentPermutation[LR][k]!=dp3aux[LR][k])    f=1;
    }
    return f;
}
/**
 * @brief do the local search
 * 
 * @param USEDP_ set to true to use dp3 as a perturbation.
 * false to use random perturbation.
 * @param USETABU true to use tabu search.
 */
void prob::localsearch(bool USEDP_,bool USETABU)
{
    int lr = 0;
    int cnt = 100;
    int presol = 0;
    ops[0][0]=ops[0][1]=make_tuple(-1,0,0);
    auto tabu=[&](int i,int j)->void{
        int randtmp = 7 + gen() % 5;
        if (op[i][j] == 1 && cnt - tabutable[lr][currentPermutation[lr][i]] < randtmp)
        {
            op[i][j] = 0;
            Delta[lr][i][j] = 0;
        }
        else if (op[i][j] == 2 && cnt - tabutable[lr][currentPermutation[lr][j]] < randtmp)
        {
            op[i][j] = 0;
            Delta[lr][i][j] = 0;
        }
        else if (op[i][j] == 3 && (cnt - tabutable[lr][currentPermutation[lr][i]] < randtmp 
        || cnt - tabutable[lr][currentPermutation[lr][j]] < randtmp))
        {
            op[i][j] = 0;
            Delta[lr][i][j] = 0;
        }
    };
    auto Delta_dp=[&]()->void{
        /////////////////////////////////////
        int num = (lr ? right_num : left_num);
        auto *g = (lr ? r2l : l2r);
        /////////////////////////////////////
        maintainpositions();
        // computeM(lr, num);
        computeDelta(lr, num);
        // here delta matrix should be non-positive.
        for (int i = 0; i < num; i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (movable[lr][currentPermutation[lr][j]] && movable[lr][currentPermutation[lr][i]])
                {
                    Delta[lr][i][j] = min(min(Deltai[lr][i][j], Deltai[lr][j][i]), Deltas[lr][i][j]);
                    if (Delta[lr][i][j] >= 0)
                        op[i][j] = 0;
                    else if (Delta[lr][i][j] == Deltai[lr][i][j])
                        op[i][j] = 1;
                    else if (Delta[lr][i][j] == Deltai[lr][j][i])
                        op[i][j] = 2;
                    else
                        op[i][j] = 3;
                }
                else if (movable[lr][currentPermutation[lr][i]])
                {
                    Delta[lr][i][j] = Deltai[lr][i][j];
                    if (Delta[lr][i][j] < 0)
                        op[i][j] = 1;
                    else
                        op[i][j] = 0, Delta[lr][i][j] = 0;
                }
                else if (movable[lr][currentPermutation[lr][j]])
                {
                    Delta[lr][i][j] = Deltai[lr][j][i];
                    if (Delta[lr][i][j] < 0)
                        op[i][j] = 2;
                    else
                        op[i][j] = 0, Delta[lr][i][j] = 0;
                }
                else
                {
                    op[i][j] = 0;
                    Delta[lr][i][j] = 0;
                }
                if(USETABU) tabu(i,j);
            }
        }
        // dp
        for (int i = 0; i < num; i++)
        {
            dp[i] = 0;
            if (i > 1)
            {
                int minn = 0, minidx = -2;
                if(Delta[lr][i][0]<0)
                {
                    minn=Delta[lr][i][0];
                    minidx=-1;
                }
                for (int j = 1; j <= i; j++)
                {
                    if (minn > Delta[lr][i][j] + dp[j - 1])
                    {
                        minn = Delta[lr][i][j] + dp[j - 1];
                        minidx = j - 1;
                    }
                }
                if (minidx == -2)
                {
                    dp[i] = 0;
                    ops[i][0] = make_tuple(0, 0, 0);
                    ops[i][1] = make_tuple(-1, 0, 0);
                }
                else
                {
                    dp[i] = minn;
                    int j = minidx;
                    ops[i][0] = make_tuple(op[i][j + 1], i, j + 1);
                    if(minn!=-1)
                        ops[i][1] = make_tuple(-1, j, 0);
                    else    ops[i][1]=make_tuple(-1,0,0);
                }
            }
            else if (i == 1)
            {
                dp[i] = Delta[lr][0][1];
                if (op[0][1] != 0)
                {
                    ops[i][0] = make_tuple(op[0][1], 1, 0);
                    ops[i][1] = make_tuple(0, 0, 0);
                }
            }
        }

    };
    maintainpositions();
    computeM(0,left_num);
    computeM(1,right_num);
    while (1)
    {
        cnt++;
        /////////////////////////////////////
        int num = (lr ? right_num : left_num);
        auto *g = (lr ? r2l : l2r);
        /////////////////////////////////////
        Delta_dp();
        // adjust Permutations.
        // O(n) because there is no intersaction between operations.
        int pos = num - 1;
        int ppos = -1;
        while (pos != ppos)
        {
            ppos = pos;

            int i, j, l, t;
            tie(l, i, j) = ops[pos][0];

            // int count[(lr?left_num:right_num)]={0};
            vector<int> count((lr?left_num:right_num),0);
            set<int> idxs;
            switch (l)
            {
            case 1:
                // insert i before j
                // maintain M matrix
                for(int ii=j;ii<i;ii++)
                {
                    // cout<<"idxsleft: "<<position[lr][ii]<<endl;
                    for(auto ridx:g[currentPermutation[lr][ii]])
                    {
                        count[ridx]++;
                        idxs.insert(ridx);
                    }
                }
                for(auto ridxa:g[currentPermutation[lr][i]])
                {
                    for(auto ridxb:idxs)
                    {
                        m[lr^1][ridxb][ridxa]+=count[ridxb];
                        m[lr^1][ridxa][ridxb]-=count[ridxb];
                    }
                }
                tabutable[lr][currentPermutation[lr][i]] = cnt;
                t = currentPermutation[lr][i];
                for (int ii = i; ii >= j + 1; ii--)
                {
                    currentPermutation[lr][ii] = currentPermutation[lr][ii - 1];
                }
                currentPermutation[lr][j] = t;
                // maintainpositions();
                // computeM(lr^1,(lr?left_num:right_num));
                break;

            case 2:
                // 2 -> insert j after i
                // maintain M
                for(int ii=j+1;ii<=i;ii++)
                {
                    for(auto ridx:g[currentPermutation[lr][ii]])
                    {
                        count[ridx]++;
                        idxs.insert(ridx);
                    }
                }
                for(auto ridxa:g[currentPermutation[lr][j]])
                {
                    for(auto ridxb:idxs)
                    {
                        m[lr^1][ridxb][ridxa]-=count[ridxb];
                        m[lr^1][ridxa][ridxb]+=count[ridxb];
                    }
                }
                tabutable[lr][currentPermutation[lr][j]] = cnt;
                t = currentPermutation[lr][j];
                for (int ii = j; ii <= i - 1; ii++)
                {
                    currentPermutation[lr][ii] = currentPermutation[lr][ii + 1];
                }
                currentPermutation[lr][i] = t;
                break;

            case 3:
            {
                // 3 -> swap
                //  maintain M
                // insert i before j
                for(int ii=j;ii<i;ii++)
                {
                    for(auto ridx:g[currentPermutation[lr][ii]])
                    {
                        count[ridx]++;
                        idxs.insert(ridx);
                    }
                }
                for(auto ridxa:g[currentPermutation[lr][i]])
                {
                    for(auto ridxb:idxs)
                    {
                        m[lr^1][ridxb][ridxa]+=count[ridxb];
                        m[lr^1][ridxa][ridxb]-=count[ridxb];
                    }
                }
                vector<int> count1((lr?left_num:right_num),0);
                set<int> idxs1;
                // insert j after i-1
                for(int ii=j+1;ii<=i-1;ii++)
                {
                    for(auto ridx:g[currentPermutation[lr][ii]])
                    {
                        count1[ridx]++;
                        idxs1.insert(ridx);
                    }
                }
                for(auto ridxa:g[currentPermutation[lr][j]])
                {
                    for(auto ridxb:idxs1)
                    {
                        m[lr^1][ridxb][ridxa]-=count1[ridxb];
                        m[lr^1][ridxa][ridxb]+=count1[ridxb];
                    }
                }
                tabutable[lr][currentPermutation[lr][i]] = cnt;
                tabutable[lr][currentPermutation[lr][j]] = cnt;
                swap(currentPermutation[lr][i], currentPermutation[lr][j]);
                break;
            }

            default:
                break;
            }
            tie(l, i, j) = ops[pos][1];

            if (l == (-1))
                pos = i;
        }
        maintainpositions();
        currentsol+=dp[num-1];
        // can not improve the current solution.
        if (lr == 1 && presol == currentsol)
        {
            // if TS is used, reset cnt and do it again.
            if(USETABU)
            {
                cnt+=5000;
                Delta_dp();
            }
            if(!USETABU||dp[num-1]==0)
            {
                // new opt
                if (opt > currentsol)
                {
                    opt = currentsol;
                    memcpy(CPopt[0], currentPermutation[0], 4 * left_num);
                    memcpy(CPopt[1], currentPermutation[1], 4 * right_num);
                }
                // same solution
                else if (opt<currentsol||(opt == currentsol && gen()%3))
                {
                    // backtracking to the previous solution.
                    memcpy(currentPermutation[0], CPopt[0], 4 * left_num);
                    memcpy(currentPermutation[1], CPopt[1], 4 * right_num);
                    maintainpositions();
                    currentsol=getcurrentsolution();
                    computeM(0,left_num);
                    computeM(1,right_num);
                    cnt+=15;    // this resets the tabu table.
                }
                if((USEDP_&&!dp3())||(!USEDP_))
                {
                    opt=min(opt,currentsol);
                    random_perturbation(0,0.1+(gen()%100)/500.0);
                    random_perturbation(1,0.1+(gen()%100)/500.0);
                    // random_perturbation(lr);
                    maintainpositions();
                    currentsol=getcurrentsolution();
                    computeM(0,left_num);
                    computeM(1,right_num);
                    // computeM(lr^1,(lr?left_num:right_num));
                    runtimes++;
                    cnt+=15;    // this resets the tabu table.
                }
            }
        }
        presol = currentsol;
        lr ^= 1;
        ed = system_clock::now();
        if (counttime(st, ed) > timelimit)
            break;
    }
}
void prob::checkcp()
{
    cout<<"left--right"<<endl;
    for(int i=0;i<max(left_num,right_num);i++)
    {
        string a=(i<left_num?to_string(currentPermutation[0][i]):"");
        string b=(i<right_num?to_string(currentPermutation[1][i]):"");
        cout<<setw(5)<<a<<' '<<b<<endl;
    }
    vector<int> checkl(left_num),checkr(right_num);
    for(int i=0;i<left_num;i++)
    {
        checkl[currentPermutation[0][i]]++;
    }
    for(int i=0;i<right_num;i++)
    {
        checkr[currentPermutation[1][i]]++;
    }
    cout<<"l ";
    for(auto e:checkl)
    {
        if(e!=1)    cout<<"is not permutation\n";
    }
    cout<<"r ";
    for(auto e:checkr)
    {
        if(e!=1)    cout<<"is not permutation\n";
    }

    cout<<"input to continue\n";
    char aaa;
    cin>>aaa;
}
int main(int argc,char **argv)
{
    string filepath(argv[1]);
    string respath=string(argv[2]);
    ofstream testresult(respath,ios::app);
    ios::sync_with_stdio(false);
    static prob a;
    a.readgraph(filepath);
    a.random_perturbation(0);
    a.random_perturbation(1);
    a.maintainpositions();
    a.currentsol = a.getcurrentsolution();
    // cout << "init number of crossings: " << currentsol << endl;
    // testresult<<filepath<<endl;
    // testresult << "init number of crossings: " << a.currentsol << endl;
    // st=system_clock::now();
    // getcurrentsolution();
    // ed=system_clock::now();
    // cout<<"getsolution time: "<<counttime()<<endl;

    st = system_clock::now();
    a.localsearch(false,true);
    a.currentsol = a.getcurrentsolution();
    a.opt=min(a.opt,a.currentsol);
    ed = system_clock::now();
    // testresult<<"main.cpp opt: "<<a.opt<<endl;
    testresult<<a.opt<<endl;
    //////////////////////////////////////////////////////////
    // dpls result                                          //
    //////////////////////////////////////////////////////////
    // ifstream inputfile("./result2.out");
    // for (int i = 0; i < a.left_num; i++)
    //     inputfile >> a.currentPermutation[0][i];
    // for (int i = 0; i < a.right_num; i++)
    //     inputfile >> a.currentPermutation[1][i];
    // a.currentsol = a.opt = 0x7fffffff;
    // a.currentsol=a.getcurrentsolution();
    // testresult<<"NewGraph.cpp opt: "<<a.currentsol<<"\n\n";
}