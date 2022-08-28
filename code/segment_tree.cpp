#include <bits/stdc++.h>
using namespace std;
// int mm, nn;
#define lson l, m, rt << 1
#define rson m + 1, r, rt << 1 | 1
// int a[100010];
class node
{
public:
    int l, r;
    int sum;
    int col;

    void init(int p,vector<int> *a)
    {
        col = 0;
        l = r = p;
        // scanf("%d", &sum);
        sum = a->at(p);
    }
};

node operator+(const node &l, const node &r)
{
    node ans;
    ans.l = l.l;
    ans.r = r.r;
    ans.sum = l.sum + r.sum;
    ans.col = 0;
    return ans;
}

class segtree
{
public:
    vector<int> *arr;
    int l, r;
    vector<node> z;
    segtree(int _l, int _r,vector<int> *_arr) : l(_l), r(_r), arr(_arr)
    {
        z.resize(4 * (_r - _l + 2));
        // segment tree of values
        build(_l,_r,1);
    }

    void update(int rt)
    {
        z[rt] = z[rt << 1] + z[rt << 1 | 1];
    }

    void color(int rt, int v)
    {
        z[rt].col += v;
        z[rt].sum += (z[rt].r - z[rt].l + 1) * v;
    }

    void push_col(int l, int r, int rt)
    {
        if (z[rt].col != 0)
        {
            color(rt << 1, z[rt].col);
            color(rt << 1 | 1, z[rt].col);
            z[rt].col = 0;
        }
    }

    void build(int l, int r, int rt)
    {
        if (l == r)
        {
            z[rt].init(l,arr);
            return;
        }
        int m = (l + r) >> 1;
        build(lson);
        build(rson);
        update(rt);
    }

    node query(int l, int r, int rt, int ql, int qr)
    {
        if (ql <= l && r <= qr)
            return z[rt];
        push_col(l, r, rt);
        int m = (l + r) >> 1;
        if (ql <= m)
        {
            if (m < qr)
                return query(lson, ql, qr) + query(rson, ql, qr);
            else
                return query(lson, ql, qr);
        }
        else
            return query(rson, ql, qr);
    }

    void modify(int l, int r, int rt, int ml, int mr, int v)
    {
        if (ml <= l && r <= mr)
        {
            color(rt, v);
            return;
        }
        push_col(l, r, rt);
        int m = (l + r) >> 1;
        if (ml <= m)
            modify(lson, ml, mr, v);
        if (m < mr)
            modify(rson, ml, mr, v);
        update(rt);
    }
};

// node z[100000 << 2];

// int x, y, k;
// int main()
// {
//     cin >> nn >> mm;
//     for (int i = 1; i <= nn; i++)
//         cin >> a[i];
//     build(1, nn, 1);
//     int opt;
//     for (int i = 1; i <= mm; i++)
//     {
//         cin >> opt;
//         if (opt == 1)
//         {
//             cin >> x >> y >> k;
//             modify(1, nn, 1, x, y, k);
//         }
//         if (opt == 2)
//         {
//             cin >> x >> y;
//             cout << query(1, nn, 1, x, y).sum << endl;
//         }
//     }
// }