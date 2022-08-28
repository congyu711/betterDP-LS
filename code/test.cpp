#include <bits/stdc++.h>
#include "segment_tree.cpp"
#include "counttime.cpp"
using namespace std;
int leftnum,rightnum;
vector<bool> movableL(1000),movableR(1000);
vector<vector<int>> lr[500],rl[500];
int main()
{
    auto st=clock();
    auto stt=system_clock::now();
    vector<int> arr={0,1,2,5,6,4,3,8,7,9};
    segtree seg(1,8,&arr);
    cout<<seg.query(1,8,1,1,8).sum<<endl;
    for(int i=1;i<=100000000;i++)
        seg.modify(1,8,1,1,2,100);
    cout<<seg.query(1,8,1,1,3).sum<<endl;
    auto ed=clock();
    auto edd=system_clock::now();
    cout<<counttime(stt,edd)<<"------"<<0.001*difftime(ed,st)<<endl;

}