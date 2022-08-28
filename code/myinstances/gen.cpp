#include<bits/stdc++.h>
using namespace std;
int num,n,m;
int main()
{
    cout<<"number of instances: ";
    cin>>num;
    cout<<"number of vertices: ";
    cin>>n;
    cout<<"number of edges(approximate): ";
    cin>>m;
    mt19937 gen(time(nullptr));
    for(int i=1;i<=num;i++)
    {
        vector<set<int>> adj(n);
        for(int i=0;i<m;i++)
        {
            adj[gen()%n].insert(gen()%n);
        }
        ofstream fout(to_string(i)+".txt");
        fout<<2<<endl;
        fout<<n<<' '<<n<<endl;
        for(int i=0;i<n;i++)
        {
            fout<<(i<5?1:0)<<' '<<i<<' ';
            int aux=0;
            for(auto e:adj[i])
            {
                fout<<e+n;
                if((++aux)!=adj[i].size())  fout<<' ';
            }
            if(adj[i].size()==0)    fout<<n+gen()%n;
            fout<<endl;
        }
        for(int i=0;i<n;i++)
        {
            fout<<(i<5?1:0)<<' '<<(i<5?i:i+n)<<endl;
        }
        fout.close();
    }
    
}