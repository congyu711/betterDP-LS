#include <chrono>
#include <iostream>
using namespace std;
using namespace chrono;
auto st=system_clock::now();
auto ed=system_clock::now();
double counttime(system_clock::time_point start=st,system_clock::time_point end=ed)
{
    auto duration = duration_cast<microseconds>(end - start);
    return double(duration.count()) * microseconds::period::num / microseconds::period::den;
}
// int main()
// {
//     auto start=system_clock::now();
//     int a=5;
//     for(int i=1;i<=100010;i++)  cout<<"sdf"<<' ';
//     auto end=system_clock::now();
//     cout<<counttime(start,end)<<endl;
// }