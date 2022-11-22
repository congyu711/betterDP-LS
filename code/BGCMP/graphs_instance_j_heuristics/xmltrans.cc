#include<bits/stdc++.h>
using namespace std;

int main(int argc, char **argv)
{
    string str(argv[1]),filename;
    filename=str.substr(0,str.length()-4);
    cout<<filename<<endl;
    ifstream fin(str+".xml");
    ofstream fout(str+".txt");
    fin.ignore(numeric_limits<streamsize>::max(),'\n');
    getline(fin,str);

    int ecnt=0;
    for(int i=0;i<str.length();i++)
    {
        if(str[i]=='=')
        {
            ecnt++;
            if(ecnt==2)
            {
                string tmp; i++;
                while (str[i+1]!='"')
                {
                    i++;
                    tmp+=str[i];
                }
                
            }
            else if(ecnt==3)
            {
                string tmp; i++;
                while (str[i+1]!='"')
                {
                    i++;
                    tmp+=str[i];
                }
            }
            else if(ecnt==4)
            {

            }
        }
    }
}