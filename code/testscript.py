import os
import subprocess

pwd=os.popen("pwd").read().strip('\n')+'/'
print(pwd)
datapath=pwd+"/BGCMP/instances/"
# datapath="/home/congyu/dpls/code/myinstances/"
# datapath="/home/congyu/dpls/instances_dblp/stallman_reduced/"
output = subprocess.check_output(["ls",datapath])
files = output.decode("utf-8")
files = files.split('\n')
for f in files:
    print(f)
    if(len(f)==0 or f=='gen' or f=='gen.cpp' or f[0]=='.'):   continue
    # os.system("./NewGraph "+datapath+f)
    # os.system("./main "+datapath+f+' '+pwd+'testresult.out')
    subprocess.Popen(pwd+"main "+datapath+f+' '+pwd+"testresult1.out",shell=True)
    subprocess.Popen(pwd+"main "+datapath+f+' '+pwd+"testresult2.out",shell=True)
    subprocess.Popen(pwd+"main "+datapath+f+' '+pwd+"testresult3.out",shell=True)
    subprocess.Popen(pwd+"main "+datapath+f+' '+pwd+"testresult4.out",shell=True).wait()
    # subprocess.Popen(pwd+"main "+datapath+f+' '+pwd+"testresult5.out",shell=True)
    # subprocess.Popen(pwd+"main "+datapath+f+' '+pwd+"testresult6.out",shell=True)
    # subprocess.Popen(pwd+"main "+datapath+f+' '+pwd+"testresult7.out",shell=True)
    # subprocess.Popen(pwd+"main "+datapath+f+' '+pwd+"testresult8.out",shell=True)
    # subprocess.Popen(pwd+"main "+datapath+f+' '+pwd+"testresult9.out",shell=True)
    # subprocess.Popen(pwd+"main "+datapath+f+' '+pwd+"testresult10.out",shell=True).wait()
tmp=[]
for i in range(1,5):
    with open("testresult"+str(i)+".out") as f:
        idx=0
        for line in f.readlines():
            if(i==1):tmp.append(line.strip('\n'))
            else:
                tmp[idx]=tmp[idx]+", "+line.rstrip('\n')
                if(i==4): tmp[idx]=tmp[idx]+'\n'
                idx=idx+1

with open("testres.csv",'a') as f:
    f.writelines(tmp)