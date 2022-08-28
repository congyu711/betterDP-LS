import os
import subprocess

workingdictionary='/home/congyu/dpls/code/'
path="/home/congyu/dpls/instances/G_21_06/"
# path="/home/congyu/dpls/code/myinstances/"
# path="/home/congyu/dpls/instances_dblp/stallman_reduced/"
output = subprocess.check_output(["ls",path])
files = output.decode("utf-8")
files = files.split('\n')
for f in files:
    print(f)
    if(len(f)==0 or f=='gen' or f=='gen.cpp' or f[0]=='.'):   continue
    # os.system("./NewGraph "+path+f)
    # os.system("./main "+path+f+' '+workingdictionary+'testresult.out')
    subprocess.Popen(workingdictionary+"main "+path+f+' '+workingdictionary+"testresult1.out",shell=True)
    subprocess.Popen(workingdictionary+"main "+path+f+' '+workingdictionary+"testresult2.out",shell=True)
    subprocess.Popen(workingdictionary+"main "+path+f+' '+workingdictionary+"testresult3.out",shell=True)
    subprocess.Popen(workingdictionary+"main "+path+f+' '+workingdictionary+"testresult4.out",shell=True)
    subprocess.Popen(workingdictionary+"main "+path+f+' '+workingdictionary+"testresult5.out",shell=True)
    subprocess.Popen(workingdictionary+"main "+path+f+' '+workingdictionary+"testresult6.out",shell=True)
    subprocess.Popen(workingdictionary+"main "+path+f+' '+workingdictionary+"testresult7.out",shell=True)
    subprocess.Popen(workingdictionary+"main "+path+f+' '+workingdictionary+"testresult8.out",shell=True)
    subprocess.Popen(workingdictionary+"main "+path+f+' '+workingdictionary+"testresult9.out",shell=True)
    subprocess.Popen(workingdictionary+"main "+path+f+' '+workingdictionary+"testresult10.out",shell=True).wait()
tmp=[]
for i in range(1,11):
    with open("testresult"+str(i)+".out") as f:
        idx=0
        for line in f.readlines():
            if(i==1):tmp.append(line.strip('\n'))
            else:
                tmp[idx]=tmp[idx]+", "+line.rstrip('\n')
                if(i==10): tmp[idx]=tmp[idx]+'\n'
                idx=idx+1

with open("testresultALL.csv",'a') as f:
    f.writelines(tmp)