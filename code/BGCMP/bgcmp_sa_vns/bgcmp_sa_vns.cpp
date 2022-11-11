/* Program: An implementation of the iterated simulated annealing (SA) and 
         variable neighborhood search (VNS) hybrid (SA-VNS) for the 
         bipartite graph crossing minimization problem (BGCMP): 
         Find an embedding the parts of a bipartite graph along two 
         parallel lines so that the number of edge crossings is minimized. 
         Each edge is drawn as the straight line connecting its endvertices. 
   Author: Gintaras Palubeckis
   Date: 2016-12-11
   Language: C++
   Some of the input data are supplied through the parameters and the rest 
     through the (input) file. An instance of the BGCMP in this file is 
     represented by the list of edges of the graph. The program terminates
     when a specified time limit is reached.
   Parameters:
     - input file name;
     - output file name;
     - seed for random number generator;
     - time limit (in seconds);
     - a pointer to structure 'Results' for writing down the solution found 
       and the performance characteristics. The structure is defined in 
       file 'bgcmp_sa_vns.h'. It is possible to have the last parameter null: 
       in this case, information is directed only to the output file.
   Examples of invocation:
       either
     Results *pres;
     char in_file_name[80],out_file_name[80];
     double seed=5000;
     pres=(Results *)calloc(1,sizeof(Results));
     strcpy(in_file_name,"D:\\data\\Mesh30x15.txt");
     strcpy(out_file_name,"D:\\temp\\Mesh30x15.res");
     bgcmp_sa_vns(in_file_name,out_file_name,seed,300,pres);
       or just
     char in_file_name[80],out_file_name[80];
     double seed=5000;
     strcpy(in_file_name,"D:\\data\\Mesh30x15.txt");
     strcpy(out_file_name,"D:\\temp\\Mesh30x15.res");
     bgcmp_sa_vns(in_file_name,out_file_name,seed,300,NULL);
   Input file contains:
     - the number of vertices in the first part;
     - the number of vertices in the second part;
     - the number of edges;
     - the edges represented by their endvertices. 
   Example of the input file:
 225  225 855
   9  288 
   9  247 
  58  288 
  58  327 
  58  324 
 164  327 
 164  298 
 164  449 
  49  298 
  49  339 
  49  399 
  62  339 
  62  287 
  62  441 
  37  287 
  37  355 
  37  359 
  48  355 
  48  354 
  48  282 
  68  354 
  68  395 
  77  247 
  77  324 
  77  288 
  77  313 
 129  324 
 129  449 
 129  327 
 129  414 
  71  449 
  71  399 
...
...
...
*/



#include <alloc.h>
#include <process.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "bgcmp_sa_vns.h"



 double random(double *seed,double coef)
 {double rd,rf;
  rd=16807*(*seed);rf=floor(rd/coef);
  *seed=rd-rf*coef;
  return(*seed/(coef+1));
 }


 double take_time(int *time_values,clock_t start)
 {int i;
  int hours,mins;
  long longsecs;
  double elapsed_sec;
  clock_t end;
  end=clock();
  elapsed_sec=(end-start)/CLK_TCK;
  longsecs=elapsed_sec;
  for (i=1;i<=4;i++) time_values[i]=0;
  hours=(int)(longsecs/3600);
  if (hours>0)     /* more than an hour  */
     {time_values[1]=hours;longsecs-=hours*3600;}
  mins=(int)(longsecs/60);
  if (mins>0)     /* more than a minute */
     {time_values[2]=mins;longsecs-=mins*60;}
  time_values[3]=(int)longsecs;
  time_values[4]=elapsed_sec*1000-(long)elapsed_sec*1000;
  return elapsed_sec;
 }


 long get_best_value(int n1,int n2,Instance *pinst,Solution *psol)
 {int i,j,m,q;
  int u,v,w;
  int pos;
  long sol_val=0;
// Calculation of the value of the best solution
  for (m=1;m<=n2;m++) (psol+(psol+m)->gbp2)->pi2=m;
  for (m=1;m<n1;m++)
     {u=(psol+m)->gbp1;
      for (q=m+1;q<=n1;q++)
         {v=(psol+q)->gbp1;
          for (i=1;i<=(psol+u)->aln1;i++)
             {w=a1(u,i);pos=(psol+w)->pi2;
              for (j=1;j<=(psol+v)->aln1;j++) {if (pos>(psol+a1(v,j))->pi2) sol_val++;}
             }
         }
     }
  return sol_val;
 }


 void random_start(int n1,int n2,double coef,double *seed,Solution *psol)
 {int i,k;
  int r;
// Randomly generating first permutation
  for (i=1;i<=n1;i++) (psol+i)->a=i;
  for (i=1;i<=n1;i++)
     {k=random(seed,coef)*(n1-i+1);k+=i;
      r=(psol+k)->a;(psol+i)->p1=r;(psol+r)->pi1=i;
      (psol+k)->a=(psol+i)->a;
     }
// Randomly generating second permutation
  for (i=1;i<=n2;i++) (psol+i)->a=i;
  for (i=1;i<=n2;i++)
     {k=random(seed,coef)*(n2-i+1);k+=i;
      r=(psol+k)->a;(psol+i)->p2=r;(psol+r)->pi2=i;
      (psol+k)->a=(psol+i)->a;
     }
 }


 long objective_value(int n1,int n2,Instance *pinst,Solution *psol)
 {int i,l,m;
  int u,v;
  int add;
  long sol_val=0;

  for (u=1;u<=n1;u++) (psol+u)->zln=0;
  for (m=1;m<=n2;m++)
     {(psol+m)->a=0;
      u=(psol+m)->p2;
      for (i=1;i<=(psol+u)->aln2;i++)
         {v=a2(u,i);
          ((psol+v)->zln)++;z(v,(psol+v)->zln)=m;
         }
     }
//
  for (m=1;m<=n1;m++)
     {u=(psol+m)->p1;
      if ((psol+u)->zln==0) continue;
      for (i=1;i<=(psol+u)->zln;i++) sol_val+=((psol+z(u,i))->a);
      i=z(u,(psol+u)->zln)-1;
      l=(psol+u)->zln-1;
      add=1;
      while (i>0)
         {(psol+i)->a+=add;
          i--;
          if (l>0 && z(u,l)>i) {add++;l--;}
         }
     }
  return sol_val;
 }


/* long objective_value(int n1,int n2,Instance *pinst,Solution *psol)
 {int i,k,l,m,q;
  int u,v;
  int cont;
  long sol_val=0;

  for (u=1;u<=n1;u++) (psol+u)->zln=0;
  for (m=1;m<=n2;m++)
     {u=(psol+m)->p2;
      for (i=1;i<=(psol+u)->aln2;i++)
         {v=a2(u,i);
          ((psol+v)->zln)++;z(v,(psol+v)->zln)=m;
         }
     }
//
  for (m=1;m<n1;m++)
     {u=(psol+m)->p1;
      for (q=m+1;q<=n1;q++)
         {v=(psol+q)->p1;
          if ((psol+u)->zln==0 || (psol+v)->zln==0) continue;
          cont=1;k=l=0;
          while (cont>0)
             {if (z(u,k+1)<z(v,l+1))
                 {sol_val+=l;
                  if (k+1==(psol+u)->zln) cont=0; else k++;
                 }
               else if (z(u,k+1)>z(v,l+1))
                 {if (l+1==(psol+v)->zln) {sol_val+=((l+1)*((psol+u)->zln-k));cont=0;}
                   else l++;
                 }
                else
                 {if (k+1==(psol+u)->zln) {sol_val+=l;cont=0;}
                   else if (l+1==(psol+v)->zln) {sol_val+=((l+1)*((psol+u)->zln-k)-1);cont=0;}
                    else {sol_val+=l;k++;l++;}
                 }
             } // end while
         }
     }
  return sol_val;
 }*/


 void sort_integers(int *weaux,int *we,long s_ind,long e_ind)
 {int i,j,k,l;
  int m_ind;
  int e;
  int gap;

// This code implements the well-known merge algorithm for the sorting problem.
// The array submitted for sorting is we(.) and is sorted in nondecreasing order
  gap=e_ind-s_ind;
  if (gap<7)
     {for (i=s_ind;i<e_ind;i++)
          for (j=i;j>s_ind && *(we+j-1)>*(we+j);j--)
             {e=*(we+j-1);*(we+j-1)=*(we+j);*(we+j)=e;}
      return;
     }
  m_ind=(s_ind+e_ind)/2;
  sort_integers(we,weaux,s_ind,m_ind);
  sort_integers(we,weaux,m_ind,e_ind);
  if (*(weaux+m_ind-1)<=*(weaux+m_ind))
     {for (i=s_ind;i<e_ind;i++) *(we+i)=*(weaux+i);
      return;
     }
  for (i=s_ind,k=s_ind,l=m_ind;i<e_ind;i++)
     {if (l>=e_ind || k<m_ind && *(weaux+k)<=*(weaux+l))
         {*(we+i)=*(weaux+k);k++;}
       else {*(we+i)=*(weaux+l);l++;}
     }
 }


int find_position_1(int ind,int right,Instance *pinst)
 {int left=1,middle;
  if (ind<h1(0,1)) return 0;
   else if (ind==h1(0,1)) {h1(0,0)=1;return 0;}
  if (ind>h1(0,right)) return right;
   else if (ind==h1(0,right)) {h1(0,0)=1;return right-1;}
  while (left<=right)
     {if ((right-left)<=1) return left;
      middle=(left+right)/2;
      if (ind==h1(0,middle)) {h1(0,0)=1;return middle-1;}
      if (ind<h1(0,middle)) right=middle; else left=middle;
     }
  return -1;
 }


int find_position_2(int ind,int right,Instance *pinst)
 {int left=1,middle;
  if (ind<h2(0,1)) return 0;
   else if (ind==h2(0,1)) {h2(0,0)=1;return 0;}
  if (ind>h2(0,right)) return right;
   else if (ind==h2(0,right)) {h2(0,0)=1;return right-1;}
  while (left<=right)
     {if ((right-left)<=1) return left;
      middle=(left+right)/2;
      if (ind==h2(0,middle)) {h2(0,0)=1;return middle-1;}
      if (ind<h2(0,middle)) right=middle; else left=middle;
     }
  return -1;
 }


 int get_gain_1(int r,int s,int k,int l,Instance *pinst,Solution *psol)
 {int i,q;
  int v,u;
  int dif=0,deg_sum=0,ad_sum1=0,ad_sum2=0;
  int count=0;
  int pos1,pos2;

// Permutation interval between swapped vertices (inclusively the right end) is processed. 
// All the neighbors of vertices in this interval are identified
  for (q=k+1;q<=l;q++)
     {v=(psol+q)->p1;deg_sum+=((psol+v)->aln1);
      for (i=1;i<=(psol+v)->aln1;i++)
         {u=a1(v,i);
          if ((psol+u)->ed<=0)
             {(psol+u)->ed=1;
              count++;
              (psol+count)->ad=u;
             }
           else ((psol+u)->ed)++;
         }
     }
// Positions of the neighbors of the vertex r are sorted from left to right
  if ((psol+r)->aln1>0)
     {for (i=1;i<=(psol+r)->aln1;i++)
         {u=a1(r,i);
          z(0,i)=(psol+u)->pi2;h1(0,i)=z(0,i);
         }
     }
  if ((psol+r)->aln1>1) sort_integers(pinst->z,pinst->h1,1,(psol+r)->aln1+1);
// Positions of the neighbors of the vertex s are sorted from left to right
  if ((psol+s)->aln1>0)
     {for (i=1;i<=(psol+s)->aln1;i++)
         {u=a1(s,i);
          z(0,i)=(psol+u)->pi2;h2(0,i)=z(0,i);
         }
     }
  if ((psol+s)->aln1>1) sort_integers(pinst->z,pinst->h2,1,(psol+s)->aln1+1);
// Running through the neighbors of the vertices in the interval
  for (i=1;i<=count;i++)
     {u=(psol+i)->ad;
// Operations for the left vertex in the pair
      if ((psol+r)->aln1>0)
         {pos1=find_position_1((psol+u)->pi2,(psol+r)->aln1,pinst);
          if (h1(0,0)==1) {h1(0,0)=0;ad_sum1+=((psol+u)->ed);}
         }
       else pos1=0;
// Operations for the right vertex in the pair
      if ((psol+s)->aln1>0)
         {pos2=find_position_2((psol+u)->pi2,(psol+s)->aln1,pinst);
          if (h2(0,0)==1) {h2(0,0)=0;ad_sum2+=(2*pos2-(psol+u)->ed);}
         }
       else pos2=0;
      dif+=((pos1-pos2)*(psol+u)->ed);
      (psol+u)->ed=0;
     }
  dif*=2;
  dif+=(((psol+s)->aln1-(psol+r)->aln1)*deg_sum+ad_sum1+ad_sum2-(psol+s)->s1);
  return dif;
 }


 int get_gain_2(int r,int s,int k,int l,Instance *pinst,Solution *psol)
 {int i,q;
  int v,u;
  int dif=0,deg_sum=0,ad_sum1=0,ad_sum2=0;
  int count=0;
  int pos1,pos2;

// This procedure is similar to get_gain_1(...). 
// It performs the same operations, but with respect to the second part of the bipartite graph
  for (q=k+1;q<=l;q++)
     {v=(psol+q)->p2;deg_sum+=((psol+v)->aln2);
      for (i=1;i<=(psol+v)->aln2;i++)
         {u=a2(v,i);
          if ((psol+u)->ed<=0)
             {(psol+u)->ed=1;
              count++;
              (psol+count)->ad=u;
             }
           else ((psol+u)->ed)++;
         }
     }
  if ((psol+r)->aln2>0)
     {for (i=1;i<=(psol+r)->aln2;i++)
         {u=a2(r,i);
          z(0,i)=(psol+u)->pi1;h1(0,i)=z(0,i);
         }
     }
  if ((psol+r)->aln2>1) sort_integers(pinst->z,pinst->h1,1,(psol+r)->aln2+1);
  if ((psol+s)->aln2>0)
     {for (i=1;i<=(psol+s)->aln2;i++)
         {u=a2(s,i);
          z(0,i)=(psol+u)->pi1;h2(0,i)=z(0,i);
         }
     }
  if ((psol+s)->aln2>1) sort_integers(pinst->z,pinst->h2,1,(psol+s)->aln2+1);
  for (i=1;i<=count;i++)
     {u=(psol+i)->ad;
// Operations for the left vertex in the pair
      if ((psol+r)->aln2>0)
         {pos1=find_position_1((psol+u)->pi1,(psol+r)->aln2,pinst);
          if (h1(0,0)==1) {h1(0,0)=0;ad_sum1+=((psol+u)->ed);}
         }
       else pos1=0;
// Operations for the right vertex in the pair
      if ((psol+s)->aln2>0)
         {pos2=find_position_2((psol+u)->pi1,(psol+s)->aln2,pinst);
          if (h2(0,0)==1) {h2(0,0)=0;ad_sum2+=(2*pos2-(psol+u)->ed);}
         }
       else pos2=0;
      dif+=((pos1-pos2)*(psol+u)->ed);
      (psol+u)->ed=0;
     }
  dif*=2;
  dif+=(((psol+s)->aln2-(psol+r)->aln2)*deg_sum+ad_sum1+ad_sum2-(psol+s)->s2);
  return dif;
 }


 void perform_interchange(int r,int s,int k,int l,int part,Solution *psol)
 {
// Exchanging positions of the selected vertices
  if (part==1)
     {(psol+k)->p1=s;(psol+l)->p1=r;
      (psol+r)->pi1=l;(psol+s)->pi1=k;
     }
   else
     {(psol+k)->p2=s;(psol+l)->p2=r;
      (psol+r)->pi2=l;(psol+s)->pi2=k;
     }
 }


 int init_temperature(int n1,int n2,int sam_size,double coef,
      double *seed,Instance *pinst,Solution *psol)
 {int i,j,k,l;
  int r,s;
  int t_max=-1,dif;
  int moves1,moves2;

  moves1=sam_size/2;moves2=sam_size-moves1;
// Pairwise interchanges of vertices of the first part
  for (j=1;j<=moves1;j++)
     {r=random(seed,coef)*n1;if (r<n1) r++;k=(psol+r)->pi1;
      s=r;
      while (s==r) {s=random(seed,coef)*n1;if (s<n1) s++;}
      l=(psol+s)->pi1;
      if (l<k) {i=k;k=l;l=i;i=r;r=s;s=i;}
      dif=get_gain_1(r,s,k,l,pinst,psol);
      if (dif<0) dif=-dif;
      if (dif>t_max) t_max=dif;
     }
// Pairwise interchanges of vertices of the second part
  for (j=1;j<=moves2;j++)
     {r=random(seed,coef)*n2;if (r<n2) r++;k=(psol+r)->pi2;
      s=r;
      while (s==r) {s=random(seed,coef)*n2;if (s<n2) s++;}
      l=(psol+s)->pi2;
      if (l<k) {i=k;k=l;l=i;i=r;r=s;s=i;}
      dif=get_gain_2(r,s,k,l,pinst,psol);
      if (dif<0) dif=-dif;
      if (dif>t_max) t_max=dif;
     }
  return t_max;
 }


 void store_best_sol_sa(int n1,int n2,Solution *psol)
 {int k;
// Best solution (found by SA) in an iteration of SA-VNS is updated
  for (k=1;k<=n1;k++) (psol+k)->bp1=(psol+k)->p1;
  for (k=1;k<=n2;k++) (psol+k)->bp2=(psol+k)->p2;
 }


 void store_best_sol_vns(int n1,int n2,Solution *psol)
 {int k;
// Best solution (found by VNS) in an iteration of SA-VNS is updated
  for (k=1;k<=n1;k++) (psol+k)->bp1=(psol+k)->p1;
  for (k=1;k<=n2;k++) (psol+k)->bp2=(psol+k)->p2;
  ((psol+8)->perf)++;
 }


 void store_best_sol_glob(int n1,int n2,int st,int *time_values_opt,
      double *time_to_opt,clock_t start_time,Solution *psol)
 {int k;
// Globally best solution (over all iterations of SA-VNS) is updated
  for (k=1;k<=n1;k++) (psol+k)->gbp1=(psol+k)->bp1;
  for (k=1;k<=n2;k++) (psol+k)->gbp2=(psol+k)->bp2;
  ((psol+1)->perf)++;(psol+2)->perf=st;
  *time_to_opt=take_time(time_values_opt,start_time);
 }


 long sa(int n1,int n2,int tred,int tlength,long best_value,long st,
      long time_limit,double time_limit_sa,double t_max,double part_prob,
      double coef,int *stop_cond,long *gl_best_value,double *seed1,
      double *seed2,int *time_values_opt,double *time_to_opt,
      clock_t start_time,Instance *pinst,Solution *psol)
 {int i,k,l,r,s;
  int tr,tl;
  int accept,succ;
  int dif;
  int part;
  int at;
  long sol_value;
  double elapsed_time;
  double t;
  double c_factor;
  double db,md;
  clock_t end_time;

  if (st==1) sol_value=best_value;
   else
     {random_start(n1,n2,coef,seed1,psol);
      sol_value=best_value=objective_value(n1,n2,pinst,psol);
      store_best_sol_sa(n1,n2,psol);
     }
  c_factor=COOLING_FACTOR;
  t=t_max;
// Simulated annealing ...
  for (tr=1;tr<=tred;tr++)
     {for (tl=1;tl<=tlength;tl++)
// Selecting a pair of vertices at random 
         {db=random(seed2,coef);
          if (db<=part_prob) part=1; else part=2;
          at=0;succ=0;
          while (at<TRIAL_LIMIT)
             {at++;
              if (part==1)
                 {r=random(seed2,coef)*n1;if (r<n1) r++;
                  s=r;
                  while (s==r) {s=random(seed2,coef)*n1;if (s<n1) s++;}
                  if ((psol+r)->aln1==0 || (psol+s)->aln1==0) continue;
                  k=(psol+r)->pi1;l=(psol+s)->pi1;
                  if (k>l) {i=r;r=s;s=i;i=k;k=l;l=i;}
                  succ=1;break;
                 }
              else
                 {r=random(seed2,coef)*n2;if (r<n2) r++;
                  s=r;
                  while (s==r) {s=random(seed2,coef)*n2;if (s<n2) s++;}
                  if ((psol+r)->aln2==0 || (psol+s)->aln2==0) continue;
                  k=(psol+r)->pi2;l=(psol+s)->pi2;
                  if (k>l) {i=r;r=s;s=i;i=k;k=l;l=i;}
                  succ=1;break;
                 }
             } // end while
          if (succ==0) {printf("The graph has too many isolated vertices! \n");exit(1);}
// Getting the gain value 
          if (part==1) dif=get_gain_1(r,s,k,l,pinst,psol);
           else dif=get_gain_2(r,s,k,l,pinst,psol);
// Deciding on acceptance or rejection of the move
          if (dif<=0) accept=2;
           else
             {md=random(seed2,coef);
              db=-dif/t;db=exp(db);
              if (md<=db) accept=1; else accept=0;
             }
          if (accept==0) continue;
// Performing interchange move
          sol_value+=dif;
          perform_interchange(r,s,k,l,part,psol);
          if (accept<2) continue;
          if (sol_value<best_value)
             {best_value=sol_value;((psol+14)->perf)++;
              store_best_sol_sa(n1,n2,psol);
              if (best_value<*gl_best_value)
                 {store_best_sol_glob(n1,n2,st,time_values_opt,time_to_opt,start_time,psol);
                  *gl_best_value=best_value;
                 }
             }
         }
      end_time=clock();
      elapsed_time=(end_time-start_time)/CLK_TCK;
      if (st==1) {if (elapsed_time>=time_limit_sa) {*stop_cond=1; break;}}
       else if (elapsed_time>=time_limit) {*stop_cond=2; break;}
// Geometrically decreasing the temperature
      t=c_factor*t;
     }
  return best_value;
 }


 void prepare_data(int n1,int n2,Instance *pinst,Solution *psol)
 {int i,j;
  int u,v,w;
// Constructing "h1" matrix (for the first part of the bipartite graph) 
  for (u=1;u<n1;u++) for (v=u+1;v<=n1;v++) h1(u,v)=h1(v,u)=0;
  for (w=1;w<=n2;w++)
     {for (i=1;i<(psol+w)->aln2;i++)
         {u=a2(w,i);
          for (j=i+1;j<=(psol+w)->aln2;j++) {v=a2(w,j);(h1(u,v))++;(h1(v,u))++;}
         }
      (psol+w)->s2=(psol+w)->aln2*((psol+w)->aln2-1);
     }
  for (u=1;u<n1;u++) for (v=u+1;v<=n1;v++)
     {h1(u,v)=(psol+u)->aln1*(psol+v)->aln1-h1(u,v);
      h1(v,u)=h1(u,v);
     }
// Constructing "h2" matrix (for the second part of the bipartite graph) 
  for (u=1;u<n2;u++) for (v=u+1;v<=n2;v++) h2(u,v)=h2(v,u)=0;
  for (w=1;w<=n1;w++)
     {for (i=1;i<(psol+w)->aln1;i++)
         {u=a1(w,i);
          for (j=i+1;j<=(psol+w)->aln1;j++) {v=a1(w,j);(h2(u,v))++;(h2(v,u))++;}
         }
      (psol+w)->s1=(psol+w)->aln1*((psol+w)->aln1-1);
     }
  for (u=1;u<n2;u++) for (v=u+1;v<=n2;v++)
     {h2(u,v)=(psol+u)->aln2*(psol+v)->aln2-h2(u,v);
      h2(v,u)=h2(u,v);
     }
 }


 long init_insertions(int n1,int n2,Instance *pinst,Solution *psol)
 {int i,k,l,m,q;
  int u,v;
  int cross,cont;
  long sol_val=0;

// Initializing "b1" matrix (for the first part of the bipartite graph) 
  for (u=1;u<=n1;u++) (psol+u)->zln=0;
  for (m=1;m<=n2;m++)
     {u=(psol+m)->p2;
      for (i=1;i<=(psol+u)->aln2;i++)
         {v=a2(u,i);
          ((psol+v)->zln)++;z(v,(psol+v)->zln)=m;
         }
     }
  for (m=1;m<n1;m++)
     {u=(psol+m)->p1;
      for (q=m+1;q<=n1;q++)
         {v=(psol+q)->p1;
          if ((psol+u)->zln==0 || (psol+v)->zln==0) {b1(u,v)=b1(v,u)=0;continue;}
          cross=0;cont=1;k=l=0;
          while (cont>0)
             {if (z(u,k+1)<z(v,l+1))
                 {cross+=l;
                  if (k+1==(psol+u)->zln) cont=0; else k++;
                 }
               else if (z(u,k+1)>z(v,l+1))
                 {if (l+1==(psol+v)->zln) {cross+=((l+1)*((psol+u)->zln-k));cont=0;}
                   else l++;
                 }
                else
                 {if (k+1==(psol+u)->zln) {cross+=l;cont=0;}
                   else if (l+1==(psol+v)->zln) {cross+=((l+1)*((psol+u)->zln-k)-1);cont=0;}
                    else {cross+=l;k++;l++;}
                 }
             } // end while
          b1(u,v)=b1(v,u)=cross;sol_val+=cross;
         }
     }
// Initializing "b2" matrix (for the second part of the bipartite graph) 
  for (u=1;u<=n2;u++) (psol+u)->zln=0;
  for (m=1;m<=n1;m++)
     {u=(psol+m)->p1;
      for (i=1;i<=(psol+u)->aln1;i++)
         {v=a1(u,i);
          ((psol+v)->zln)++;z(v,(psol+v)->zln)=m;
         }
     }
  for (m=1;m<n2;m++)
     {u=(psol+m)->p2;
      for (q=m+1;q<=n2;q++)
         {v=(psol+q)->p2;
          if ((psol+u)->zln==0 || (psol+v)->zln==0) {b2(u,v)=b2(v,u)=0;continue;}
          cross=0;cont=1;k=l=0;
          while (cont>0)
             {if (z(u,k+1)<z(v,l+1))
                 {cross+=l;
                  if (k+1==(psol+u)->zln) cont=0; else k++;
                 }
               else if (z(u,k+1)>z(v,l+1))
                 {if (l+1==(psol+v)->zln) {cross+=((l+1)*((psol+u)->zln-k));cont=0;}
                   else l++;
                 }
                else
                 {if (k+1==(psol+u)->zln) {cross+=l;cont=0;}
                   else if (l+1==(psol+v)->zln) {cross+=((l+1)*((psol+u)->zln-k)-1);cont=0;}
                    else {cross+=l;k++;l++;}
                 }
             } // end while
          b2(u,v)=b2(v,u)=cross;
         }
     }
  return sol_val;
 }


 void perform_insertion(int n1,int n2,int r,int l,int side,
      Instance *pinst,Solution *psol)
 {int i,k,m;
  int u,v;
  int count;
  int pos;

// Performing insertion
// Updating matrices b1 and b2
// Considering the first part of the bipartite graph
  if (side==1)
     {k=(psol+r)->pi1;
      if (l<k)
// The case of inserting to the left
         {for (m=k;m>l;m--)
             {u=(psol+m-1)->p1;
              (psol+m)->p1=u;(psol+u)->pi1=m;
              b1(r,u)=h1(r,u)-b1(r,u);b1(u,r)=b1(r,u);
             }
         }
       else
// The case of inserting to the right
         {for (m=k;m<l;m++)
             {u=(psol+m+1)->p1;
              (psol+m)->p1=u;(psol+u)->pi1=m;
              b1(r,u)=h1(r,u)-b1(r,u);b1(u,r)=b1(r,u);
             }
         }
      (psol+l)->p1=r;(psol+r)->pi1=l;
      for (m=1;m<=n2;m++)
         {count=0;
          u=(psol+m)->p2;
          for (i=1;i<=(psol+u)->aln2;i++)
             {pos=(psol+a2(u,i))->pi1;
              if (l<k) {if (pos>l && pos<=k) count++;}
               else if (pos>=k && pos<l) count++;
             }
          for (i=1;i<=(psol+r)->aln1;i++)
             {v=a1(r,i);if (v==u) continue;
              if (l<k) {if ((psol+v)->pi2<(psol+u)->pi2) b2(u,v)-=count; else b2(u,v)+=count;}
               else if ((psol+v)->pi2<(psol+u)->pi2) b2(u,v)+=count; else b2(u,v)-=count;
              b2(v,u)=b2(u,v);
             }
         }
     }
// Considering the second part of the bipartite graph
   else
     {k=(psol+r)->pi2;
      if (l<k)
// The case of inserting to the left
         {for (m=k;m>l;m--)
             {u=(psol+m-1)->p2;
              (psol+m)->p2=u;(psol+u)->pi2=m;
              b2(r,u)=h2(r,u)-b2(r,u);b2(u,r)=b2(r,u);
             }
         }
       else
// The case of inserting to the right
         {for (m=k;m<l;m++)
             {u=(psol+m+1)->p2;
              (psol+m)->p2=u;(psol+u)->pi2=m;
              b2(r,u)=h2(r,u)-b2(r,u);b2(u,r)=b2(r,u);
             }
         }
      (psol+l)->p2=r;(psol+r)->pi2=l;
      for (m=1;m<=n1;m++)
         {count=0;
          u=(psol+m)->p1;
          for (i=1;i<=(psol+u)->aln1;i++)
             {pos=(psol+a1(u,i))->pi2;
              if (l<k) {if (pos>l && pos<=k) count++;}
               else if (pos>=k && pos<l) count++;
             }
          for (i=1;i<=(psol+r)->aln2;i++)
             {v=a2(r,i);if (v==u) continue;
              if (l<k) {if ((psol+v)->pi1<(psol+u)->pi1) b1(u,v)-=count; else b1(u,v)+=count;}
               else if ((psol+v)->pi1<(psol+u)->pi1) b1(u,v)+=count; else b1(u,v)-=count;
              b1(v,u)=b1(u,v);
             }
         }
     }
 }


 long local_search_insertions(int n1,int n2,Instance *pinst,Solution *psol)
 {int k,l;
  int r,s;
  int rbest,lbest,side;
  int best_change=-1;
  int delta;
  long value_change=0;

  while (best_change<0)
     {best_change=0;
// All pairs of type (vertex, its new position in the permutation) are searched
// Considering the first part of the bipartite graph
      for (r=1;r<=n1;r++) 
         {k=(psol+r)->pi1;
// Evaluating insertions to the left
          delta=0;
          for (l=k-1;l>=1;l--)
             {s=(psol+l)->p1;
              delta+=(h1(r,s)-2*b1(r,s));
              if (delta<best_change) {best_change=delta;rbest=r;lbest=l;side=1;}
             }
// Evaluating insertions to the right
          delta=0;
          for (l=k+1;l<=n1;l++)
             {s=(psol+l)->p1;
              delta+=(h1(r,s)-2*b1(r,s));
              if (delta<best_change) {best_change=delta;rbest=r;lbest=l;side=1;}
             }
         }
// Considering the second part of the bipartite graph
      for (r=1;r<=n2;r++) 
         {k=(psol+r)->pi2;
// Evaluating insertions to the left
          delta=0;
          for (l=k-1;l>=1;l--)
             {s=(psol+l)->p2;
              delta+=(h2(r,s)-2*b2(r,s));
              if (delta<best_change) {best_change=delta;rbest=r;lbest=l;side=2;}
             }
// Evaluating insertions to the right
          delta=0;
          for (l=k+1;l<=n2;l++)
             {s=(psol+l)->p2;
              delta+=(h2(r,s)-2*b2(r,s));
              if (delta<best_change) {best_change=delta;rbest=r;lbest=l;side=2;}
             }
         }
// If best_change<0, then an improving insertion is detected
      if (best_change<0)
// Performing insertion
         {perform_insertion(n1,n2,rbest,lbest,side,pinst,psol);
          value_change+=best_change;
         }
     }
  return value_change;
 }


 void shake(int n1,int n2,int interch_req,double share1,double share2,
      double coef,double *seed,Solution *psol)
 {int k,l,r,s;
  int req1,req2;
  int interch_count=0;
  int ind1,ind2,rem;

  req1=interch_req*share1;if (req1>n1/2) req1=n1/2;
  req2=interch_req*share2;if (req2>n2/2) req2=n2/2;
// Taking the best first permutation found thus far 
  for (k=1;k<=n1;k++) {(psol+k)->p1=(psol+k)->bp1;(psol+(psol+k)->p1)->pi1=k;(psol+k)->av=k;}
// Considering the first part of the bipartite graph
  rem=n1;
  while (interch_count<req1)
// Selecting a pair of vertices
     {ind1=random(seed,coef)*rem+1;
      ind2=random(seed,coef)*(rem-1)+1;if (ind2>=ind1) ind2++;
      r=(psol+ind1)->av;s=(psol+ind2)->av;
      interch_count++;
// Exchanging positions of the selected vertices
      k=(psol+r)->pi1;l=(psol+s)->pi1;
      (psol+k)->p1=s;(psol+l)->p1=r;(psol+r)->pi1=l;(psol+s)->pi1=k;
// Removing the selected pair from the set of available vertices
      if (ind2<rem)
         {(psol+ind1)->av=(psol+rem)->av;rem--;
          (psol+ind2)->av=(psol+rem)->av;rem--;
         }
       else
         {rem--;(psol+ind1)->av=(psol+rem)->av;rem--;}
     }
// Taking the best second permutation found thus far 
  interch_count=0;
  for (k=1;k<=n2;k++) {(psol+k)->p2=(psol+k)->bp2;(psol+(psol+k)->p2)->pi2=k;(psol+k)->av=k;}
// Considering the second part of the bipartite graph
  rem=n2;
  while (interch_count<req2)
// Selecting a pair of vertices
     {ind1=random(seed,coef)*rem+1;
      ind2=random(seed,coef)*(rem-1)+1;if (ind2>=ind1) ind2++;
      r=(psol+ind1)->av;s=(psol+ind2)->av;
      interch_count++;
// Exchanging positions of the selected vertices
      k=(psol+r)->pi2;l=(psol+s)->pi2;
      (psol+k)->p2=s;(psol+l)->p2=r;(psol+r)->pi2=l;(psol+s)->pi2=k;
// Removing the selected pair from the set of available vertices
      if (ind2<rem)
         {(psol+ind1)->av=(psol+rem)->av;rem--;
          (psol+ind2)->av=(psol+rem)->av;rem--;
         }
       else
         {rem--;(psol+ind1)->av=(psol+rem)->av;rem--;}
     }
 }


 long sa_vns(FILE *out,int n1,int n2,int dist_min,
      long time_limit,double time_limit_sa,int *time_values_opt,
      double *time_to_opt,double *seed1,clock_t start_time,
      Instance *pinst,Solution *psol)
 {int i;
  int ntot,nbig;
  int pred_st_count;
  int dist;
  int st=0;
  int outer_stop_cond=0,stop_cond=0;
  int tred;
  int tlength;
  int dist_max,dist_step,limit;
  long best_value,sol_value,gl_best_value;
  double seed2,seed3,seed4,coef;
  double t_max;
  double share1,share2;
  double fct,db;
  double time_limit_vns,time_limit_vns_restart;
  double elapsed_time,sa_elapsed_time,vns_elapsed_time;
  clock_t sa_start_time,vns_start_time,end;

  coef=2048;coef*=1024;coef*=1024;coef-=1;
  seed2=2*(*seed1);seed3=3*(*seed1);seed4=4*(*seed1);
  (psol+1)->perf=-1;(psol+10)->perf=0;
  (psol+8)->perf=0;(psol+14)->perf=0;(psol+3)->perf=0;

// Initialization
  prepare_data(n1,n2,pinst,psol);
  ntot=n1+n2;if (n1>=n2) nbig=n1; else nbig=n2;
  share1=((double)n1)/((double)ntot);share2=((double)n2)/((double)ntot);
  limit=n1/2+n2/2;
  time_limit_vns=time_limit-time_limit_sa;
  random_start(n1,n2,coef,seed1,psol);
  best_value=objective_value(n1,n2,pinst,psol);
  store_best_sol_sa(n1,n2,psol);
  store_best_sol_glob(n1,n2,0,time_values_opt,time_to_opt,start_time,psol);
  gl_best_value=best_value;
  sa_start_time=clock();
// Preparing SA parameters
  t_max=init_temperature(n1,n2,SAMPLE_SIZE,coef,seed1,pinst,psol);
  if (t_max==0) return gl_best_value;
  tred=(log(MIN_TEMPERATURE)-log(t_max))/log(COOLING_FACTOR);(psol+11)->perf=tred;
  tlength=TEMPERATURE_COEF*ntot;(psol+12)->perf=tlength;
  for (i=1;i<=nbig;i++) (psol+i)->ed=0;
  h1(0,0)=0;h2(0,0)=0;
// Iterations of SA_VNS ...
  while (outer_stop_cond==0)
     {st++;stop_cond=0;
      best_value=sa(n1,n2,tred,tlength,best_value,st,time_limit,
          time_limit_sa,t_max,share1,coef,&stop_cond,&gl_best_value,
          seed1,&seed2,time_values_opt,time_to_opt,start_time,pinst,psol);
      if (st==1)
         {(psol+7)->perf=best_value;
          if (stop_cond==0)
             {end=clock();
              sa_elapsed_time=(end-sa_start_time)/CLK_TCK;
// The number of iterations (restarts) is predicted
              pred_st_count=ceil(time_limit_sa/sa_elapsed_time);if (pred_st_count==0) pred_st_count=1;
             }
           else pred_st_count=1;
          time_limit_vns_restart=time_limit_vns/pred_st_count;
         }
       else if (best_value<(psol+7)->perf) (psol+7)->perf=best_value;
      if (stop_cond==1)
         {fprintf(out,"Too small CPU time limit for simulated annealing! \n");fflush(out);
          if (time_limit_vns<0.001) stop_cond=2;
         }
      if (stop_cond==2)
         {if (best_value<gl_best_value)
             {store_best_sol_glob(n1,n2,st,time_values_opt,time_to_opt,start_time,psol);
              gl_best_value=best_value;
             }
          (psol+3)->perf=st;
          return gl_best_value;
         } 
      if (time_limit_vns<0.001) continue;
      for (i=1;i<=n1;i++) {(psol+i)->p1=(psol+i)->bp1;(psol+(psol+i)->p1)->pi1=i;}
      for (i=1;i<=n2;i++) {(psol+i)->p2=(psol+i)->bp2;(psol+(psol+i)->p2)->pi2=i;}
      vns_start_time=clock();
// Initial execution of local search procedure
      sol_value=init_insertions(n1,n2,pinst,psol);
      if (sol_value!=best_value)
         {fprintf(out,
          "!!! some discrepancy in solution values:  SA_value=%12ld  calc_value=%12ld\n",best_value,sol_value);
          exit(1);
         }
      sol_value+=local_search_insertions(n1,n2,pinst,psol);
      ((psol+10)->perf)++;
      if (sol_value<best_value)
         {store_best_sol_vns(n1,n2,psol);
          best_value=sol_value;
          if (best_value<gl_best_value)
             {store_best_sol_glob(n1,n2,st,time_values_opt,time_to_opt,start_time,psol);
              gl_best_value=best_value;
             }
         }
      stop_cond=0;
      while (stop_cond==0)
         {dist=dist_min;
// Parameters used in the inner loop of variable neighborhood search
          if (dist_min>limit) dist_min=limit;
          fct=random(&seed3,coef)*(DISTANCE_MAX_HIGH-DISTANCE_MAX_LOW)+DISTANCE_MAX_LOW;
          db=fct*ntot;dist_max=db;
          if (dist_max>limit) dist_max=limit; else if (dist_max<dist_min) dist_max=dist_min;
          dist_step=db/DISTANCE_STEP_COEF;if (dist_step<1) dist_step=1;
          while (dist<=dist_max && stop_cond==0)
// Shaking
             {shake(n1,n2,dist,share1,share2,coef,&seed4,psol);
// Local search
              sol_value=init_insertions(n1,n2,pinst,psol);
              sol_value+=local_search_insertions(n1,n2,pinst,psol);
              ((psol+10)->perf)++;
// Applying neighborhood change rule
              if (sol_value<best_value)
                 {store_best_sol_vns(n1,n2,psol);
                  best_value=sol_value;
                  if (best_value<gl_best_value)
                     {store_best_sol_glob(n1,n2,st,time_values_opt,time_to_opt,start_time,psol);
                      gl_best_value=best_value;
                     }
                  dist=dist_min;
                 }
               else dist+=dist_step;
// Checking variable neighborhood search termination rule
              end=clock();
              elapsed_time=(end-start_time)/CLK_TCK;
              if (elapsed_time>=time_limit) stop_cond=2;
              if (stop_cond==0)
                 {vns_elapsed_time=(end-vns_start_time)/CLK_TCK;
                  if (vns_elapsed_time>=time_limit_vns_restart) stop_cond=1;
                 }
             } // inner while
         } // outer while
      if (stop_cond==2) outer_stop_cond=1;
     } // outmost while
// Storing the total number of SA-VNS starts executed
  (psol+3)->perf=st;
  return gl_best_value;
 }




// "bgcmp_sa_vns" is an SA and VNS hybrid for the bipartite graph crossing minimization problem (BGCMP)
 void bgcmp_sa_vns(char *in_file_name,char *out_file_name,double seed,
      long time_limit,Results *pres)
 {FILE *out,*in;
  Instance *pinst;
  Solution *psol;

  int i;
  int n1,n2,nbig,ec;
  int v1,v2;
  int dist_min;
  int time_values[5],time_values_opt[5];
  long best_value,sol_value;
  double time_limit_sa;
  double time_in_seconds,time_to_opt;
  double seed_saved;
  double doub;
  clock_t start_time;

  if ((in=fopen(in_file_name,"r"))==NULL)
     {printf("  fopen failed for input");exit(1);}
// n1 (respect., n2) is the number of vertices in the first (respect., second) part. 
// ec is the number of edges
  fscanf(in,"%d %d %d",&n1,&n2,&ec);
  if ((out=fopen(out_file_name,"w"))==NULL)
     {printf("  fopen failed for output  %s",out_file_name);exit(1);}
  seed_saved=seed;
  if (n1>=n2) nbig=n1; else nbig=n2;
// Allocation of core memory for an array of Instance structure objects 
  ALS(pinst,Instance,nbig+1)
  for (i=0;i<=n1;i++) ALI((pinst+i)->a1,n2+1)
  for (i=0;i<=n2;i++) ALI((pinst+i)->a2,n1+1)
  for (i=0;i<=nbig;i++) ALI((pinst+i)->z,nbig+1)
  for (i=0;i<=n1;i++) ALI((pinst+i)->h1,nbig+1)
  for (i=0;i<=n2;i++) ALI((pinst+i)->h2,nbig+1)
  for (i=0;i<=n1;i++) ALI((pinst+i)->b1,n1+1)
  for (i=0;i<=n2;i++) ALI((pinst+i)->b2,n2+1)
// Allocation of core memory for an array of Solution structure objects 
// that contains data used in various methods. 
// ST_COUNT is the length of performance statistics array 'perf'
  if (nbig>ST_COUNT) i=nbig; else i=ST_COUNT;
  ALS(psol,Solution,i+1)

// Input of the graph
  for (i=1;i<=n1;i++) (psol+i)->aln1=0;
  for (i=1;i<=n2;i++) (psol+i)->aln2=0;
  for (i=1;i<=ec;i++)
     {fscanf(in,"%d %d",&v1,&v2);
      v2-=n1;
      ((psol+v1)->aln1)++;a1(v1,(psol+v1)->aln1)=v2;
      ((psol+v2)->aln2)++;a2(v2,(psol+v2)->aln2)=v1;
     }

  dist_min=DISTANCE_MIN;
  doub=SA_TIME_QUOTA;
  if (doub>=1) time_limit_sa=time_limit; else time_limit_sa=time_limit*doub;
  start_time=clock();
// Running simulated annealing and variable neighborhood search applied iteratively 
  best_value=sa_vns(out,n1,n2,dist_min,time_limit,time_limit_sa,
      time_values_opt,&time_to_opt,&seed,start_time,pinst,psol);
  time_in_seconds=take_time(time_values,start_time);
  (psol+6)->perf=best_value;
  sol_value=get_best_value(n1,n2,pinst,psol);
  if (best_value != sol_value) fprintf(out,
      "!!! some discrepancy in solution values: %8ld   %8ld\n",best_value,sol_value);

// Printing parameters, data characteristics and some statistics 
// regarding the performance of the algorithm
  fprintf(out,"   parameters:                                 \n");
  fprintf(out,"      SA_TIME_QUOTA                          = %5lf\n",SA_TIME_QUOTA);
  fprintf(out,"      DISTANCE_MIN                           = %5d\n",DISTANCE_MIN);
  fprintf(out,"      DISTANCE_MAX_LOW                       = %6lf\n",DISTANCE_MAX_LOW);
  fprintf(out,"      DISTANCE_MAX_HIGH                      = %6lf\n",DISTANCE_MAX_HIGH);
  fprintf(out,"      DISTANCE_STEP_COEF                     = %5d\n",DISTANCE_STEP_COEF);
  fprintf(out,"      COOLING_FACTOR                         = %6lf\n",COOLING_FACTOR);
  fprintf(out,"      MIN_TEMPERATURE                        = %8lf\n",MIN_TEMPERATURE);
  fprintf(out,"      TEMPERATURE_COEF                      = %3d\n",TEMPERATURE_COEF);
  fprintf(out,"      SAMPLE_SIZE                            = %5d\n",SAMPLE_SIZE);
  fprintf(out,"      seed for random number generator       = %8lf\n",seed_saved);
  fprintf(out,"   number of vertices in the first part      = %5d\n",n1);
  fprintf(out,"   number of vertices in the second part     = %5d\n",n2);
  fprintf(out,"   number of edges                           = %5d\n",ec);
  fprintf(out,"   time limit                                = %5ld\n",time_limit);
  fprintf(out,"   total of LS starts                        = %4ld\n",(psol+10)->perf);
  fprintf(out,"   value found by SA                         = %8ld\n",(psol+7)->perf);
  fprintf(out,"   number of SA-VNS starts executed          = %4ld\n",(psol+3)->perf);
  fprintf(out,"   number of improvements during SA-VNS      = %4ld\n",(psol+1)->perf);
  fprintf(out,"   last improvement at SA-VNS start no.      = %3ld\n",(psol+2)->perf);
  fprintf(out,"   value found by SA-VNS                     = %8ld\n",(psol+6)->perf);
  (psol+4)->perf=(long)time_to_opt;
  fprintf(out,"   time to the best solution: %d : %d : %d.%3d  (=%4ld seconds)\n",
      time_values_opt[1],time_values_opt[2],time_values_opt[3],time_values_opt[4],(long)time_to_opt);
  (psol+5)->perf=(long)time_in_seconds;
  fprintf(out,"   total time of the algorithm: %d : %d : %d.%3d  (=%4ld seconds)\n",
      time_values[1],time_values[2],time_values[3],time_values[4],(long)time_in_seconds);
  fflush(out);

// Saving results
  if (pres!=NULL)
     {pres->value=best_value;
      pres->total_time=time_in_seconds;
      pres->time_to_best=time_to_opt;
      pres->characts[0]=n1;
      pres->characts[1]=n2;
      pres->characts[2]=ec;
      if (pres->sol1==NULL) {ALI(pres->sol1,n1+1) ALI(pres->sol2,n2+1)}
      for (i=1;i<=n1;i++) sol1(i)=(psol+i)->gbp1;
      for (i=1;i<=n2;i++) sol2(i)=(psol+i)->gbp2;
     }

// Releasing the memory
  if (psol!=NULL) free(psol);
  for (i=0;i<=n2;i++) if ((pinst+i)->b2!=NULL) free((pinst+i)->b2);
  for (i=0;i<=n1;i++) if ((pinst+i)->b1!=NULL) free((pinst+i)->b1);
  for (i=0;i<=n2;i++) if ((pinst+i)->h2!=NULL) free((pinst+i)->h2);
  for (i=0;i<=n1;i++) if ((pinst+i)->h1!=NULL) free((pinst+i)->h1);
  for (i=0;i<=nbig;i++) if ((pinst+i)->z!=NULL) free((pinst+i)->z);
  for (i=0;i<=n2;i++) if ((pinst+i)->a2!=NULL) free((pinst+i)->a2);
  for (i=0;i<=n1;i++) if ((pinst+i)->a1!=NULL) free((pinst+i)->a1);
  if (pinst!=NULL) free(pinst);
// Closing input and output files
  fclose(out);fclose(in);
 }
