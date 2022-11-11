#include <alloc.h>
#include <process.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef struct
     {int *sol1;            /* permutation of the first part vertices  */
      int *sol2;            /* permutation of the second part vertices */
      long value;           /* solution value                          */
      double total_time;    /* total time, secs                        */
      double time_to_best;  /* time to the best solution, secs         */
      int characts[10];     /* some characteristics:                   */
                            /*   characts[0] - size of the first part  */
                            /*   characts[1] - size of the second part */
                            /*   characts[2] - number of edges         */
     }Results;

void bgcmp_sa_vns(char *,char *,double,long,Results *);


 void main(int argc,char **argv)
 {FILE *out,*out2;
  Results *pres;
  char in_file_name[80];
  char out_file_name[80];
  char out_file_name_init[80];
  char summary_file_name[80];
  char result_file_name[80];
  char instance_name[80];
  char bkv_string[80];
  char param_string[80];
  int i,j,k;
  int count;
  long time_limit;
  long min_value=0,best_known_value=0;
  long dif_best;
  double seeds [21]={0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 
      8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 
      17000, 18000, 19000, 20000};
  char numbs[31][2]={ {'0'}, {'1'}, {'2'}, {'3'}, {'4'}, {'5'},
      {'6'}, {'7'}, {'8'}, {'9'}, {'1','0'}, {'1','1'}, {'1','2'}, 
      {'1','3'}, {'1','4'}, {'1','5'}, {'1','6'}, {'1','7'}, 
      {'1','8'}, {'1','9'}, {'2','0'} };
  char res_ext[4]={'.', 'r', 'e', 's'};
  double av_value=0.,av_time=0.,av_time_tot=0.;
  double dif_aver;

  if (argc<=4) {printf("  specify data, output, summary and results files");exit(1);}
  strcpy(in_file_name,argv[1]);
  strcpy(out_file_name_init,argv[2]);
  strcpy(summary_file_name,argv[3]);
  strcpy(result_file_name,argv[4]);
  strcpy(param_string,argv[5]);time_limit=atol(param_string);
  strcpy(param_string,argv[6]);count=atoi(param_string);
  if (argc==8)
     {strcpy(bkv_string,argv[7]);best_known_value=atol(bkv_string);}
  pres=(Results *)calloc(1,sizeof(Results));
  if ((out=fopen(summary_file_name,"w"))==NULL)
     {printf("  fopen failed for output  %s",summary_file_name);exit(1);}
  if ((out2=fopen(result_file_name,"a"))==NULL)
     {printf("  fopen failed for output  %s",result_file_name);exit(1);}
  for (i=0;i<strlen(in_file_name);i++)
     {if (in_file_name[i]=='.') break;
      instance_name[i]=in_file_name[i];
     }
  instance_name[i]='\0';
  for (i=1;i<=count;i++)
     {strcpy(out_file_name,out_file_name_init);
      k=strlen(out_file_name);
      out_file_name[k]='_';k++;
      if (i<10) out_file_name[k]=numbs[i][0];
       else {out_file_name[k]=numbs[i][0];k++;out_file_name[k]=numbs[i][1];}
      for (j=0;j<=3;j++) out_file_name[k+j+1]=res_ext[j];
      out_file_name[k+j+1]='\0';

      bgcmp_sa_vns(in_file_name,out_file_name,seeds[i],time_limit,pres);
      fprintf(out," %8ld       %8lf       %8lf\n",pres->value,
          pres->time_to_best,pres->total_time);
      av_value+=pres->value;
      av_time+=pres->time_to_best;av_time_tot+=pres->total_time;
      if (i==1 || pres->value<min_value) min_value=pres->value;
     }
  av_value/=count;av_time/=count;av_time_tot/=count;
  fprintf(out,"----------------------------------------\n");
  fprintf(out,"%11.3lf    %11.3lf    %11.3lf\n",av_value,av_time,av_time_tot);
  fprintf(out,"minimum value over all runs = %8ld \n",min_value);
  fprintf(out2,"%s   %8ld   %11.3lf   %11.3lf",instance_name,min_value,av_value,av_time);
  if (argc==8)
     {dif_best=min_value-best_known_value;dif_aver=av_value-best_known_value;
      fprintf(out,"best known value = %8ld \n",best_known_value);
      fprintf(out,"dif_best = %8ld  dif_aver = %8lf\n",dif_best,dif_aver);
      fprintf(out2,"  %8ld  %8lf\n",dif_best,dif_aver);
     }
   else fprintf(out2,"\n");
  free(pres);
  fclose(out);fclose(out2);
 }
