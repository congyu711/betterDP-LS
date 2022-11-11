#include <stdio.h>

// Variable neighborhood search
#define DISTANCE_MIN                3
#define DISTANCE_MAX_LOW         0.05
#define DISTANCE_MAX_HIGH         0.2
#define DISTANCE_STEP_COEF         25

// Simulated annealing
#define COOLING_FACTOR           0.95
#define MIN_TEMPERATURE        0.0001
#define TEMPERATURE_COEF          100
#define SAMPLE_SIZE              5000
#define TRIAL_LIMIT               100
#define SA_TIME_QUOTA             0.5

#define ST_COUNT                   20

#define ALS(X,Y,Z) if ((X=(Y *)calloc(Z,sizeof(Y)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}
#define ALI(X,Z) if ((X=(int *)calloc(Z,sizeof(int)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}

#define a1(X,Y) *((pinst+X)->a1+Y)
#define a2(X,Y) *((pinst+X)->a2+Y)
#define z(X,Y) *((pinst+X)->z+Y)
#define h1(X,Y) *((pinst+X)->h1+Y)
#define h2(X,Y) *((pinst+X)->h2+Y)
#define b1(X,Y) *((pinst+X)->b1+Y)
#define b2(X,Y) *((pinst+X)->b2+Y)
#define sol1(Y) *(pres->sol1+Y)
#define sol2(Y) *(pres->sol2+Y)


typedef struct
     {int *a1;        /*  */
      int *a2;        /*  */
      int *z;         /*  */
      int *h1;        /*  */
      int *h2;        /*  */
      int *b1;        /*  */
      int *b2;        /*  */
     }Instance;

typedef struct
     {int gbp1;       /*  */
      int gbp2;       /*  */
      int bp1;        /*  */
      int bp2;        /*  */
      int p1;         /*  */
      int p2;         /*  */
      int pi1;        /*  */
      int pi2;        /*  */
      int aln1;       /*  */
      int aln2;       /*  */
      int a;          /*  */
      int av;         /*  */
      int zln;        /*  */
      int ad;         /*  */
      int ed;         /*  */
      int s1;         /*  */
      int s2;         /*  */
      long perf;      /*  */
     }Solution;

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
