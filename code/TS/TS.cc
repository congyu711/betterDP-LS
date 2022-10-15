#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "cstring"
#define N 1000
#define S1 1
#define tenure 100
#define S2 5000
#define ELEMENT 10
#define S3 5000   // number of allowed no-improving generations 
#define D 20
#define SHUFFLEDIV 1
#define RANGE 1
#define S4 0.5
// #define CLOCKS_PER_SEC ((clock_t)1000)
#define TimePermission 60 * CLOCKS_PER_SEC
typedef struct NODE
{
    int NO;
    struct NODE *next;
} NODE;
int ordero[N][2];
int lseqo[N];
int rseqo[N];
int lcrosso[N][N];
int rcrosso[N][N];
int supremelseq[N];
int supremerseq[N];
int supremeorder[N][2];
int childlcross[N][N], childrcross[N][N], childbest, childorder[N][2], childlseq[N], childrseq[N];
NODE *graph[N];
int vertex;
int right;
int left;
int freeright;
int freeleft;
int aggr;
int watershed;
int shuffles;
int candidate[N];
clock_t start, finish;
void read(const char *target);
void startmatrix(int lseq[], int rseq[], int lcross[][N], int rcross[][N], int order[][2]);
void initialize(int lcross[][N], int rcross[][N], int lseq[], int rseq[], int order[][2], int freeleft, int freeright, int opt);
void changematrix(int cross[][N], int vertex1, int vertex2);
void generate(int lseq[], int rseq[], int lcross[][N], int rcross[][N], int lstart, int rstart, int lcrossgen[][N], int rcrossgen[][N]);
int crosscal(int lcross[][N], int leq[]);
void tabusearch_with_shuffle(int lseq[], int rseq[], int rcross[][N], int lcross[][N], int *best, int order[][2], int ok);
int order[ELEMENT][N][2];
int lseq[ELEMENT][N];
int rseq[ELEMENT][N];
int best[ELEMENT];
int lcross[ELEMENT][N][N];
int rcross[ELEMENT][N][N];
int lcrossgen[N][N];
int rcrossgen[N][N];
int tabu[N][N];
int bestorder[N][2];
int bestlseq[N], bestrseq[N];
int bestlcross[N][N], bestrcross[N][N];
int reckon(int lseq[], int rseq[], int order[][2]);
int main(int argc,char **argv)
{
    srand(time(0));
    start = clock();
    // char filename[100] = "/home/congyu/dpls/instances/GB_1_rnd1_01/GB_1_rnd1_01_0001_30.txt";
    char *filename=argv[1];
    read(filename);
    startmatrix(lseqo, rseqo, lcrosso, rcrosso, ordero);
    for (int k = 0; k < ELEMENT; k++)
    {
        for (int i = 0; i < vertex; i++)
        {
            lseq[k][i] = lseqo[i];
            rseq[k][i] = rseqo[i];
            for (int j = 0; j < 2; j++)
                order[k][i][j] = ordero[i][j];
            for (int j = 0; j < vertex; j++)
            {
                lcross[k][i][j] = lcrosso[i][j];
                rcross[k][i][j] = rcrosso[i][j];
            }
        }
    }
    for (int k = 0; k < ELEMENT; k++)
    {
        initialize(lcross[k], rcross[k], lseq[k], rseq[k], order[k], freeleft, freeright, 0);
        best[k] = crosscal(lcross[k], lseq[k]);
    }
    aggr = left - freeleft + vertex - freeright;
    shuffles = aggr / SHUFFLEDIV;
    watershed = left - freeleft;
    int x = 0;
    for (int i = 0; i < vertex; i++)
    {
        if (ordero[i][0] == 0)
        {
            candidate[x++] = i;
        }
    }
    int generation = 0;
    int minlastgen = best[0];
    int noimproving = 0;
    int adventgen = 0;
    int supremetime;
    int minslastgen;

    while (finish - start < TimePermission)
    {
        generation++;
        for (int k = 0; k < ELEMENT; k++)
        {
            if (generation == 1)
                tabusearch_with_shuffle(lseq[k], rseq[k], lcross[k], rcross[k], &(best[k]), order[k], 0);
            else
            {
                if (k == minslastgen)
                    tabusearch_with_shuffle(lseq[k], rseq[k], lcross[k], rcross[k], &(best[k]), order[k], 0);
                else
                    tabusearch_with_shuffle(lseq[k], rseq[k], lcross[k], rcross[k], &(best[k]), order[k], 1);
            }
        }
        int max = best[0];
        int min = best[0];
        int spot = 0;
        int spotmin = 0;
        for (int k = 1; k < ELEMENT; k++)
        {
            if (best[k] > max)
            {
                max = best[k];
            }
            else if (best[k] < min)
            {
                min = best[k];
                spotmin = k;
            }
        }
        minslastgen = spotmin;
        if (generation == 1)
            minlastgen = min;
        int seedl = rand() % ELEMENT;
        // int seedl=spot;
        int seedr;
        do
        {
            seedr = rand() % ELEMENT;
        } while (best[seedr] - best[seedl] > D || -best[seedr] + best[seedl] > D);
        // seedr=spotmin;
        for (int i = 0; i < left; i++)
            childlseq[i] = lseq[seedl][i];
        for (int i = 0; i < right; i++)
            childrseq[i] = rseq[seedr][i];
        // childrseq[i]=rseq[seedl][i]; //�ڶ��ֽ���
        for (int i = 0; i < left; i++)
            for (int j = 0; j < 2; j++)
                childorder[i][j] = order[seedl][i][j];
        for (int i = left; i < vertex; i++)
            for (int j = 0; j < 2; j++)
                childorder[i][j] = order[seedr][i][j];
        // childorder[i][j]=order[seedl][i][j];//�ڶ��ֽ���
        for (int i = 0; i < left; i++)
            for (int j = 0; j < left; j++)
                childlcross[i][j] = 0;
        for (int i = left; i < vertex; i++)
            for (int j = left; j < vertex; j++)
                childrcross[i][j] = 0;
        /* �ڶ��ֽ���
        int lordering[left-freeleft],rordering[vertex-freeright];
        int pointl=0,pointr=0;
        for(int i=0;i<left;i++)
        if(lseq[seedr][i]>=freeleft)
        lordering[pointl++]=lseq[seedr][i];
        for(int i=0;i<right;i++)
        if(rseq[seedr][i]>=freeright)
        rordering[pointr++]=rseq[seedr][i];
        pointl=0;pointr=0;
        for(int i=0;i<left;i++)
        if(childlseq[i]>=freeleft){
            childlseq[i]=lordering[pointl];
            childorder[lordering[pointl++]][1]=i;
        }
        for(int i=0;i<right;i++)
        if(childrseq[i]>=freeright){
            childrseq[i]=rordering[pointr];
            childorder[rordering[pointr++]][1]=i;
        }
        */
        startmatrix(childlseq, childrseq, childlcross, childrcross, childorder);
        childbest = crosscal(childlcross, childlseq);
        tabusearch_with_shuffle(childlseq, childrseq, childlcross, childrcross, &childbest, childorder, 0);
        // printf("\n%d",childbest);
        int differences[ELEMENT], lists[ELEMENT];
        for (int i = 0; i < ELEMENT; i++)
        {
            for (int j = 0; j < vertex; j++)
            {
                if (childorder[j][1] != order[i][j][1])
                    differences[i]++;
            }
        }
        int minimum = vertex;
        int maximum = 0;
        for (int i = 0; i < ELEMENT; i++)
        {
            if (differences[i] < minimum)
                minimum = differences[i];
            if (differences[i] > maximum)
                maximum = differences[i];
        }
        int theta[ELEMENT];
        int currentmin = 1;
        int timeshitting;
        for (int i = 0; i < ELEMENT; i++)
        {
            theta[i] = S4 * ((differences[i] - minimum) / maximum - minimum + 1) + (1 - S4) * ((max - best[i]) / max - min + 1);
            if (theta[i] < currentmin)
            {
                spot = i;
                currentmin = theta[i];
                timeshitting = 2;
            }
            else if (theta[i] == currentmin)
            {
                if (rand() % timeshitting == 0)
                {
                    spot = i;
                    timeshitting++;
                }
            }
        }
        for (int i = 0; i < left; i++)
            lseq[spot][i] = childlseq[i];
        for (int i = 0; i < right; i++)
            rseq[spot][i] = childrseq[i];
        for (int i = 0; i < left; i++)
            for (int j = 0; j < 2; j++)
                order[spot][i][j] = childorder[i][j];
        for (int i = left; i < vertex; i++)
            for (int j = 0; j < 2; j++)
                order[spot][i][j] = childorder[i][j];
        for (int i = 0; i < left; i++)
        {
            for (int j = 0; j < left; j++)
                lcross[spot][i][j] = childlcross[i][j];
        }
        for (int i = left; i < vertex; i++)
        {
            for (int j = left; j < vertex; j++)
                rcross[spot][i][j] = childrcross[i][j];
        }
        best[spot] = childbest;
        if (childbest < min)
        {
            min = childbest;
            spotmin = spot;
            minslastgen = spot;
        }
        finish = clock();
        if (min >= minlastgen && generation != 1)
        {
            noimproving++;
            if (noimproving >= S3)
                break;
        }
        else
        {
            minlastgen = min;
            noimproving = 0;
            adventgen = generation;
            for (int i = 0; i < left; i++)
                supremelseq[i] = lseq[spotmin][i];
            for (int i = 0; i < right; i++)
                supremerseq[i] = rseq[spotmin][i];
            for (int i = 0; i < left; i++)
                for (int j = 0; j < 2; j++)
                    supremeorder[i][j] = order[spotmin][i][j];
            for (int i = left; i < vertex; i++)
                for (int j = 0; j < 2; j++)
                    supremeorder[i][j] = order[spotmin][i][j];
            supremetime = finish - start;
        }
        if (generation % 10 == 0)
            printf("\n%d generation,best so far:%d", generation, minlastgen);
    }
    printf("\nbest possible:%d\nat:%d generation %fseconds", minlastgen, adventgen, (float)supremetime / (float)CLOCKS_PER_SEC);
    printf("\nlseq:");
    for (int i = 0; i < left; i++)
        printf("%d ", supremelseq[i]);
    printf("\nrseq:");
    for (int i = 0; i < right; i++)
        printf("%d ", supremerseq[i]-left);
    // printf("\n%d!",reckon(supremelseq,supremerseq,supremeorder));
    FILE *fp;
    fp = fopen("Resultr6.txt", "a");
    fprintf(fp, "S4 %f\n%s:\nbest possible:%d\nat:%d generation %fseconds with RANGE %d", S4, filename, minlastgen, adventgen, (float)supremetime / (float)CLOCKS_PER_SEC, RANGE);
    fprintf(fp, "\nlseq:");
    for (int i = 0; i < left; i++)
        fprintf(fp, "%d ", supremelseq[i]);
    fprintf(fp, "\nrseq:");
    for (int i = 0; i < right; i++)
        fprintf(fp, "%d ", supremerseq[i]);
    fprintf(fp, "\n");
    fclose(fp);
    exit(1);
}

void read(const char *target)
{
    FILE *fp = NULL;
    int mark = 0;
    if ((fp = fopen(target, "r")) != NULL)
    {
        int temp;
        fscanf(fp, "%d", &temp);
        fscanf(fp, "%d", &left);
        fscanf(fp, "%d", &right);
        vertex = left + right;
        int v = -1;
        int loop = 0;
        for (int i = 0; i < vertex; i++)
        {
            graph[i] = (NODE *)malloc(sizeof(NODE));
            graph[i]->next = NULL;
        }
        while (fscanf(fp, "%d", &temp) != EOF)
        {
            loop++;
            if (temp == 0 || temp == 1)
            {
                v++;
                if (temp == 1)
                {
                    if (mark == 1)
                        mark = 0;
                    int x;
                    fscanf(fp, "%d", &x);
                    ordero[v][0] = 1;
                    ordero[v][1] = x;
                }
                else
                {
                    if (mark == 0)
                    {
                        if (v < left)
                            freeleft = v;
                        else
                            freeright = v;
                        mark = 1;
                    }

                    int x;
                    fscanf(fp, "%d", &x);
                    ordero[v][0] = 0;
                    if (v < left)
                        ordero[v][1] = x;
                    else
                        ordero[v][1] = x - left;
                }
            }
            else
            {
                NODE *nd = graph[v];
                NODE *nd1 = graph[temp];
                while (nd->next != NULL)
                    nd = nd->next;
                nd->next = (NODE *)malloc(sizeof(NODE));
                nd->next->NO = temp;
                nd->next->next = NULL;
                while (nd1->next != NULL)
                    nd1 = nd1->next;
                nd1->next = (NODE *)malloc(sizeof(NODE));
                nd1->next->NO = v;
                nd1->next->next = NULL;
            }
        }
        for (int i = 0; i < left; i++)
        {
            if (ordero[i][0] == 1)
            {
                int note = ordero[i][1];
                lseqo[note] = i;
            }
            else
                lseqo[i] = i;
        }
        for (int i = left; i < vertex; i++)
        {
            if (ordero[i][0] == 1)
            {
                int note = ordero[i][1];
                rseqo[note] = i;
            }
            else
                rseqo[i - left] = i;
        }
    }
    else
        printf("open failed!\n");
}

void startmatrix(int lseq[], int rseq[], int lcross[][N], int rcross[][N], int order[][2])
{

    for (int i = 0; i < left; i++)
    {
        for (int j = 0; j < i; j++)
        {
            NODE *P = graph[i];
            while (P->next != NULL)
            {
                int r1 = P->next->NO;
                NODE *Q = graph[j];
                while (Q->next != NULL)
                {
                    int r2 = Q->next->NO;
                    if (order[r1][1] > order[r2][1])
                        lcross[i][j]++;
                    else if (order[r1][1] < order[r2][1])
                        lcross[j][i]++;
                    Q = Q->next;
                }
                P = P->next;
            }
        }
    }
    for (int i = left; i < vertex; i++)
    {
        for (int j = left; j < i; j++)
        {
            NODE *P = graph[i];
            while (P->next != NULL)
            {
                int l1 = P->next->NO;
                NODE *Q = graph[j];
                while (Q->next != NULL)
                {
                    int l2 = Q->next->NO;
                    if (order[l1][1] > order[l2][1])
                        rcross[i][j]++;
                    else if (order[l1][1] < order[l2][1])
                        rcross[j][i]++;
                    Q = Q->next;
                }
                P = P->next;
            }
        }
    }
}

void changematrix(int cross[][N], int vertex1, int vertex2)
{
    NODE *p = graph[vertex1];
    int i = 0;
    while (p->next != NULL)
    {
        NODE *q = graph[vertex2];
        while (q->next != NULL)
        {
            int nei1 = p->next->NO;
            int nei2 = q->next->NO;
            cross[nei1][nei2]++;
            cross[nei2][nei1]--;
            q = q->next;
        }
        p = p->next;
    }
}

void generate(int lseq[], int rseq[], int lcross[][N], int rcross[][N], int lstart, int rstart, int lcrossgen[][N], int rcrossgen[][N])
{
    for (int i = 0; i < left - lstart; i++)
        for (int j = 0; j < left; j++)
            lcrossgen[i][j] = 0;
    for (int i = 0; i < vertex - rstart; i++)
        for (int j = left; j < vertex; j++)
            rcrossgen[i][j] = 0;

    for (int i = lstart; i < left; i++)
    {
        int j = 1;
        for (int k = 0; k < lstart; k++)
            lcrossgen[i - lstart][0] += lcross[lseq[i]][lseq[k]];
        while (j != lstart + 1)
        {
            lcrossgen[i - lstart][j] += lcrossgen[lseq[i - lstart]][j - 1];
            lcrossgen[i - lstart][j] -= lcross[lseq[i]][lseq[j - 1]];
            lcrossgen[i - lstart][j] += lcross[lseq[j - 1]][lseq[i]];
            j++;
        }
        while (j != left)
        {
            lcrossgen[i - lstart][j] = lcrossgen[i - lstart][lstart];
            j++;
        }
    }
    for (int i = rstart; i < vertex; i++)
    {
        int j = 1 + left;
        for (int k = 0; k < rstart - left; k++)
            rcrossgen[i - rstart][left] += rcross[rseq[i]][rseq[k]];
        while (j != rstart + 1)
        {
            rcrossgen[i - rstart][j] += rcrossgen[rseq[i - rstart]][j - 1];
            rcrossgen[i - rstart][j] -= rcross[rseq[i]][rseq[j - left - 1]];
            rcrossgen[i - rstart][j] += rcross[rseq[j - left - 1]][rseq[i]];
            j++;
        }
        while (j != vertex)
        {
            rcrossgen[i - rstart][j] = rcrossgen[i - rstart][rstart];
            j++;
        }
    }
}

void initialize(int lcross[][N], int rcross[][N], int lseq[], int rseq[], int order[][2], int freeleft, int freeright, int opt)
{
    int lremainder = left - freeleft;
    int rremainder = vertex - freeright;
    int remainders = lremainder + rremainder;
    int x = 0;
    while (remainders != 0)
    {
        generate(lseq, rseq, lcross, rcross, left - lremainder, vertex - rremainder, lcrossgen, rcrossgen);
        /*for(int i=0;i<lremainder;i++){
            for(int j=0;j<left;j++)
            printf("%d ",lcrossgen[i][j]);
            printf("/\n");
        }
        printf("\n");
        for(int i=0;i<rremainder;i++){
            for(int j=left;j<vertex;j++)
            printf("%d ",rcrossgen[i][j]);
            printf("/\n");
        }*/
        int lminseq[left];
        int rminseq[vertex];
        if (opt == 0)
        {
            for (int i = 0; i < lremainder; i++)
            {
                lminseq[i] = lcrossgen[i][0];
                for (int j = 1; j <= left - lremainder; j++)
                {
                    if (lcrossgen[i][j] < lminseq[i])
                        lminseq[i] = lcrossgen[i][j];
                }
            }
            for (int i = 0; i < rremainder; i++)
            {
                rminseq[i] = rcrossgen[i][left];
                for (int j = left + 1; j <= vertex - rremainder; j++)
                {
                    if (rcrossgen[i][j] < rminseq[i])
                        rminseq[i] = rcrossgen[i][j];
                }
            }
        }
        else if (opt == 1)
        {
            for (int i = 0; i < lremainder; i++)
            {
                lminseq[i] = lcrossgen[i][0];
                int chooose = rand() % (left - lremainder + 1);
                lminseq[i] = lcrossgen[i][chooose];
            }
            for (int i = 0; i < rremainder; i++)
            {
                rminseq[i] = rcrossgen[i][left];
                int chooose = rand() % (left - lremainder + 1 - left) + left;
                rminseq[i] = rcrossgen[i][chooose];
            }
        }
        int lmin = lminseq[0], lmax = lminseq[0];
        int rmin = rminseq[0], rmax = rminseq[0];
        for (int i = 0; i < lremainder; i++)
        {
            if (lminseq[i] < lmin)
                lmin = lminseq[i];
            if (lminseq[i] > lmax)
                lmax = lminseq[i];
        }
        for (int i = 0; i < rremainder; i++)
        {
            if (rminseq[i] < rmin)
                rmin = rminseq[i];
            if (rminseq[i] > rmax)
                rmax = rminseq[i];
        }
        int ltao = lmin + S1 * (lmax - lmin), rtao = rmin + S1 * (rmax - rmin);
        // printf("go:%d %d %d %d\n",lmax,lmin,rmax,rmin);
        int lrestricted[lremainder], rrestricted[rremainder];
        int lreco = 0, rreco = 0;
        for (int i = 0; i < lremainder; i++)
        {
            if (lminseq[i] <= ltao)
                lrestricted[lreco++] = i;
        }
        for (int i = 0; i < rremainder; i++)
        {
            if (rminseq[i] <= rtao)
                rrestricted[rreco++] = i;
        }
        int reco = lreco + rreco;
        int choose = rand() % reco;
        if (choose < lreco)
        {
            int moving = lrestricted[choose] + left - lremainder;
            int pmoving[left - lremainder];
            x = 0;
            for (int i = 0; i <= left - lremainder; i++)
            {
                if (lcrossgen[moving - left + lremainder][i] == lminseq[moving - left + lremainder])
                    pmoving[x++] = i;
            }
            int position = pmoving[rand() % x];
            int element = lseq[moving];
            // printf("%d:%d->%d\n",lrestricted[choose],moving+1,position+1);
            for (int i = moving; i > position; i--)
                changematrix(rcross, lseq[i - 1], lseq[moving]);
            for (int i = moving - 1; i >= position; i--)
            {
                lseq[i + 1] = lseq[i];
                order[lseq[i + 1]][1] = i + 1;
            }
            lseq[position] = element;
            order[lseq[position]][1] = position;
            lremainder--;
            remainders--;
        }
        else
        {
            choose -= lreco;
            int moving = rrestricted[choose] + vertex - rremainder - left;
            int pmoving[vertex - rremainder];
            x = 0;
            for (int i = left; i <= vertex - rremainder; i++)
            {
                if (rcrossgen[moving - vertex + rremainder + left][i] == rminseq[moving - vertex + rremainder + left])
                    pmoving[x++] = i - left;
            }
            int position = pmoving[rand() % x];
            int element = rseq[moving];
            for (int i = moving; i > position; i--)
            {
                changematrix(lcross, rseq[i - 1], rseq[moving]);
            }
            // printf("%d:%d->%d\n",rrestricted[choose],moving+left+1,position+1+left);
            for (int i = moving - 1; i >= position; i--)
            {
                rseq[i + 1] = rseq[i];
                order[rseq[i + 1]][1] = i + 1;
            }
            rseq[position] = element;
            order[rseq[position]][1] = position;
            rremainder--;
            remainders--;
        }
        // int crosses=crosscal();
    }
}

int crosscal(int lcross[][N], int lseq[])
{
    int crosses = 0;
    for (int i = 0; i < left; i++)
    {
        for (int j = 0; j < i; j++)
            crosses += lcross[lseq[j]][lseq[i]];
    }
    return crosses;
}

void tabusearch_with_shuffle(int lseq[], int rseq[], int lcross[][N], int rcross[][N], int *best, int order[][2], int ok)
{
    int round = 0;
    int count = 0;
    int times = 0;
    for (int i = 0; i < aggr; i++)
        for (int j = 0; j < vertex; j++)
            tabu[i][j] = 0;
    int historybest = *best + 1;
    int crosses = *best;
    int markling;
    while (1)
    {
        int dvall = -crosses, moveall, pointall, place;
        int dvtabu = -crosses, movetabu, pointabu, tabucount = 1;
        int dvntabu = -crosses, moventabu, pointnabu, ntabucount = 1;
        for (int i = 0; i < aggr; i++)
        {
            int mark, pos, dvcurrent = -crosses, movecurrentup = 0, movecurrentdown = 0, movecurrent, poscurrent;
            mark = candidate[i];
            pos = order[mark][1];
            int upjudge[RANGE + 1], downjudge[RANGE + 1];
            upjudge[0] = 0;
            downjudge[0] = 0;
            if (i < watershed)
            {
                for (int i = 1; i <= RANGE; i++)
                {
                    if (pos - i >= 0)
                        upjudge[i] = upjudge[i - 1] + lcross[lseq[pos - i]][lseq[pos]] - lcross[lseq[pos]][lseq[pos - i]];
                    if (pos + i < left)
                        downjudge[i] = downjudge[i - 1] + lcross[lseq[pos]][lseq[pos + i]] - lcross[lseq[pos + i]][lseq[pos]];
                    if (pos - i >= 0)
                        if (upjudge[i] > upjudge[i - 1])
                            movecurrentup = i;
                    if (pos + i < left)
                        if (downjudge[i] > downjudge[i - 1])
                            movecurrentdown = i;
                }
                int judge = (upjudge[movecurrentup] > downjudge[movecurrentdown]) ? 1 : 0;
                dvcurrent = (judge == 1) ? (upjudge[movecurrentup]) : (downjudge[movecurrentdown]);
                movecurrent = (judge == 1) ? (-movecurrentup) : movecurrentdown;
                poscurrent = (judge == 1) ? pos - movecurrentup : pos + movecurrentdown;
            }
            else
            {
                for (int i = 1; i <= RANGE; i++)
                {
                    if (pos - i >= 0)
                        upjudge[i] = upjudge[i - 1] + rcross[rseq[pos - i]][rseq[pos]] - rcross[rseq[pos]][rseq[pos - i]];
                    if (pos + i < right)
                        downjudge[i] = downjudge[i - 1] + rcross[rseq[pos]][rseq[pos + i]] - rcross[rseq[pos + i]][rseq[pos]];
                    if (pos - i >= 0)
                        if (upjudge[i] > upjudge[i - 1])
                            movecurrentup = i;
                    if (pos + i < right)
                        if (downjudge[i] > downjudge[i - 1])
                            movecurrentdown = i;
                }
                int judge = (upjudge[movecurrentup] > downjudge[movecurrentdown]) ? 1 : 0;
                dvcurrent = (judge == 1) ? (upjudge[movecurrentup]) : (downjudge[movecurrentdown]);
                movecurrent = (judge == 1) ? (-movecurrentup) : movecurrentdown;
                poscurrent = (judge == 1) ? pos - movecurrentup : pos + movecurrentdown;
            }

            if (round < tabu[i][poscurrent])
            {
                if (dvcurrent > dvtabu)
                {
                    dvtabu = dvcurrent;
                    movetabu = movecurrent;
                    pointabu = mark;
                    tabucount = 2;
                }
                else if (dvcurrent == dvtabu)
                {
                    if (rand() % tabucount == 0)
                    {
                        movetabu = movecurrent;
                        pointabu = mark;
                    }
                    tabucount++;
                }
                else
                    continue;
            }
            else
            {
                if (dvcurrent > dvntabu)
                {
                    dvntabu = dvcurrent;
                    moventabu = movecurrent;
                    pointnabu = mark;
                    ntabucount = 2;
                }
                else if (dvcurrent == dvntabu)
                {
                    if (rand() % tabucount == 0)
                    {
                        moventabu = movecurrent;
                        pointnabu = mark;
                    }
                    ntabucount++;
                }
                else
                    continue;
            }
        }
        if (crosses - dvtabu < historybest && dvtabu > dvntabu)
        {
            dvall = dvtabu;
            moveall = movetabu;
            pointall = pointabu;
        }
        else if (dvntabu != -crosses)
        {
            dvall = dvntabu;
            moveall = moventabu;
            pointall = pointnabu;
        }
        else
        {
            dvall = dvtabu;
            moveall = movetabu;
            pointall = pointabu;
        }
        crosses -= dvall;
        place = order[pointall][1];
        int oldplace = place;
        int newplace = place + moveall;
        if (pointall < left)
        {
            if (moveall < 0)
            {
                for (int i = oldplace; i > newplace; i--)
                {
                    int element = lseq[i - 1];
                    lseq[i] = element;
                    order[element][1]++;
                    changematrix(rcross, lseq[i], pointall);
                }
                lseq[newplace] = pointall;
                order[pointall][1] = newplace;
            }
            else
            {
                for (int i = oldplace; i < newplace; i++)
                {
                    int element = lseq[i + 1];
                    lseq[i] = element;
                    order[element][1]--;
                    changematrix(rcross, pointall, lseq[i]);
                }
                lseq[newplace] = pointall;
                order[pointall][1] = newplace;
            }
            int hint = pointall - freeleft;
            tabu[hint][newplace] = tenure + round;
        }
        else
        {
            if (moveall < 0)
            {
                for (int i = oldplace; i > newplace; i--)
                {
                    int element = rseq[i - 1];
                    rseq[i] = element;
                    order[element][1]++;
                    changematrix(lcross, rseq[i], pointall);
                }
                rseq[newplace] = pointall;
                order[pointall][1] = newplace;
            }
            else
            {
                for (int i = oldplace; i < newplace; i++)
                {
                    int element = rseq[i + 1];
                    rseq[i] = element;
                    order[element][1]--;
                    changematrix(lcross, pointall, rseq[i]);
                }
                rseq[newplace] = pointall;
                order[pointall][1] = newplace;
            }
            int hint = pointall - freeright + watershed;
            tabu[hint][newplace] = tenure + round;
        }
        if (crosses < historybest)
        { // markling=round;
            historybest = crosses;
            for (int i = 0; i < vertex; i++)
                for (int j = 0; j < 2; j++)
                    bestorder[i][j] = order[i][j];
            for (int i = 0; i < left; i++)
                bestlseq[i] = lseq[i];
            for (int i = 0; i < right; i++)
                bestrseq[i] = rseq[i];
            for (int i = 0; i < left; i++)
            {
                for (int j = 0; j < left; j++)
                    bestlcross[i][j] = lcross[i][j];
            }
            for (int i = left; i < vertex; i++)
            {
                for (int j = left; j < vertex; j++)
                    bestrcross[i][j] = rcross[i][j];
            }
            count = 0;
        }
        round++;
        count++;
        if (count == S2)
        {
            // printf("\n%d",markling);
            if (times == ok)
                break;
            if (shuffles != 0)
            {
                for (int i = 0; i < left; i++)
                {
                    lseq[i] = bestlseq[i];
                }
                for (int i = 0; i < right; i++)
                {
                    rseq[i] = bestrseq[i];
                }
                for (int i = 0; i < vertex; i++)
                    for (int j = 0; j < 2; j++)
                        order[i][j] = bestorder[i][j];
                for (int i = 0; i < left; i++)
                    for (int j = 0; j < left; j++)
                        lcross[i][j] = bestlcross[i][j];
                for (int i = left; i < vertex; i++)
                    for (int j = left; j < vertex; j++)
                        rcross[i][j] = bestrcross[i][j];
                int shuffled[aggr];
                int markingl = left, markingr = vertex;
                for (int i = 0; i < aggr; i++)
                    shuffled[i] = 1;
                for (int i = 0; i < shuffles; i++)
                {
                    int selection = rand() % aggr;
                    do
                    {
                        selection = rand() % aggr;
                    } while (shuffled[selection] == 0);
                    shuffled[selection] = 0;
                    int targetp = candidate[selection];
                    int targets = order[targetp][1];
                    if (selection < watershed)
                    {
                        for (int i = targets; i < left - 1; i++)
                        {
                            int element = lseq[i + 1];
                            lseq[i] = element;
                            order[element][1]--;
                            changematrix(rcross, targetp, lseq[i]);
                        }
                        lseq[left - 1] = targetp;
                        order[targetp][1] = left - 1;
                        markingl--;
                    }
                    else
                    {
                        for (int i = targets; i < right - 1; i++)
                        {
                            int element = rseq[i + 1];
                            rseq[i] = element;
                            order[element][1]--;
                            changematrix(lcross, targetp, rseq[i]);
                        }
                        rseq[right - 1] = targetp;
                        order[targetp][1] = right - 1;
                        markingr--;
                    }
                }
                initialize(lcross, rcross, lseq, rseq, order, markingl, markingr, 0);
                crosses = crosscal(lcross, lseq);
                count = 0;
                times++;
            }
            else
                break;
        }
        // if(round%1000==0){
        //	printf("\n%drounds:%d(%d:%d->%d)",round,crosses,pointall,oldplace,newplace);
        //	printf("\n%d",crosscal(lcross,lseq));
        //	}
    }
    // printf("\n%drounds:%d",round,historybest);
    // printf("\nlseq:");
    for (int i = 0; i < left; i++)
    {
        // printf("%d ",bestlseq[i]);
        lseq[i] = bestlseq[i];
    }
    // printf("\nrseq:");
    for (int i = 0; i < right; i++)
    {
        //	printf("%d ",bestrseq[i]);
        rseq[i] = bestrseq[i];
    }
    for (int i = 0; i < vertex; i++)
        for (int j = 0; j < 2; j++)
            order[i][j] = bestorder[i][j];
    for (int i = 0; i < left; i++)
        for (int j = 0; j < left; j++)
            lcross[i][j] = bestlcross[i][j];
    for (int i = left; i < vertex; i++)
        for (int j = left; j < vertex; j++)
            rcross[i][j] = bestrcross[i][j];
    *best = historybest;
}

int reckon(int lseq[], int rseq[], int order[][2])
{
    int crossed = 0;
    for (int i = 0; i < left; i++)
    {
        int pointf = lseq[i];
        for (int j = i + 1; j < left; j++)
        {
            NODE *p = graph[pointf];
            int pointr = lseq[j];
            while (p->next != NULL)
            {
                NODE *q = graph[pointr];
                while (q->next != NULL)
                {
                    if (order[p->next->NO][1] > order[q->next->NO][1])
                        crossed++;
                    q = q->next;
                }
                p = p->next;
            }
        }
    }
    return crossed;
}
