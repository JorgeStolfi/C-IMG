#include "/tmp_mnt/home/spec/glasov/prog/dot/os.c"

void main()
{
OS x,y;
int i,ii,m,n,limit,krit;
cmp *fltmsk;
FILE *f;
cmp cmpz;
int intz; 
cmp *cmpp;
cmp *prev;
cmp *cur;
char cz;
int *intp;

printf("\n Let's look for it...");
printf("\n how long?");
scanf("%i",&limit);

m=8;n=8;

fltmsk=(cmp*)calloc(m*n,sizeof(cmp));
ReadMask("/tmp_mnt/home/spec/glasov/dat/diam-low.dat" ,fltmsk);

prev=(cmp*)calloc(m*n,sizeof(cmp));
cur=(cmp*)calloc(m*n,sizeof(cmp));

x.obj=(cmp*)calloc(m*n,sizeof(cmp));
x.fur=(cmp*)calloc(m*n,sizeof(cmp));
y.obj=(cmp*)calloc(m*n,sizeof(cmp));
y.fur=(cmp*)calloc(m*n,sizeof(cmp));
x.T=(int*)calloc(m*n,sizeof(int));
y.T=(int*)calloc(m*n,sizeof(int));
x.m=m; x.n=n;
y.m=m; y.n=n;
x.number=0;

/* Putting knuth's matrix to x */
f=fopen("/tmp_mnt/home/spec/glasov/dat/knuth8.dat","r");
fscanf(f,"%i",&intz);
for(i=0;i<m*n;i++)
{
fscanf(f,"%i",&intz);
cmpz.re=intz*1.0;
cmpz.im=0.0; 
*(x.obj+i)=cmpz;
*(x.T+i)=intz;
}
fclose(f);
printf("\n To start with - x.T");
PrT(m,n,x.T);

printf("\n *** POYEHALY... "); 

i=0; 

/************* ITERATION BODY *****************/

for(krit=1;krit>0 & i<limit;i++)
{
OneStep(x,y,fltmsk,cur);
/*for(ii=0;ii<m*n;ii++) */
/* printf("\n %.2f %.2f %.2f %.2f",(prev+ii)->re,(prev+ii)->im,
(y.fur+ii)->re,(y.fur+ii)->im); */

if(i>1) krit=Comparator(m,n,prev,cur);
PrT(m,n,y.T);
printf("\n Iteration N %i ...Pronto! TRAB = %i ",i,krit);
y.number=i;

cmpp=x.obj;  
x.obj=y.obj;
y.obj=cmpp;

cmpp=cur;  
cur=prev;
prev=cmpp;

intp=x.T;
x.T=y.T;
y.T=intp;
}

printf("\n Before...");
PrT(m,n,y.T);
printf("\n \n *** The result is \n");
PrT(m,n,x.T);

SaveT("result",m,n,y.T);

f=fopen("resultM.dat","w");
fprintf(f,"%i %i",m,n);
for(i=0;i<m*n;i++)
fprintf(f,"\n%i",*(y.T+i));
fclose(f);

free(fltmsk);
free(x.fur); free(y.fur);
free(x.obj); free(y.obj);

}
