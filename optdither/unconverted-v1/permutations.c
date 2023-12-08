#include "/tmp_mnt/home/spec/glasov/ready/TOTAL.c"

void analyse(int m,int n,int *x,int *y,float *crit,
int *fltmsk,int *trig)
{
cmp *obj,*fur;
int i;
float energy; plr plrz;

obj=(cmp*)calloc(m*n,sizeof(cmp));
fur=(cmp*)calloc(m*n,sizeof(cmp));

for(i=0;i<m*n;i++)
{
(obj+i)->re=*(x+i)*1.0;
(obj+i)->im=0.0;
}

FFT2d(1,m,n,obj,fur);

energy=0.0;
for(i=0;i<m*n;i++)
{
plrz=mCmp(*(fur+i));
energy=energy+*(fltmsk+i)*plrz.r;
}

if(energy<*crit)
{
*crit=energy;
*trig=*trig+1;

/* for(i=0;i<m*n;i++)
*(y+i)=*(x+i); */

SaveT("SUPER",m,n,x);
}

printf("*%i",*trig);

free(obj); free(fur);

}

void change(int i,int j,int *x,int *beginning)
{
int z;

z=*(x+i);
*(x+i)=*(x+j);
*(x+j)=z;

}

void circulation(int *what,int tail,int *beginning)
{
int i;char cz;
/*printf("\n Circulation call tail %i",tail);*/
for(i=0;i<tail-1;i++)
change(i,i+1,what+i,beginning);
/*scanf("%c",&cz);*/
}

void permutation(int *what,int tail,int m,int n,
int *counter1,int *counter2,int *beginning,
int *where,int *fltmsk,float *crit,int *trig)
{
int i;char cz;

if(tail==1) 
{ 
if(*counter1==9999) {*counter2=*counter2+1; *counter1=0;}; 
*counter1=*counter1+1; 
printf("\n%i",*counter1);
analyse(m,n,beginning,where,crit,fltmsk,trig);
}

else        { /* printf("\n Permutation call tail=%i",tail);*/
               for(i=1;i<tail;i++)
{
permutation(what,tail-1,m,n,counter1,counter2,beginning,
where,fltmsk,crit,trig);
circulation(what,tail,beginning);
}
}
}

void main()
{
int i,counter1,counter2,trig,m,n;
float crit;
int *x,*result,*fltmsk;
char cz;
FILE *f;
int intz;

counter1=0; counter2=0; trig=0; crit=10000.0; m=4; n=4;

/* Putting knuth's matrix to start with */

x=(int*)calloc(m*n,sizeof(int));
result=(int*)calloc(m*n,sizeof(int));

f=fopen("/tmp_mnt/home/spec/glasov/dat/knuth4.dat","r");
fscanf(f,"%i",&intz);
for(i=0;i<m*n;i++)
{ 
fscanf(f,"%i",&intz);
*(x+i)=intz; 
}
fclose(f);

PrT(m,n,x);

fltmsk=(int*)calloc(m*n,sizeof(int));
ReadMaskInt("/tmp_mnt/home/spec/glasov/dat/for_perm.dat" ,fltmsk);

analyse(m,n,x,result,&crit,fltmsk,&trig);
trig=0;
printf("\n KNUTH analyzed criterium=%f",crit);
scanf("%c",&cz);
crit=crit-1.0;

printf("\n POYEHALY !!!");

permutation(x,m*n,m,n,&counter1,&counter2,x,result,fltmsk,
&crit,&trig);

printf("\n PRONTO");


}

               
