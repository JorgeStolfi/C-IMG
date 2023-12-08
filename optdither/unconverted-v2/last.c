#include "/tmp_mnt/home/spec/glasov/prog/dot/os.c"

void main()
{
OS x,y;
int m,n,krit,limit;
cmp *fltmsk;
FILE *f;
int intz;
int i;
cmp cmpz;
char cz;

m=8; n=8;

printf("\n\n   *** WELLCOME TO CIRCUS ***\n\n   Limit ? ");
scanf("%i",&limit);

bore(m,n,&x); bore(m,n,&y);

fltmsk=(cmp*)calloc(m*n,sizeof(cmp));
ReadMask("/tmp_mnt/home/spec/glasov/dat/diam-low8.dat" ,fltmsk);
printf("\n Mask read");

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

/* Initialization for the first step */
FFT2d(1,m,n,x.obj,x.fur); i=0;

/* Main loop. */
do {
printf("\n   STEP %i",i);

OneStep(&x,&y,fltmsk);

SaveCT("dat/0Xfur",m,n,x.fur,'p',7);
SaveCT("dat/0Xreb",m,n,x.reb,'p',7);
SaveCT("dat/0Yreb",m,n,y.reb,'p',7);
SaveCT("dat/0Yobj",m,n,y.obj,'p',7);
SaveCT("dat/0Yfur",m,n,y.fur,'p',7);   

krit=Comparator(m,n,x.fur,y.fur);

troca(&x,&y); i++;

/*
SaveCT("dat/1Xfur",m,n,x.fur,'p',7);
SaveCT("dat/1Xreb",m,n,x.reb,'p',7);
SaveCT("dat/1Yreb",m,n,y.reb,'p',7);
SaveCT("dat/1Yobj",m,n,y.obj,'p',7);
SaveCT("dat/1Yfur",m,n,y.fur,'p',7);   
scanf("%c",&cz);
scanf("%c",&cz); */

}  while(krit>0 & i<limit);

printf("\n\n   *** The result is ... ***\n");
PrT(m,n,x.T);
printf("\n\n Does it look like ? ");

SaveT("result",m,n,x.T);

f=fopen("resultM.dat","w");
fprintf(f,"%i %i",m,n);
for(i=0;i<m*n;i++)
fprintf(f,"\n%i",*(x.T+i));
fclose(f);
}
