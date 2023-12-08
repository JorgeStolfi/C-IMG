 typedef struct { int m; int n; float* a; } pgmDscr; 
 typedef struct { int m; int n; int *num; } MDscr;   


/***** Reads pgm-file and returns pgm-descriptor ******************************/

pgmDscr Readpgm(char *name)
{
int i,j;
int m,n,mx,z;
char c;
FILE* f;
char *rname;
pgmDscr res;

rname=calloc(20,sizeof(char));
strcat(rname,name);
strcat(rname,".pgm");
f=fopen(rname,"r");

for(i=0;i<3;i++)
fscanf(f,"%c",&c);

fscanf(f,"%i %i %i",&n,&m,&mx);

res.a=(float*)calloc(m*n,sizeof(float));
res.m=m;
res.n=n;

for(i=0;i<m;i++)
for(j=0;j<n;j++)
{
fscanf(f,"%i",&z);
*(res.a+i*n+j)=1.0*z/mx;
}

return(res);
}


/***** Makes a pbm-file from array x ****************************************/

void pbmer(char *name, int m,int n,int *x)
{
int i,j;
char *nname;
FILE *f;

nname=calloc(20,sizeof(char));
strcat(nname,name);
strcat(nname,".pbm");
f=fopen(nname,"w");

fprintf(f,"P1 %i %i",n,m);

for(i=0;i<m;i++)
for(j=0;j<n;j++)
fprintf(f,"\n%i",*(x+i*n+j));

fclose(f);
free(nname);

}



/***** Dithers array x using matrix matr, output to array y *****************/

void Dith(int m,int n,float *x,int *y,int mm,int nn,int *matr)
{
int i,j,ii,jj;
float dd;

dd=1.0/m/n;

ii=0; jj=0;
for(i=0;i<m;i++)
{
for(j=0;j<n;j++)
{
if(ii==mm) ii=0;
if(jj==nn) jj=0;
if( *(x+i*n+j)>dd*(*(matr+ii*nn+jj)+0.5) ) *(y+i*n+j)=1; else *(y+i*n+j)=0;
jj++;
};
ii++;
}
}

MDscr ReadM(char *name)
{
char *rname;
FILE *f;
MDscr res;
int i,m,n;

rname=calloc(20,sizeof(char));
strcat(rname,name);
strcat(rname,"M.dat");
f=fopen(rname,"r");

fscanf(f,"%i %i",&m,&n);
res.m=m;
res.n=n;

for(i=0;i<m*n;i++)
fscanf(f,"%i",(res.num+i));

fclose(f);
free(rname);

return(res);
}
