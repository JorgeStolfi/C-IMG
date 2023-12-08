
void Trabalheador(int m,int n,cmp *x,cmp *y,plr *whereto)
{
int i,j;
plr plrx,plry,plrz;


for(i=0;i<m;i++)
for(j=0;j<n;j++)
{
plrx=mCmp(*(x+i*n+j));
plry=mCmp(*(y+i*n+j));
plrz.r=plry.r/plrx.r;
plrz.an=plry.an-plrx.an;
*(whereto+i*n+j)=plrz;
}
  
}

int Comparator(int m,int n,cmp *x,cmp *y,float *qr,float *qa)
{
int i,j;
plr *trab;
int res;
plr plrz;
float quadroan,quadror;

trab=(plr*)calloc(m*n,sizeof(plr));
Trabalheador(m,n,x,y,trab);
res=0;
quadroan=0.0;
quadror=0.0;

for(i=0;i<m;i++)
for(j=0;j<n;j++)
{
plrz=*(trab+i*n+j);
if((plrz.r>0.0 & plrz.r>1.05) | 
   (plrz.r<0.0 & plrz.r<0.95) | 
   (plrz.an<0.0 & plrz.an<-0.05) | 
   (plrz.an>0.0 & plrz.an>0.05 )) res++;
quadroan=quadroan+sqrt(plrz.an*plrz.an);
quadror=quadror+sqrt((plrz.r-1.0)*(plrz.r-1.0));

}

printf("\n quadro %f  %f and %i dismatches",quadror,quadroan,res);
free(trab);
 
*qa=quadroan; *qr=quadror;
return(res);
}

