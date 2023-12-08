void dither(int m,int n,int mm,int nn,cmp *x,float mx,cmp *y,plr *z,int *matr)
{
int i,j,ii,jj;
float dd;
plr plrz;
cmp zero;

zero.re=0.0; zero.im=0.0;
dd=mx/mm/nn;

ii=0; jj=0;
for(i=0;i<m;i++)
{
for(j=0;j<n;j++)
{
if(ii==mm) ii=0;
if(jj==nn) jj=0;
plrz=*(z+i*n+j);
if(plrz.r>dd*(*(matr+ii*nn+jj)+0.5) ) *(y+i*n+j)=*(x+i*n+j);
 else *(y+i*n+j)=zero;
jj++;
};
ii++;
}

}
