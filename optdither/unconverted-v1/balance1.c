#include "/tmp_mnt/home/spec/glasov/ready/TOTAL.c"

void balance(int m,int n,cmp *x)
{
int i,j;
cmp cmpz;


for(i=0;i<m;i++)
for(j=i;j<n;j++)
{
cmpz=A(*(x+i*n+j),*(x+i+j*n));
cmpz.re=cmpz.re/2.0;
cmpz.im=cmpz.im/2.0;
*(x+i*n+j)=cmpz;
*(x+i+j*n)=cmpz;
}
}

void main()
{
int i,j;
cmp x[8][8];

for(i=0;i<8;i++)
for(j=0;j<8;j++)
{
x[i][j].re=(10*i+j)*1.0;
x[i][j].im=(i+10*j)*1.0;
}

SaveCT("before",8,8,x,'c',5);

balance(8,8,x);

SaveCT("after",8,8,x,'c',5);

}
