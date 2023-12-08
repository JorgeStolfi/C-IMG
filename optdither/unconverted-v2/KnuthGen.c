#include <stdio.h>

void FillKnuth(int *x,int n,int size)
{
int i,j,az;

if(n==1) *x=0; 
  else
{
for(i=0;i<n/2;i++)
for(j=0;j<n/2;j++) 
{
az=(*(x+i*size+j))*4;
*(x+i*size+j)=az;
*(x+i*size+j+n/2)=az+2;
*(x+(i+n/2)*size+j)=az+3;
*(x+(i+n/2)*size+n/2+j)=az+1;
};

};
}

void main()
{
int *knuth;
char *name;
FILE *f;
int i,j,n;

name=(char*)calloc(20,sizeof(char));
printf("\n Hallo, going to built Knuth...");
printf("\n Enter size... n=");
scanf("%i",&n);
printf("\n Wellcome, what file would you like to get it in?\n");
scanf("%s",name);
f=fopen(name,"w");
free(name);

knuth=(int*)calloc(n*n,sizeof(int));

for(i=1;i<=n;i=i*2)
FillKnuth(knuth,i,n);

fprintf(f,"%i",n);
for(i=0;i<n;i++)
for(j=0;j<n;j++)
fprintf(f,"\n%i",*(knuth+i*n+j));

fclose(f);
free(knuth);
}

