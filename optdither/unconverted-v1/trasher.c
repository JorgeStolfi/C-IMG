#include <stdio.h>
#include <math.h>
void main()
{
FILE *f;
int i,j,z1,z2;
float r;
f=fopen("trash.pgm","w");
fprintf(f,"P2 256 256 128");
for(i=0;i<256;i++)
for(j=0;j<256;j++)
{
z1=(j-128)*(j-128);
z2=(i-128)*(i-128);
r=(float)(z1+z2);
fprintf(f,"\n%i",(int)(63*(1+sin(r/250.0))));
}
fclose(f);

}
