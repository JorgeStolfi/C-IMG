/* Last edited on 2003-09-25 16:18:37 by stolfi */

#include <js.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char **argv)
{
  FILE *f;
  int i,j,z1,z2;
  float r;
  f=fopen("test-image.pgm","w");
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
  return 0;
}
