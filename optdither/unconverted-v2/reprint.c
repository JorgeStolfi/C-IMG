#include <stdio.h>
#include <malloc.h>
#include <strings.h>
#include <math.h>

void reprint(char *name,int length)
{ 
char *rname;
int i,jm,jn,m,n,l,ileft,iright,k,jj,trig;
char cz;
char *czp,*left,*right;
FILE *f1,*f2;

/****** Openning files, allocating memory. ******/ 
rname=calloc(32,sizeof(char));
strcat(rname,name);
strcat(rname,"CA.dat");
f1=fopen(rname,"r");
free(rname);
rname=calloc(32,sizeof(char));
strcat(rname,name);
strcat(rname,"CT.dat");
f2=fopen(rname,"w");
free(rname);

left=calloc(32,sizeof(char));
right=calloc(32,sizeof(char));

/****** Reading dimention of the array ******/
fscanf(f1,"%i %i\n",&m,&n);

/****** jj is used to CR when line is completed ******/
jj=0;

/****** Main loop ******/
for(jm=0;jm<m;j++)
{
fprintf(f2,"\n%.2i>",jm)
for(jn=0;jn<n;jn++)
{

if(jj==2*n) jj=0;

trig=0;
fscanf(f1,"%c",&cz);

/****** Inner loop, until CR is reached ******/
ileft=0; iright=0;
for(i=0;cz!='\n';i++)
{
if(trig==0) { *(left+i)=cz; ileft++; }
       else { *(right+iright)=cz; iright++; };
if(cz=='.') trig=1; 
fscanf(f1,"%c",&cz);
}

*(right+iright)='\0';
*(left+ileft)='\0';

if(jj==0) fprintf(f2,"\n%.2i*",j);
fprintf(f2," ");

/****** Printing digits before . & dot itself ******/
for(k=0;k<ileft;k++)
fprintf(f2,"%c",*(left+k)); 

/****** length-k-1 digits left to print ******/

/****** Printing digits left ******/
for(i=0;i<length-k & *(right+i)!='\0';i++) 
fprintf(f2,"%c",*(right+i));

/****** Putting spaces if neccessary ******/
for(l=0;l<length-k-i;l++)
fprintf(f2," ");

jj++;

}
}
/****** End of main loop ******/

fclose(f1); fclose(f2);
free(left); free(right);

}

void main()
{
int i;
FILE *f;

f=fopen("trCA.dat","w");

fprintf(f,"4 4\n");
for(i=0;i<16;i++)
{
fprintf(f,"%f\n",100.0*sin(i*1.0));
fprintf(f,"-6.66\n");
}

fclose(f);

reprint("tr",7);


}
