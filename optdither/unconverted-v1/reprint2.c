
#include "/tmp_mnt/home/spec/glasov/ready/TOTAL.c"

void reprint(char *name,int length)
{ 
char *rname;
int i,jm,jn,m,n,l,k,trig;
char cz;
int ileft[64],iright[64];
char *left[64],*right[64];
FILE *f1,*f2;

/****** Openning files, allocating memory. ******/ 
rname=calloc(32,sizeof(char));
strcat(rname,name);
strcat(rname,".dat");
f1=fopen(rname,"r");
free(rname);
rname=calloc(32,sizeof(char));
strcat(rname,name);
strcat(rname,"CT.dat");
f2=fopen(rname,"w");
fprintf(f2,"\n\n   *** THIS IS %s ***\n",name);
free(rname);


/****** Reading dimention of the array ******/
fscanf(f1,"%i %i\n",&m,&n);

for(i=0;i<n*2;i++)
{
left[i]=calloc(8,sizeof(char));
right[i]=calloc(8,sizeof(char));
}

/****** Main loop ******/
for(jm=0;jm<m;jm++) /* Loop N0 */
{
fprintf(f2,"\n");
for(i=0;i<n*(length+2)+3;i++)
fprintf(f2,"-");
fprintf(f2,"\n%.2i>",jm); 

/* Loop N 1 */
for(jn=0;jn<n*2;jn++)
{
trig=0;
fscanf(f1,"%c",&cz);

/****** Inner loop, N1, until CR is reached ******/
ileft[jn]=0; iright[jn]=0;
for(i=0;cz!='\n';i++)
{

if(trig==0) { *(left[jn]+ileft[jn])=cz; 
               ileft[jn]++; }
       else { *(right[jn]+iright[jn])=cz; iright[jn]++; };
if(cz=='.') trig=1; 
fscanf(f1,"%c",&cz);
}
/* End of inner loop */

*(right[jn]+iright[jn])='\0';
*(left[jn]+ileft[jn])='\0';

}
/* End of N1 loop */

/****************************************************************/

/* Loop N2 */
for(jn=0;jn<2*n;jn=jn+2)
{

fprintf(f2,"| ");

/****** Printing digits before . & dot itself ******/
for(k=0;k<ileft[jn];k++)
fprintf(f2,"%c",*(left[jn]+k)); 


/****** length-k digits left to print ******/

/****** Printing digits left ******/
for(i=0;i<length-k & *(right[jn]+i)!='\0';i++) 
fprintf(f2,"%c",*(right[jn]+i));

/****** Putting spaces if neccessary ******/
for(l=0;l<length-k-i;l++)
fprintf(f2," ");

}
/* End of loop N 2 */

fprintf(f2,"| ");

fprintf(f2,"\n   ");
/* Loop N3 */
for(jn=1;jn<2*n;jn=jn+2)
{

fprintf(f2,"| ");

/****** Printing digits before . & dot itself ******/
for(k=0;k<ileft[jn];k++)
fprintf(f2,"%c",*(left[jn]+k)); 

/****** length-k-1 digits left to print ******/

/****** Printing digits left ******/
for(i=0;i<length-k & *(right[jn]+i)!='\0';i++) 
fprintf(f2,"%c",*(right[jn]+i));

/****** Putting spaces if neccessary ******/
for(l=0;l<length-k-i;l++)
fprintf(f2," ");

}
/* End of Loop N 3 */

fprintf(f2,"| ");

}
/****** End of main loop ******/

fprintf(f2,"\n");
for(i=0;i<n*(length+2)+3;i++)
fprintf(f2,"-");

fclose(f1); fclose(f2);
free(left); free(right);

}

void SaveCT(char *name,int m,int n,cmp *x,char mode,int length)
/* modes : 'p' - savein polar form, 'c' - save in complex form */
{
char *rname;
cmp cmpz;
plr plrz;
int i;
FILE *f;
char cz;

rname=calloc(32,sizeof(char));
strcat(rname,name);
strcat(rname,".dat");
f=fopen(rname,"w");
free(rname);

fprintf(f,"%i %i\n",m,n);

for(i=0;i<m*n;i++)
{
cmpz=*(x+i);
if (mode=='c')
fprintf(f,"%f\n%f\n",cmpz.re,cmpz.im);
else { plrz=mCmp(cmpz);
fprintf(f,"%f\n%f\n",plrz.r,plrz.an); };
} /* End for */       

fclose(f);

reprint(name,length);

}

void main()
{
cmp x[8][8];
int i,j;
cmp cmpz;

for (i=0;i<8;i++)
for (j=0;j<8;j++)
{
cmpz.re=10.0*sin(i/2.5);
cmpz.im=-6.66-10*j;
x[i][j]=cmpz;
}
printf("\n Array created...");

SaveCT("tr",8,8,x,'c',7);
printf("\n 1 Save passed");
SaveCT("ttrr",8,8,x,'p',6);

}
