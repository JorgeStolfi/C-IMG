/*         Table of content        
FFT
FFT2d
Lanco
LiveAgua 
Pgmer
PrT 
ReadM
ReadMask
Show 
SaveCT
SaveM
SaveT */

#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <strings.h>
typedef struct { float re; float im; } cmp ;
typedef struct { float r ; float an; } plr ;
typedef struct { int m; int n; float* a; } pgmDscr;
typedef struct { int m; int n; int *num; } MDscr;


cmp A(cmp x,cmp y)
{
cmp res;
res.re=x.re+y.re;
res.im=x.im+y.im;
return (res);
}
cmp M(cmp x,cmp y)
{
cmp res;
res.re=x.re*y.re-x.im*y.im;
res.im=x.im*y.re+x.re*y.im;
return(res);
}
cmp mPolar(plr x)
{
cmp res;
res.re=cos(x.an)*x.r;
res.im=sin(x.an)*x.r;
return(res);
}
float mAbs(cmp x)
{
return(sqrt(x.re*x.re+x.im*x.im));
}
plr mExp_p(cmp x)
{
plr res;
res.r=exp(x.re);
res.an=x.im;
return(res);
}

plr mCmp(cmp x)
{
plr res;
res.r=mAbs(x);
res.an=(float)atan(((double)x.im)/x.re);
return(res);
}

float mCmpdd(int m,int n,cmp *wfrom,plr *wto)
{
int i,j;
float maxampl;
plr z;

maxampl=0.0;

for(i=0;i<m;i++)
for(j=0;j<n;j++)
{
z=mCmp(*(wfrom+i*n+j));
if(z.r>maxampl) maxampl=z.r;
*(wto+i*n+j)=z;
}

return(maxampl);
}

float Power(int m,int n,cmp *x)
{
int i,j;
float power;
plr z;

power=0.0;

for (i=0;i<m;i++)
for (j=0;j<n;j++)
{
z=mCmp(*(x+i*n+j));
power=power+z.r;
};

return(power);
}


/* reprint is used by SaveCT */

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
strcat(rname,".CAdat");
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

/* Saves cmp array as a table in CT.dat file */

void SaveCT(char *name,int m,int n,cmp *x,char mode,int length)
/* modes : 'p' - savein polar form, 'c' - save in complex form */
{
char *rname;
cmp cmpz;
plr plrz;
int i;
FILE *f;

rname=calloc(32,sizeof(char));
strcat(rname,name);
strcat(rname,".CAdat");
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


/* Records plr array to 2 pgm-files : ..A & ..F - ampl & phase,   */
/* dividing every member by factor maxampl. Phase range +-pi/2    */

void SShow(char *name,int m,int n,plr *ax,float maxampl)
{
FILE *f1,*f2;
char *nname;
plr z;
int i,j;
float trace1;
int trace2;


nname=calloc(20,sizeof(char));

strcat(nname,name);
strcat(nname,"A.pgm");
f1=fopen(nname,"w");
free(nname);

nname=calloc(20,sizeof(char));
strcat(nname,name);
strcat(nname,"F.pgm");
f2=fopen(nname,"w");

fprintf(f1,"P2 %i %i 128",m,n);
fprintf(f2,"P2 %i %i 128",m,n);
free(nname);

for (i=0;i<m;i++)
for (j=0;j<n;j++)
{
z=*(ax+i*n+j);
trace1=z.r*126.0/maxampl;
trace2=(int)trace1;
fprintf(f1,"\n%i",trace2);
trace1=63.0+z.an*40.0;
trace2=(int)trace1;
fprintf(f2,"\n%i",trace2);
}

fclose(f1);
fclose(f2);
}


/****** Records cmp array to 2 pgm-files ******************************/

void Show(char *name,int m,int n,cmp *x)
{
plr *y;
float mx;
/* cmp *extd;
int mm,nn,llm,lln;

llm=log2i(m); lln=log2i(n); */

y=(plr*)calloc(m*n,sizeof(plr));

mx=mCmpdd(m,n,x,y);

SShow(name,m,n,y,mx);
}


/****** Prints int array to stdout ********************************************/

void PrT(int m,int n,int *x)
{
int i,j;

for (i=0;i<m;i++)
{
printf("\n");
for (j=0;j<n;j++)
printf(" %.2i",*(x+i*n+j));
}

}


/***** Saves matrix to dat-file ***********************************************/

void SaveT(char *name,int m,int n,int *x)
{
int i,j;
FILE *f;
char *rname;

rname=calloc(20,sizeof(char));
strcat(rname,name);
strcat(rname,"T.dat");
f=fopen(rname,"w");

for (i=0;i<m;i++)
{
fprintf(f,"\n");
for (j=0;j<n;j++)
fprintf(f," %.2i",*(x+i*n+j));
};

fclose(f);
free(rname);

}

/* SaveM saves int matrix to M.dat file */

void SaveM(char *name,int m,int n,int *x)
{
int i,j;
FILE *f;
char *rname;

rname=calloc(20,sizeof(char));
strcat(rname,name);
strcat(rname,"M.dat");
f=fopen(rname,"w");

fprintf(f,"%i %i",m,n);

for (i=0;i<m*n;i++)
fprintf(f,"\n%.3i",*(x+i*n+j));

fclose(f);
free(rname);

}
 
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
free(rname);

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

printf("\n pgmRead message... Tudo bem.");
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

dd=1.0/mm/nn;

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

/* ReadM reads matrix from M.dat-file and returns M-descriptor */

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
free(rname);

fscanf(f,"%i %i",&m,&n);
res.m=m;
res.n=n;

res.num=(int*)calloc(m*n,sizeof(int));

for(i=0;i<m*n;i++)
fscanf(f,"%i",(res.num+i));

fclose(f);
printf("\n ReadM message... Tudo bem. %i %i %i *",res.m,res.n,*(res.num));
return(res);
}
/* rebuild Restores usual frequency order ***************************************/

void rebuild(int m,int n,cmp *x,cmp *y)
{
int i,j;

for(i=0;i<m/2;i++)
for(j=0;j<n/2;j++)
{
*(y+i*n+j)=*(x+(i+m/2)*n+j+n/2);
*(y+(i+m/2)*n+j+n/2)=*(x+i*n+j);
*(y+(i+m/2)*n+j)=*(x+i*n+j+n/2);
*(y+i*n+j+n/2)=*(x+(i+m/2)*n+j);
}
}


/***** Do not use it *******************************************************/

void Filter(int m,int n,cmp *ax,cmp *ay,int par)
{
int i,j,krit1,krit2;
cmp zero;

zero.re=0.0;zero.im=0.0;
printf("\n Par= %i Zero: re %f im %f",par,zero.re,zero.im);

for(i=0;i<m;i++)
for(j=0;j<n;j++)
{
krit1=n*(i-m/2)+m*(j-n/2);
krit2=n*(i-m/2)-m*(j-n/2);
if((krit1>par) | (krit1<-par) | (krit2>par) | (krit2<-par))
*(ay+i*n+j)=zero; else *(ay+i*n+j)=*(ax+i*n+j);
}

}


/***** Restores a dot matrix according to groth of values ******************/

void LiveAgua(int m,int n,cmp *x, int *y)
{
int i,j,k,im,jm;
float mx;
cmp z;
int *oo;

oo=(int*)calloc(m*n,sizeof(int));
for(i=0;i<m;i++)
for(j=0;j<n;j++)
*(oo+i*n+j)=0;

for(k=m*n;k>0;k--)
{

mx=-20000.0;

for(i=0;i<m;i++)
for(j=0;j<n;j++)
{
z=*(x+i*n+j);

if (*(oo+i*n+j)==0)
  if(z.re>=mx) 
  {
  mx=z.re;
  im=i;jm=j;
  };
};

*(y+im*n+jm)=k-1;
*(oo+im*n+jm)=1000;
};

}

/***** Reads manually typed pattern for MskFlt ******************************/

void ReadMask(char *name,cmp *ay)
{
int i,j,m,n;
cmp z;
float az;
FILE *f;

f=fopen(name,"r");

fscanf(f,"%i %i",&m,&n);

for(i=0;i<m;i++)
for(j=0;j<n;j++)
{
fscanf(f,"%f",&az);
z.re=az;
z.im=0.0;
*(ay+i*n+j)=z;
}

fclose(f);
}


/***** Filts putting 0 where 0 are in the pattern ***************************/

void MskFlt(int m,int n,cmp *ax,cmp *ay,cmp *az)
{
int i,j;


for(i=0;i<m;i++)
for(j=0;j<n;j++)
*(ay+i*n+j)=M(*(ax+i*n+j),*(az+i*n+j));

}


/***** Multiplies every member after flt to restore total energy ************/

void Justifier(int m,int n,cmp *x,cmp *y,cmp *res)
{
float jst;
int i,j;
cmp z;

jst=Power(m,n,x)/Power(m,n,y);

for(i=0;i<m;i++)
for(j=0;j<n;j++)
{
z=*(y+i*n+j);

z.re=z.re*jst;
z.im=z.im*jst;

*(res+i*n+j)=z;
}

}
int Rev(int n,int x)
{
int i; /* counter */
int m; /* mask    */
int res,cur,shift;
m=1;
res=0;
for(i=1;i<=n;i++)
   { cur=x&m;
     shift=n+1-2*i;
     if ( shift>0 ) cur=cur<<shift; else cur=cur>>(-shift);
     res=res+cur;
       m=m<<1;
    } ;
return(res);
}

int log2i(int x)
{
int res;
for(res=0;x>1;res++) x=x/2;
return(res);
};

void RevSort( int n, cmp *(wherefrom),cmp *(whereto))
{
int i,k;
for (i=0;i<n;i++)
{
k=Rev(log2i(n),i);
whereto[i]=wherefrom[k];
};
}
cmp WWW(int way,int n,int k)
{
cmp res;
float angle;
angle=way*k*6.28/n;
res.re=cos(angle);
res.im=sin(angle);
return(res);
}

void Lanco(int way,int n,cmp *wherefrom,cmp *whereto)
{
int i,j;
for(j=0;j<2;j++)
for(i=0;i<n;i++)
whereto[n*j+i]=A(wherefrom[i],M(WWW(way,2*n,n*j+i),wherefrom[n+i]));
}


void FFT(int way,int n,cmp *ax,cmp *ay)
{
int i,j;
cmp *az,*as;
as=ay;
RevSort(n,ax,ay);

for(i=1;i<n;i=i*2)
{
az=ax;
ax=ay;
ay=az;

for(j=0;j<n/i/2;j++)
Lanco(way,i,ax+j*i*2,ay+j*i*2);
}
for(i=0;i<n;i++) as[i]=ay[i];
}

void TR(int n,int m,cmp *ax)
{
int i,j;
cmp z;
for(i=0;i<m;i++)
for(j=i+1;j<n;j++)
{
z=*(ax+i*n+j);
*(ax+i*n+j)=*(ax+j*m+i);
*(ax+j*m+i)=z;
}
}

void TR1(int m,int n,cmp *ax,cmp *ay)
{
int i,j;
for(i=0;i<m;i++)
for(j=0;j<n;j++)
*(ay+i*n+j)=*(ax+j*m+i);
}

void FFT2d(int way,int m,int n,cmp *ax,cmp *ay)
{
int i,j;
cmp *as;
cmp z,z1;

as=(cmp*)calloc(m*n,sizeof(cmp));
for(i=0;i<m*n;i++)*(as+i)=*(ax+i);

for(i=0;i<m;i++) FFT(way,n,as+i*n,ay+i*n);
TR1(m,n,ay,as);
for(i=0;i<n;i++) FFT(way,m,as+i*m,ay+i*m);
TR(n,m,ay);
if(way<0)
{
z.re=1.0/m/n;
z.im=0.0; 
for(i=0;i<m;i++)
for(j=0;j<n;j++)
{
z1=*(ay+i*n+j);
*(ay+i*n+j)=M(z,z1);
}
}
}



void ReadMaskInt(char *name,int *ay)
{
int i,j,m,n;
int z;
FILE *f;

f=fopen(name,"r");

fscanf(f,"%i %i",&m,&n);

for(i=0;i<m;i++)
for(j=0;j<n;j++)
{
fscanf(f,"%i",&z);
*(ay+i*n+j)=z;
}

fclose(f);
}
