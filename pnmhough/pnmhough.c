#define PROG_NAME "pnmhough"
#define PROG_DESC "compute the Hough transform of an image"
#define PROG_VERS "0.0"

/* Last edited on 2021-07-17 23:44:49 by jstolfi*/

/* Copyright © 2003 by the State University of Campinas (UNICAMP).
** See the copyright, authorship, and warranty notice at end of file.
*/

#define PROG_HELP \
  PROG_NAME " [ INFILE ] [ OUTPREFIX ] [ OPTION]\n" \
  "where OPTION is:\n" \
  "1 - Copy input image to output file;\n" \
  "2 - Convert PPM to PGM;\n" \
  "3 - Draw a line on the image;\n" \
  "4 - Draw the lines recognized automatically (does ot work).\n" \
  "5 - Creates PGM file with local detector output (Not good, try option 6).\n" \
  "6 - Creates PGM file with weighted local detector output.\n" \
  "7 - Detects an approximate line crossing with 1/2 pixel accuracy.\n" \
  "8 - Paints lines automatically recognized (under test).\n" \
  "9 - Ditto with the detector of option 6.\n" \
  "H - Help"

#include <affirm.h>
#include <imagem.h>

#define N_LINHAS_SAIDA 40

int main(int argc, char **argv)
{
  int *R=NULL, *G=NULL, *B=NULL, *D=NULL, *E=NULL;
  int largura, altura, maxcor;
  int i=0;/*contador*/
  int *m;
  double *mm;
  char *opcao;
  int ld=50,cd=40;
  double mi,mj;
  char *arq_in, *pref_out;

  int nlinhas = 0;

  if (argc!=4)
    {
      printf("%s version %s, usage:\n", PROG_NAME, PROG_VERS);
      printf("%s\n", PROG_HELP);
      return 0;
    }
  
  arq_in = argv[1];
  pref_out = argv[2];
  opcao = argv[3];
  
  switch(opcao[0])
    {
      
    case '1':
      carregaimagem(&R, &G, &B, arq_in, &largura, &altura, &maxcor);      
      salvaimagem(6, R, G, B, txtcat(pref_out, ".ppm"), largura, altura, maxcor);
      break;
      
    case '2':
      carregaimagem(&R, &G, &B, arq_in, &largura, &altura, &maxcor);
      for (i=0; i<(largura*altura); i++) R[i] = (R[i]+G[i]+B[i])/3;
      salvaimagem(5, R, NULL, NULL, txtcat(pref_out, ".pgm"), largura, altura, maxcor);
      break;

    case '3':
      carregaimagem(&R, &G, &B, arq_in, &largura, &altura, &maxcor);
      pinta(100,90,R,G,B,largura,altura,maxcor/4);
      pinta(-100,45,R,G,B,largura,altura,maxcor/4);
      pinta(200,45,R,G,B,largura,altura,maxcor/4);
      pinta(100,60,R,G,B,largura,altura,maxcor/4);
      pinta(100,63,R,G,B,largura,altura,maxcor/4);
      salvaimagem(5, R, NULL, NULL, txtcat(pref_out, ".pgm"), largura, altura, maxcor);
      break;

    case '4':
      carregaimagem(&R, &G, &B, arq_in, &largura, &altura, &maxcor);
      m = extrailinhas(10,R,largura,altura);
      for (i=0; i<10; i++)
	{
	  pinta(m[i],m[2*i],R,NULL,NULL,largura,altura,maxcor);
	}
      free(m);
      salvaimagem(5, R, NULL, NULL, txtcat(pref_out, ".pgm"), largura, altura, maxcor);
      break;

    case '5':
      carregaimagem(&R, &G, &B, arq_in, &largura, &altura, &maxcor);
      alocaimagemcinza(&D, largura, altura);
      detectalinhas(R,largura,altura,D);
      salvaimagem(5, D, NULL, NULL, txtcat(pref_out, ".pgm"), largura, altura, maxcor);
      break;
      
    case '6':
      carregaimagem(&R, &G, &B, arq_in, &largura, &altura, &maxcor);
      alocaimagemcinza(&D, largura, altura);
      alocaimagemcinza(&E, largura, altura);
      detecta_linhas_usando_pesos(R,largura,altura,D,E,NULL);
      salvaimagem(5, D, NULL, NULL, txtcat(pref_out, ".pgm"), largura, altura, maxcor);
      salvaimagem(5, E, NULL, NULL, txtcat(pref_out, "-t.pgm"), largura, altura, maxcor);
      break;
      
    case '7':
      carregaimagem(&R, &G, &B, arq_in, &largura, &altura, &maxcor);
      detectacruzamentos(R,largura,altura,ld,cd,&mi,&mj);
      printf("Ponto aproximado: (%d, %d)\n", ld, cd);
      printf("Ponto exato: (%g, %g)\n\n", mi, mj);
      break;

    case '8':
      carregaimagem(&R, &G, &B, arq_in, &largura, &altura, &maxcor);
      mm = detectalinhas2(10,R,largura,altura);
      for (i=1; i<=10; i++)
	{
	  //	  printf("teta: (%d)  L: (%d)\n", mm[i], mm[i+10]);
	  pinta(mm[i+10],mm[i],R,NULL,NULL,largura,altura,maxcor);
	}
      free(mm);
      salvaimagem(5, R, NULL, NULL, txtcat(pref_out, ".pgm"), largura, altura, maxcor);
      break;
      
    case '9':
      carregaimagem(&R, &G, &B, arq_in, &largura, &altura, &maxcor);
      alocaimagemcinza(&D, largura, altura);
      alocaimagemcinza(&E, largura, altura);
      mm = extrai_linhas_com_pesos(N_LINHAS_SAIDA,R,D,E,largura,altura,pref_out);
      salvaimagem(5, D, NULL, NULL, txtcat(pref_out, "-local.pgm"), largura, altura, maxcor);
      salvaimagem(5, E, NULL, NULL, txtcat(pref_out, "-teta.pgm"), largura, altura, maxcor);
      nlinhas = (int)(mm[0] + 0.5);
      printf("melhores retas encontradas:\n");
      for (i=1; i<=mm[0]; i++)
	{
	  double teta = mm[i];
          double dist = mm[i+nlinhas];
          double peso = mm[i+2*nlinhas];
          printf("teta = %f  L = %f  peso = %f\n", teta, dist, peso);
	  pinta(dist,teta,D,NULL,NULL,largura,altura,maxcor);
	  pinta(dist,teta,D,NULL,NULL,largura,altura,maxcor);
	}
      if (mm) free(mm);
      salvaimagem(5, D, NULL, NULL, txtcat(pref_out, ".pgm"), largura, altura, maxcor);

      break;

    case 'H':
      printf("---------HELP--------\n");
      printf("Opção");
      printf(" 1 - Um novo arquivo é criado, com o mesmo conteúdo do arquivo original\n");
      printf(" 2 - Cria uma nova imagem PGM na qual cada pixel é a média dos três canais da imagem de entrada PPM\n");
      printf(" 3 - Pinta uma linha branca na imagem de entrada, a partir de um ângulo(teta) e uma distância (L) em relação à origem\n");
      printf(" 4 - Pinta as linhas reconhecidas automaticamente(ñ funciona).\n");
      printf(" 5 - A partir da entrada PPM, gera uma imagem PGM com resultado do detetor local, isto é,\n");
      printf("     dá uma nota para cada pixel. Esta nota é maior, quanto mais provável que este pixel \n");
      printf("     pertença a alguma linha do grid. (ñ está bom. Opção 6 tenta corrigir.\n");
      printf(" 6 - Produz imagem PGM com resultado do detetor local com pesos.\n");
      printf(" 7 - Detecta um cruzamento aproximado com precisão de meio pixel.\n");
      printf(" 8 - Pinta as linhas reconhecidas automaticamente(em testes).\n");
      printf(" 9 - Pinta as linhas reconhecidas automaticamente(usando máscaras de pesos).\n");
      printf(" H - Mostra esta mini Ajuda.\n\n\n");
      break;

    }      

  if (R) free(R);
  if (G) free(G);
  if (B) free(B);
  if (D) free(D);
  if (E) free(E);

  return 0;
}

char *txtcat (char *a, char *b)
  { char *r = malloc(strlen(a)+strlen(b)+1);
    if (r == NULL) { printf("txtcat - Não foi possivel alocar cadeia"); }
    strcpy(r, a);
    strcat(r, b);
    return(r);
  }

/*ROTINA QUE PULA L LINHAS COM O PONTEIRO FP*/
int pulalinha(FILE *fp, int l)
{
  char c;
  int i=0;

  fscanf(fp, "%c", &c);
  for (i=0; i<l; i++)
    {
      while (c!='\n') fscanf(fp, "%c", &c);
    }

  return l;
}

/*CARREGA O ARQUIVO DE IMAGEM NOMEARQ NOS VETORES RGB(ALOCADOS AQUI)*/
int carregaimagem(int **RR, int **GG, int **BB, char *nomearq, int *llargura, int *aaltura, int *mmaxcor)
{
  FILE *fp;
  int i=0, c=0;
  int largura, altura, maxcor;
  int n, tipo;
  char p;
  int *R=*RR, *G=*GG, *B=*BB;

  if ( !(fp = fopen(nomearq, "rb")) )
    {
      printf("%s : arquivo inexistente.\n", nomearq);
      return -1;
    }
  
  fscanf(fp, "%c", &p);
  if ( p!='P' )
    {
      printf("O arquivo %s não é um arquivo de imagem válido.\n", nomearq);
      return -2;
    }
  fscanf(fp, "%d\n", &tipo);/*le o tipo de arquivo: P3, P5, P6 etc.*/

  //pulalinha(fp, 1);
  fscanf(fp, "%d %d\n%d\n", &largura, &altura, &maxcor);
  if (llargura) { *llargura = largura; }
  if (aaltura) { *aaltura = altura; }
  if (mmaxcor) { *mmaxcor = maxcor; }

  /*aloca espaço para os 3 vetores RGB*/
  n=largura*altura;/*numero de pixels na imagem*/

  R=G=B=NULL;
  R = (int *)notnull(calloc(n, sizeof(int)), "no mem"); 
  if (tipo==6)
    {
      G = (int *)notnull(calloc(n, sizeof(int)),"no mem");
      B = (int *)notnull(calloc(n, sizeof(int))),"no mem");
    }

 *RR=R; *GG=G; *BB=B;
  
  i=0;  
  if ( R && G && B )
    while(c!=EOF)
      {
	if ((c=getc(fp))!=EOF) { R[i]= c; } else { break; }
	if ((c=getc(fp))!=EOF) { G[i]= c; } else { break; }
	if ((c=getc(fp))!=EOF) { B[i]= c; } else { break; }
	i++;
      }
  else
    if (R)
      {
	if (G && B)
	  {
	    while(c!=EOF)
	      {
		if ((c=getc(fp))!=EOF) { R[i]= c; } else { break; }
		if ((c=getc(fp))!=EOF) { G[i]= c; } else { break; }
		if ((c=getc(fp))!=EOF) { B[i]= c; } else { break; }
		i++;
	      }
	  }
	else
	  {
	    while(c!=EOF)
	      {
		if ((c=getc(fp))!=EOF) { R[i]= c; } else { break; }
		i++;
	      }
	  }
      }
    else
      {
	printf("O primeiro vetor deve ser não nulo em carregaimagem.\n");
	return -4;
      }
  
  fclose(fp);

  if ( i != n )
    {
      printf("Talvez o arquivo tenha acabado antes do término da leitura.\n");
      printf("Isso significa que os valores de altura e largura estão inconsitentes, ou que o arquivo está truncado.\n");
    }
  
  return i;
}


/* ALOCA UMA IMAGEM EM TONS DE CINZA */
int alocaimagemcinza(int **GG, int largura, int altura)
{
  (*GG) = (int *)notnull(calloc(largura*altura, sizeof(int)), "no mem");
  return (largura*altura);
}

/*SALVA OS VETORES RGB NO ARQUIVO NOMEARQ(CRIADO AQUI)*/
int salvaimagem(int tipo, int *R, int *G, int *B, char *nomearq, int largura, int altura, int maxcor)
{
  int i;
  int tam = largura*altura;
  FILE *fp;
  
  switch (tipo)
    {
    case 5:/*PGM*/
      if ( (!nomearq) || (!R) )
	{printf("Inconsistência de dados em salvaimagem.\n"); return -2;}
      if ( !(fp=fopen(nomearq, "wb")) )
	{printf("O arquivo %s não pode ser criado.\n", nomearq); return -1;}
      fprintf(fp, "P5\n");
      //fprintf(fp, "# Created by Amandio Sena\n");
      fprintf(fp, "%d %d\n%d\n", largura, altura, maxcor);
      for (i=0; i<tam; i++)
	{
	  fputc(R[i], fp);
	}
      fclose(fp);
      
      break;      
      
    case 6:/*PPM*/
      if ( (!nomearq) || (!R) || (!G) || (!B) )
	{printf("Inconsistência de dados em salvaimagem.\n"); return -2;}
      if ( !(fp=fopen(nomearq, "wb")) )
	{printf("O arquivo %s não pode ser criado.\n", nomearq); return -1;}
      fprintf(fp, "P6\n");
      fprintf(fp, "# Created by Amandio Sena\n");
      fprintf(fp, "%d %d\n%d\n", largura, altura, maxcor);
      for (i=0; i<tam; i++)
	{
	  fputc(R[i], fp);
	  fputc(G[i], fp);
	  fputc(B[i], fp);
	}
      fclose(fp);
      
      break;
      
    }

  return i;
}

/* O sentido horário é positivo.
   A reta horizontal tem ângulo zero.
   O sentido positivo de y é para baixo */
int pinta(double l, double teta, int *R, int *G, int *B, int largura, int altura, int maxcor)
{
  int tam=largura*altura;
  int i;
  double d, r=1;
  double cteta = cos((teta)/180.0*M_PI);
  double steta = sin((teta)/180.0*M_PI);
	  
  if ( R || G || B )
    {
      for (i=0; i<tam; i++)
	{
	  d=distancia(linha(i,largura),coluna(i,largura),cteta,steta);
	  if ((d<(l+r)) && (d>l-r))
	    {
	      if (R) { R[i]=maxcor; }
	      if (G) { G[i]=maxcor; }
	      if (B) { B[i]=maxcor; }
	    }
	}
    }
  return 1;  
}

/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
** 
**   Copyright © 2003 by the State University of Campinas (UNICAMP).
**
** Created on 2003 by Amandio Sena Jr., IC-UNICAMP (ra007998, PIBIC/CNPq).
** Modified by Jorge Stolfi, IC-UNICAMP.       
**
** Permission to use, copy, modify, and redistribute this software and
** its documentation for any purpose and without fee is hereby
** granted, provided that: (1) the copyright notice at the top of this
** file and this copyright, authorship, and warranty notice is retained
** in all derived source files and documentation; (2) no executable
** code derived from this file is published or distributed without the
** corresponding source code; and (3) these same rights are granted to
** any recipient of such code, under the same conditions.
** This software is provided "as is", WITHOUT ANY EXPLICIT OR IMPLICIT
** WARRANTIES, not even the implied warranties of merchantibility and
** fitness for a particular purpose. END OF NOTICE.
*/
