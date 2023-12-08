#define PROG_NAME "mavica-calib"
#define PROG_DESC "process Mavica CD-300 images of the Kodak calibration chart"
#define PROG_VERS "1.1"

/* Copyright © 2003 by the Fluminense Federal University (UFF), Niterói.
** See the copyright, authorship, and warranty notice at end of file.
** Last edited on 2017-07-28 02:36:49 by jstolfi
*/

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  { -r | --run } RATEFILE" 

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define NUM_DEGRAUS 20

typedef struct rgb  /* Cor RGB. */
  { 
    double r, g, b;
  } RGB;

typedef struct degrau /* Dados referentes a um degrau da escala Kodak. */
  { 
    /* Valores observados (extraídos da imagem): */
    RGB valor;            /* Cor do degrau, [0 _ 255]. */
    RGB papel;            /* Cor do fundo branco próximo ao degrau, [0 _ 255]. */
    /* Valores nominais (segundo documentação da escala Kodak): */
    double densidade;     /* Densidade fotográfica do degrau, [0.0 _ 2.0]. */
    RGB refletancia;      /* Refletância equivalente, [0.0 _ 1.0]. */
    RGB cor_corrigida;    /* Cor esperada na imagem corrigida, [0 _ 255]. */
  } Degrau;

/* PROTOTYPES */

int main (int argc, char** argv);

void avalia_args (char** arg, char **nome_P);
  /* Analisa as opções e argumentos fornecidos na linha de comando.
    Caso a opção usada seja "-r", devolve em {*nome_P} o nome
    da imagem (sem a extensão). */

void le_arq_med (FILE* arq, Degrau* esc);
  /* Lê as cores observadas na imagem da escala Kodak. 
    Para cada linha, lê o número do degrau (em 0..19), sua densidade
    fotográfica nominal (no intervalo [0.0 _ 2.0]), sua cor R,G,B
    observada (em [0.0 _ 255.0]^3), e a cor R,G,B do fundo branco
    adjacente ao degrau (idem). */

void escreve_arq_res (FILE* arq, Degrau* esc);
  /* Grava arquivo com valores RGB ideais de cada degrau. */

RGB calcula_refletancia (double densidade);
  /* Calcula a refletância nominal de um degrau da escala de cinza, a
    partir da sua densidade fotográfica. */

RGB calcula_cor_corrigida (RGB refletancia);
  /* Calcula a cor RGB ideal, na imagem corrigida, para um degrau da
    escala de cinza com a refletância dada. */
    
RGB calcula_refletancia_aparente(RGB valor, RGB papel, RGB preto);
  /* Devolve a razão {(valor-preto)/(papel-preto)}. O parâmetro
    {preto} deve ser o valor RGB registrado pela câmera ao fotografar
    uma cena completamente negra. */

RGB estima_gama(Degrau *e0, Degrau *e1, RGB preto);
  /* Estima o expoente {gama} a partir da variação de cor aparente
    observada entre os degraus dados {e0} e {e1}. O parâmetro {preto}
    deve ser o valor RGB registrado pela câmera numa cena totalmente
    negra. */

RGB log10_rgb (RGB x);
  /* Devolve o logaritmo decimal de cada componente. */

double refl_aparente(double valor, double papel, double preto);
  /* Devolve a razão {(valor-preto)/(papel-preto)}, para um único
    canal de cor (R, G, ou B). O parâmetro {preto} deve ser o valor
    registrado pela câmera ao fotografar uma cena completamente
    negra. */

void escreve_arq_plo (FILE* arq, Degrau* esc);
  /* Grava arquivo para gnuplot, com log(refletancia nominal) e
    log(refletância medida), para determinação do gama
    da câmera. */

FILE* abre_arquivo (char *dir, char* nome, char *ext, char *modo);
  /* Abre o arquivo "{dir}/{nome}.{ext}", com modo {modo}. */

/* IMPLEMENTATIONS */

int main (int argc, char** argv)
{
  Degrau escala[NUM_DEGRAUS];
  char *nome;

  avalia_args(argv, &nome);

  FILE *arq_med = abre_arquivo("in", nome, "-med.txt", "rt");
  FILE *arq_plo = abre_arquivo("out", nome, "-plo.txt", "wt");
  FILE *arq_res = abre_arquivo("out", nome, "-res.txt", "wt");

  le_arq_med(arq_med, escala);
  fclose(arq_med);

  escreve_arq_plo(arq_plo, escala);
  fclose(arq_plo);

  escreve_arq_res(arq_res, escala);
  fclose(arq_res); 

  return 0;
}

void le_arq_med (FILE* arq, Degrau* esc)
{
  int i, lin;

  for (i=0; i<NUM_DEGRAUS; i++)
    { 
      RGB *vi = &(esc[i].valor);
      RGB *pi = &(esc[i].papel);
      fscanf(arq, 
        "%d %lf %lf %lf %lf %lf %lf %lf", 
        &lin,                   /* Número do degrau, 00..19. */     
        &esc[i].densidade,      /* Densidade nominal do degrau. */
        &vi->r, &vi->g, &vi->b, /* RGB do degrau na imagem. */
        &pi->r, &pi->g, &pi->b  /* RGB do fundo na imagem. */
      );
      esc[i].refletancia = calcula_refletancia(esc[i].densidade);
      esc[i].cor_corrigida = calcula_cor_corrigida(esc[i].refletancia);
      assert(lin == i);
    }
}

void escreve_arq_plo (FILE* arq, Degrau* esc)
{
  /* Estima o RGB nulo {preto} da câmera olhando para o degrau [19]: */
  RGB preto;
  preto.r = 0.8*esc[19].valor.r;
  preto.g = 0.8*esc[19].valor.g;
  preto.b = 0.8*esc[19].valor.b;
  fprintf(arq, "# Preto (im.orig.) = %6.2f %6.2f %6.2f\n", preto.r, preto.g, preto.b);

  /* Plota refletância aparente contra refletância nominal: */
  fprintf(arq, "# \n");
  fprintf(arq, "# Cada linha se refere a um degrau da escala Kodak:\n");
  fprintf(arq, "# log(refl_nominal.{r,g,b}) | log(refl_aparente.{r,g,b})\n");
  int i;
  for (i=0; i<NUM_DEGRAUS; i++)
    {
      RGB vi = esc[i].valor;        /* Cor obs. do degrau {i}. */
      RGB pi = esc[i].papel;        /* Cor obs. do fundo branco adjacente. */
      RGB rni = esc[i].refletancia; /* Refletancia nominal do degrau. */
      RGB rai = calcula_refletancia_aparente(vi, pi, preto);
      RGB x = log10_rgb(rni); /* log10(refl. nominal). */
      RGB y = log10_rgb(rai); /* log10(refl. aparente). */
      fprintf(arq, "%7.4f %7.4f %7.4f  %7.4f %7.4f %7.4f\n", 
        x.r, x.g, x.b, 
        y.r, y.g, y.b);
    }
}

void escreve_arq_res (FILE* arq, Degrau* esc)
{
  /* Estima o RGB nulo {preto} da câmera olhando para o degrau [19]: */
  RGB preto;
  preto.r = 0.8*esc[19].valor.r;
  preto.g = 0.8*esc[19].valor.g;
  preto.b = 0.8*esc[19].valor.b;
  fprintf(arq, "# Preto (im.orig.) = %6.2f %6.2f %6.2f\n", preto.r, preto.g, preto.b);

  /* Estima o expoente {gama} olhando para os degraus [0] e [5]: */
  RGB gama = estima_gama(&(esc[0]), &(esc[5]), preto);
  fprintf(arq, "# Gama (câmera) = %6.2f %6.2f %6.2f \n", gama.r, gama.g, gama.b);

  /* Calcula a cor esperada {pap_corr} do papel na imagem corrigida: */
  RGB alfa = (RGB){0.95, 0.95, 0.96}; /* Refletância estimada do papel. */
  RGB pap_corr = calcula_cor_corrigida(alfa);
  fprintf(arq, "# Papel (im.corr.) = %f %f %f\n", pap_corr.r, pap_corr.g, pap_corr.b);

  /* Grava cores corrigidas no cabeçalho: */
  fprintf(arq, "# \n");
  fprintf(arq, "# Cada linha se refere a um degrau da escala Kodak:\n");
  fprintf(arq, "# densidade | refletancia.{r,g,b} | cor_corrigida.{r,g,b}\n");
  int i;
  for (i=0; i<NUM_DEGRAUS; i++)
    fprintf(arq, "%4.2f  %7.4f %7.4f %7.4f  %6.2f %6.2f %6.2f\n", 
      esc[i].densidade, 
      esc[i].refletancia.r,   esc[i].refletancia.g, esc[i].refletancia.b,  
      esc[i].cor_corrigida.r, esc[i].cor_corrigida.g, esc[i].cor_corrigida.b
    );
}

void avalia_args (char** arg, char **nome_P)
{
  if ((strcmp(arg[1], "--help") == 0) || (strcmp(arg[1], "-h") == 0))
    {
      fprintf(stderr, "  %s - %s\n", PROG_NAME, PROG_DESC);
      fprintf(stderr, "  Version %s\n\n", PROG_VERS);
      fprintf(stderr, "  Copyright © 2003 Diego Araujo, Gilberto Medeiros\n\n");
      fprintf(stderr, "  Usage: %s\n", PROG_HELP);
      exit(0);
    }
  else if ((strcmp(arg[1], "--run") == 0) || (strcmp(arg[1], "-r") == 0))
    { 
      (*nome_P) = arg[2]; 
      return;
    }
  else
    { 
      fprintf(stderr, "Invalid command line.\n\n");
      exit(1);
    }
}

FILE* abre_arquivo (char *dir, char* nome, char *ext, char *modo)
{
  char *nome_arq = NULL;
  asprintf(&nome_arq, "%s/%s%s", dir, nome, ext);
  FILE *arq = fopen(nome_arq, modo);
  if (arq == NULL)
    {
      fprintf(stderr, "File %s could not be opened.\n", nome_arq);
      exit(0);
    }
  free(nome_arq);
  return arq;
}

RGB calcula_refletancia (double densidade)
{
  double y = pow(10,-densidade);
  return (RGB){y, y, y};
}

RGB calcula_refletancia_aparente (RGB cor_med, RGB papel_med, RGB preto)
{
  RGB ra;
  ra.r = refl_aparente(cor_med.r, papel_med.r, preto.r);
  ra.g = refl_aparente(cor_med.g, papel_med.g, preto.g); 
  ra.b = refl_aparente(cor_med.b, papel_med.b, preto.b);
  return ra;
}

double refl_aparente(double valor, double papel, double preto)
{
  valor = valor - preto;
  if (valor < 1.0) { valor = 1.0; }
  papel = papel - preto;
  if (papel < 1.0) { papel = 1.0; }
  return valor/papel;
}

RGB calcula_cor_corrigida (RGB refletancia)
{
  RGB cor;
  cor.r = refletancia.r*255;
  cor.g = refletancia.g*255;
  cor.b = refletancia.b*255;
  return cor;
}

RGB log10_rgb (RGB x)
{
  RGB y;
  y.r = log10(x.r);
  y.g = log10(x.g);
  y.b = log10(x.b);
  return y;
}

RGB estima_gama(Degrau *e0, Degrau *e1, RGB preto)
{
  RGB v0 = e0->valor;        /* Cor do degrau {e0}. */
  RGB p0 = e0->papel;        /* Cor do fundo branco adjacente. */
  RGB rn0 = e0->refletancia; /* Refletancia nominal. */
  RGB ra0 = calcula_refletancia_aparente(v0, p0, preto); /* Refl. aparente. */
  RGB x0 = log10_rgb(rn0);
  RGB y0 = log10_rgb(ra0);
      
  RGB v1 = e1->valor;        /* Cor do degrau {e1}. */
  RGB p1 = e1->papel;        /* Cor do fundo branco adjacente. */
  RGB rn1 = e1->refletancia; /* Refletancia nominal. */
  RGB ra1 = calcula_refletancia_aparente(v1, p1, preto); /* Refl. aparente. */
  RGB x1 = log10_rgb(rn1);
  RGB y1 = log10_rgb(ra1);
      
  RGB gama;
  gama.r = (y1.r - y0.r)/(x1.r - x0.r);
  gama.g = (y1.g - y0.g)/(x1.g - x0.g);
  gama.b = (y1.b - y0.b)/(x1.b - x0.b);
  
  return gama;
}

/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
** 
**   Copyright © 2003 by the Fluminense Federal University (UFF), Niterói.
**
** Created 2003 by Diego Araújo, Gilberto Medeiros, UFF.
** Modified 2004-2005 by Jorge Stolfi, IC-UNICAMP.       
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
