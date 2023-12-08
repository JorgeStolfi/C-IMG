/* Last edited on 2006-04-14 09:56:09 by stolfi */

#ifndef linhas_H
#define linhas_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <listas.h>
#include <imagem.h>
#include <linhas_static.h>


double calculacurvacasamento(
  int *R, 
  int largura, 
  int v, int h, 
  double cteta, 
  double steta, 
  double *ci
);

void calcula_melhor_teta(
  int *R,              /* Imagem de entrada */ 
  int largura,         /* Largura da imagem */
  int v, int h,        /* Indices do pixel a analisar */
  double *jp,          /* Matriz de pesos dos pixels */
  double **janelas,    /* Lista de janelas */
  int NJ,              /* Número de janelas */
  double *melhor_prod, /* Melhor produto escalar dentre todas as janelas */
  int *melhor_i        /* Indice da janela que deu melhor casamento */
);

double calcula_mascara_casamento(
  int *R,             /* Imagem de entrada */ 
  int largura,        /* Largura da imagem */
  int v, int h,       /* Indices do pixel a analisar */
  double *jp,         /* Matriz de pesos dos pixels */
  double *janela      /* Mascara */
 );

/*********************************************************************/
/* Dá um erro muito esquisito quando estes dois headers estão ativos 
linhas.h:45: parse error before '*' token
linhas.h:47: parse error before '*' token
linhas.h:47: warning: type defaults to `int' in declaration of `extraipicos'
linhas.h:47: warning: data definition has no type or storage class


  E DE REPENTE TUDO FUNCIONA DE NOVO!!!!
  E DE REPENTE DEIXA DE FUNCIONAR!!!

*********************************************************************

//int *bucketing(lista_t *lista, int H, int V);

//lista_t *extraipicos(int n, int *m, int largura, int altura);

*/

void detectacruzamentos(int *R, int largura, int altura, int ld, int cd, double *pmi, double *pmj);

double *detectalinhas2(int nlinhas, int *R, int largura, int altura);

void detecta_linhas_usando_pesos(
  int *R,          /* Imagem original */
  int largura,     /*  */
  int altura,      /*  */
  int *G,          /* Imagem com valor maximo do prod. escalar (saida). */
  int *T,           /* Imagem com o teta correspondente (saida). */
  double *Tindice
);

void detectalinhas(int *R, int largura, int altura, int *G);

double *extrai_linhas_com_pesos(int n, int *R, int *G, int *T, int largura, int altura, char *pref_out);

int *extrailinhas(int n, int *R, int largura, int altura);

#endif
