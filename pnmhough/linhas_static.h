/* Last edited on 2008-07-14 20:39:48 by stolfi */

#ifndef linhas_static_H
#define linhas_static_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <listas.h>
#include <imagem.h>
#include <linhas.h>

/* Raio da janela: */
#define JRAIO 9

/* Raio do perfil da curva (JRAIO/CRAIO é a resolucao da amostragem): */
#define CRAIO 4*JRAIO

#define TAMMAX 500

/* Sao consideradas NTETAS orientacoes */
/* variando de TETAMIN graus (inclusive) a TETAMIN + 180 graus (exclusive) */

#define NTETAS 90

#define TETAMIN -90.0

/* Desvio padrão (largura nominal) do perfil da curva, em pixels da imagem: */
#define SIGMA 1

#define linha(i,largura) (i)/largura
#define coluna(i,largura) (i)%largura
#define indice(l,c,largura) ((int)((l)*largura+(c)))

/*calcula a distância do ponto (x,y) a uma reta que passa por (0,0) de inclinação teta)*/
  /*eu acho que é assim: */
  /*d = x*sin(teta)-y*cos(teta) */
  /* o stolfi disse que era -x.sen + y.cos */
  /* ele tbm disse que tem que dividir d por um fator o ... */
  /*tô fazendo do jeito do stolfi*/
/* #define distancia(y,x,cteta,steta) ((x)*(cteta)-(y)*(steta)) */
#define distancia(y,x,cteta,steta) ((x)*(steta)-(y)*(cteta))

#define cos_graus(teta_graus) cos((teta_graus)/180.0*M_PI)
#define sin_graus(teta_graus) sin((teta_graus)/180.0*M_PI)
#define cos_rad(teta_rad) cos(teta_rad)
#define sin_rad(teta_rad) sin(teta_rad)

int arredonda(double d);

void normaliza_curva(double *curva, int n);

void normaliza_mascara(double *mascara, int n, double *jp);

double *curvaideal(double d, double o, int *R, int largura);

double *janelaideal(double teta, double *ci, double *jp);

double **cria_janelas_ideais(
  int ntetas,      /* Número de angulos distintos */
  double *ci,      /* Perfil ideal da linha (gaussiana unidimensional) */
  double *jp       /* Pesos dos pixels na janela */
);

double calcula_max_prod_janelas(int njanelas, double **janelas, double *jp);

double calcula_max_prod_janela(double *janela, double *jp);

void grava_janelas_ideais(int njanelas, double **janelas, double *jp);

int libera_janelas_ideais(int njanelas, double **janelas);

double *janela_de_pesos(double o, int *R, int largura);

#endif
