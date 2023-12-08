#ifndef listas_H
#define listas_H

#include <stdio.h>
#include <malloc.h>
#include <linhas.h>

/* In�cio das defini��es de Tipos e Estruturas */

/*um n� representa um pixel(i,j) e a curva de intensidade a partir de um �ngulo teta*/
typedef struct struct_no_t
{
  double c;
  double teta;
  int v;
  int h;
  /* int largura; */
  struct struct_no_t *prox;
  struct struct_no_t *ant;
} no_t;

/*uma lista cont�m um ponteiro para o primeiro n�, o tamanho atual da lista e seu tamanho m�ximo poss�vel(definido arbitrariamente).*/
/*cmin � o pior valor de c(semelhan�a da curva representada com a curva ideal) encontrado, ou seja, o valor que est�h no �ltimo n�.*/
typedef struct struct_lista_t
{
  int tam;
  int tammax;
  double cmin;
  no_t *primno;
} lista_t;


/* Fim das defini��es de Tipos e Estruturas */


/* In�cio das defini��es de fun��es */

lista_t *novalista(int tammax);

no_t *novono();

int liberano(no_t *no);

int libno(no_t *no);

int liberalista(lista_t *lista);

int insereno(no_t *no, lista_t *lista);

/* Fim das defini��es de fun��es */

#endif
