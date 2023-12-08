#ifndef listas_H
#define listas_H

#include <stdio.h>
#include <malloc.h>
#include <linhas.h>

/* Início das definições de Tipos e Estruturas */

/*um nó representa um pixel(i,j) e a curva de intensidade a partir de um ângulo teta*/
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

/*uma lista contém um ponteiro para o primeiro nó, o tamanho atual da lista e seu tamanho máximo possível(definido arbitrariamente).*/
/*cmin é o pior valor de c(semelhança da curva representada com a curva ideal) encontrado, ou seja, o valor que estáh no último nó.*/
typedef struct struct_lista_t
{
  int tam;
  int tammax;
  double cmin;
  no_t *primno;
} lista_t;


/* Fim das definições de Tipos e Estruturas */


/* Início das definições de funções */

lista_t *novalista(int tammax);

no_t *novono();

int liberano(no_t *no);

int libno(no_t *no);

int liberalista(lista_t *lista);

int insereno(no_t *no, lista_t *lista);

/* Fim das definições de funções */

#endif
