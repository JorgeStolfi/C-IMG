#ifndef imagem_H
#define imagem_H

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <listas.h>
#include <linhas.h>
#include <linhas_static.h>

char *txtcat (char *a, char *b);

int pulalinha(FILE *fp, int l);

int carregaimagem(int **RR, int **GG, int **BB, char *nomearq, int *llargura, int *aaltura, int *mmaxcor);

int salvaimagem(int tipo, int *R, int *G, int *B, char *nomearq, int largura, int altura, int maxcor);

int pinta(double l, double teta, int *R, int *G, int *B, int largura, int altura, int maxcor);

#endif
