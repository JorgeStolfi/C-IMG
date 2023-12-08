/* Last edited on 2023-02-25 15:59:37 by stolfi */
#include <linhas_static.h>
#include <affirm.h>
#include <rmxn.h>

int arredonda(double d)
{
  return (int)floor(d + 0.5); 
}

/* normaliza o vetor *curva de tamanho n */
void normaliza_curva(double *curva, int n)
{
  double s = 0.0;
  int i;

  /* Faz com que a m�dia seja 0: */
  s = 0.0;
  for (i = 0; i < n; i++) { s += curva[i]; }
  s /= (double)n;
  for (i = 0; i < n; i++) { curva[i] -= s; }
  /* Faz com que a vari�ncia seja 1: */
  s = 0.0;
  for (i = 0; i < n; i++) { s += curva[i]*curva[i]; }
  s /= (double)n;
  s = sqrt(s);
  for (i = 0; i < n; i++) { curva[i] /= s; }
}

/* normaliza o vetor *mascara de tamanho n */
void normaliza_mascara(double *mascara, int n, double *jp)
{
  double s, sp;
  int i;

  /* Torna a mascara ortogonal � func�o constante: */
  s = 0.0; sp = 0.0;
  for (i = 0; i < n; i++) { s += jp[i]*mascara[i]; sp += jp[i]; }
  s /= sp;
  for (i = 0; i < n; i++) { mascara[i] -= s; }
  /* Faz com que o produto escalar da mascara com ela mesma seja 1: */
  s = 0.0;
  for (i = 0; i < n; i++) { s += jp[i]*mascara[i]*mascara[i]; }
  s = sqrt(s);
  for (i = 0; i < n; i++) { mascara[i] /= s; }
}


/* Fun��o que devolve um vetor de double, contendo uma curva gaussiana ideal,
   normalizada.
   O valor de c[i] � a intensidade ideal da linha, medida num ponto
   cuja dist�ncia a partir do eixo da linha � (i-CRAIO)*JRAIO/CRAIO
   Para cada elemento c[i], i representa o eixo X e c[i] seria o valor de Y para este X */
/* R cont�m a imagem;
   largura = largura da imagem em pixels;
   o = fator de corre��o da curva (sigma). Ajustado como 1 */
double *curvaideal(double d, double sigma, int *R, int largura)
{
  double *ci;
  int i;
  int n = 2*CRAIO + 1; 
  double passo = ((double)JRAIO)/((double)CRAIO);

  if (sigma==0)
    {
      printf("Erro de divis�o por zero em curvaideal.\n");
      exit (1);
    }

  ci = rn_alloc(n);
  
  for (i=-CRAIO; i<=CRAIO; i++)
    {
      double y = ((double)(passo*i-d))/sigma;
      ci[i+CRAIO] = exp(-y*y);
    }
  normaliza_curva(ci,n);
  return ci;
}


/* Esta fun��o deve, a partir da curva ideal, gerar uma m�scara, isto �, uma janela ideal para determinado teta */
/* Deve rodar somente uma vez pra cada teta. Isto j� � feito por causa da chamada em cria_janelas_ideais. */
/* Esta fun��o pode ser resumida assim:
   pega uma curva ideal e a projeta em toda uma janela */
/*   Deve ser do mesmo tamanho que a janeladepesos: nXn (n=2*JRAIO+1) */
/* ci = curva ideal (gerado por outra fun��o) */
double *janelaideal(double teta, double *ci, double *jp)
{
  int n=2*JRAIO+1; /* largura e altura da janela */
  int tam = (2*JRAIO+1)*(2*JRAIO+1);
  int i=0; /* contador */
  double cteta = cos(teta/180.0*M_PI);
  double steta = sin(teta/180.0*M_PI);
  int indice=0;/* s� para tratar erro */
  double passo = ((double)JRAIO)/((double)CRAIO);

  /*  gera uma janela nova toda em branco;
      pra cada pixel da janela:
         verifica a dist�ncia desse pixel � reta que passa pelo centro com inclina��o teta;
	 arredonda esta dist�ncia para o inteiro mais pr�ximo;
	 esse valor � o �ndice da curva a ser buscado(respeitando os limites);
	 o valor em curva[indice] dever� ser pintado no pixel em quest�o.  */

  double *janela = rn_alloc(tam);
  for (i=0; i<tam; i++) { janela[i]=0; }

  for (i=0; i<tam; i++)
    {
      int l = i/n - JRAIO;
      int c = i%n - JRAIO;
      if (l*l + c*c <= JRAIO*JRAIO)
        {
	  double dist = distancia((double)l,(double)c,cteta,steta);
          indice = arredonda(dist/passo) + CRAIO;
          if (indice>(2*CRAIO+1)) { printf("Erro em janela ideal.(%d)\n", indice); }
          janela[i] = ci[indice];
        }
      else
        { janela[i] = 0.0; }
    }

  normaliza_mascara(janela,tam, jp);
  return janela;
}

/* Acrescentar a consistencia de intervalos, isto �, caso o intervalo desejado n�o seja de um em um, deve-se fazer algumas altera��es.
   Na verdade, se o intervalo for qualquer inteiro, n�o s�o necess�rias mudan�as,
   simplesmente haver� um pequeno desperd�cio de mem�ria, pois alguns elementos do vetor de janelas ficar�o vazios.
   Muda a implementa��o se os intervalos puderem ser com valores decimais de teta.
   Uma sugest�o seria criar um tipo janela, no qual o teta seria armazenado,
   junto com outros par�metros, como os pr�pios pixels e coordenadas reais do ponto central,
   algo que seria inutilizado nesta fun��o de cria janelas, mas seria �til ao percorrer a imagem. */

/* devolve um vetor, no qual cada elemento guarda uma janelaideal (um vetor de double) */
/* considera-se que os intervalos entre �ngulos s�o inteiros */
/* ci cont�m uma curva ideal */
double **cria_janelas_ideais(
  int ntetas,      /* N�mero de angulos distintos */
  double *ci,      /* Perfil ideal da linha (gaussiana unidimensional) */
  double *jp       /* Pesos dos pixels na janela */
)
{
  double **janelas;
  double passo_teta = 180.0/((double)ntetas);
  
  int i=0;/*contador*/

  /* aloca uma janela para cada teta */
  if ( (janelas = (double **)calloc(ntetas, sizeof(double *))) == NULL )
    printf("Problemas ao alocar o vetor de janelas.\n");
    
  for (i=0; i<ntetas; i++) { janelas[i]=NULL; }
  /* confirma que todas t�m NULL */
  /* esta garantia pode ser de alguma forma �til fora desta fun��o */

  /* para cada teta pendura uma janelaideal em janelas */
  /* janelas[teta] aponta para uma janelaideal */
  /* contando tetas de 1 em 1 */
  /* ficou fora do for de inicializa��o para permitir altera��es no intervalo de tetas */
  /* como esta fun��o � chamada somente uma vez, o desempenho n�o � mais importante que a capacidade de manuten��o e extens�o */
  for (i=0; i<ntetas; i++)
    {
      janelas[i] = janelaideal(TETAMIN + i*passo_teta, ci, jp);
    }

  return janelas;
}

double calcula_max_prod_janelas(int njanelas, double **janelas, double *jp)
{
  double max_prod = 0.0;
  int i;
  for (i = 0; i < njanelas; i++)
    { double mp = calcula_max_prod_janela(janelas[i], jp);
      if (mp > max_prod) { max_prod = mp; }
    }
  return max_prod;
}

double calcula_max_prod_janela(double *janela, double *jp)
{
  double maxp = 0.0, minp = 0.0;
  int h, v;
  int n = 2*JRAIO + 1;
  for (v = 0; v < n; v++)
    for (h = 0; h < n; h++)
      { int ind = indice(v,h,n);
        double jan = janela[ind];
        double pes = jp[ind];
        double t = jan*pes;
        if (t > 0)
          { maxp += t*255; }
        else
          { minp += t*255; }
      }
  return (maxp > -minp ? maxp : -minp);
}

void grava_janelas_ideais(int njanelas, double **janelas, double *jp)
{
  /* Aloca imagem com todas as m�scaras: */

  int n = 2*JRAIO + 1;
  int M_janelas_h = (int)ceil(sqrt(njanelas));
  int M_janelas_v = (int)ceil(((double)njanelas)/M_janelas_h);
  int M_altura = (n+1)*M_janelas_v;
  int M_largura = (n+1)*M_janelas_h;
  int *M = (int *)malloc(M_altura*M_largura*sizeof(int));
  int i, h, v;

  if (M == NULL) printf("Problemas ao alocar imagem de janelas.\n");
  for (i = 0; i < njanelas; i++)
    {
      int Mi = i / M_janelas_h;
      int Mj = i % M_janelas_h;
      /* Calcula valor maximo da janela para ajuste de tom: */
      double tmax = 0.0;
      for (v = 0; v < n; v++)
        for (h = 0; h < n; h++)
          {
	    double t = fabs(janelas[i][indice(v,h,n)]);
            if (t > tmax) { tmax = t; }
          }
      
      /* Converte janela[i] para tons de cinza e guarda em M: */
      for (v = 0; v < n; v++)
        for (h = 0; h < n; h++)
          { 
            int ind = indice(v,h,n);
            double t = janelas[i][ind]*jp[ind];
            double cinza = 0.5*(t/tmax) + 0.5;
            if (cinza < 0.0) { cinza = 0.0; }
            if (cinza > 1.0) { cinza = 1.0; }
            M[indice(Mi*(n+1)+v, Mj*(n+1)+h, M_largura)] = (int)(255.0*cinza + 0.5);
          }
    }
  salvaimagem(5, M, NULL, NULL, "mascaras.pgm", M_largura, M_altura, 255);
  free(M);
}


int libera_janelas_ideais(int njanelas, double **janelas)
{
  int i=0;

  if (janelas)
    {
      for (i=0; i<njanelas; i++)
	{
	  if (janelas[i]) { free(janelas[i]); }
	}
    }
/*   else */
/*     { */
/*       printf("Passagem de argumentos incorreta em libera_janelas_ideais.\n"); */
/*       return -1; */
/*     } */

  return 1;
}


/* atribui um peso para cada pixel da janela */
/* janela criada � nXn (n=2*JRAIO+1) */
double *janela_de_pesos(double sigma, int *R, int largura)
{
  
  int i,j;
  int n = 2*JRAIO+1;
  double sigma2 = sigma*sigma;

  if (sigma==0)
    {
      printf("Erro de divis�o por zero em curvaideal.\n");
      exit (1);
    }

  double *jp = rmxn_alloc(n,n);
  
  /* pesos atribu�dos conforme uma distribui��o gaussiana */
  /* fica como um chap�u de mexicano */
  for (i=-JRAIO; i<=JRAIO; i++)
    for (j=-JRAIO; j<=JRAIO; j++)
      { int ind = n*(i+JRAIO) + (j+JRAIO);
        if (i*i + j*j <= JRAIO*JRAIO)
          {
            double x = (double)j;
            double y = (double)i;
            double r2 = x*x + y*y;
            jp[ind] = exp(-r2/sigma2);
          }
        else
          { jp[ind] = 0.0; }
      }
  return jp;
}
