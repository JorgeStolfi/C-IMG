/* Last edited on 2023-02-25 16:05:10 by stolfi */

#include <linhas.h>
#include <affirm.h>
#include <rn.h>

static double *curva = NULL;

/* Esta fun��o funciona. Meio precariamente, mas funciona. */
double calculacurvacasamento(
  int *R,
  int largura, 
  int v, int h, 
  double cteta, 
  double steta, 
  double *ci
)
{
  double d=0;
  int i, j, a, b;
  double c=0, Rmin, Rmax;
  int n=2*CRAIO+1;
  
  /* Calcula a curva observada `curva' (proje��o da janela na dire��o teta) */
  if (curva == NULL)
    {
      curva = (double *)calloc(n, sizeof(double));
      if (curva == NULL) { printf("Problemas ao alocar nova curva em calculacurvacasamento.\n"); }
    }
  for (i = 0; i < n; i++) { curva[i] = 0.0; }

  /* Encontra o maximo e minimo na janela: */
  Rmax = 0.0; Rmin = 255.0;
  for (i=v-JRAIO; i<=v+JRAIO; i++)
    for (j=h-JRAIO; j<=h+JRAIO; j++)
      {
	int Rij = R[indice(i,j,largura)];
        if (Rij < Rmin) { Rmin = Rij; }
        if (Rij > Rmax) { Rmax = Rij; }
      }
  Rmax += 3;
  Rmin -= 3;
      
  /*estes dois 'for' percorrem a janela quadrada*/
  for (i=v-JRAIO; i<=v+JRAIO; i++)
    for (j=h-JRAIO; j<=h+JRAIO; j++)
      {
        /*seja r a reta que passa pelo centro com inclina��o teta*/
        /*calcula a dist�ncia do ponto (i,j) � reta r */
        a=i-v;
        b=j-h;
        d = distancia(a,b,cteta,steta)/SIGMA;

        /*verifica se est� dentro do c�rculo*/
        if (a*a+b*b <= JRAIO*JRAIO)
          {
            /* Calcula o �ndice do ponto projetado no vetor `curva' */
            int di = arredonda((d*CRAIO)/JRAIO);

            double rij = (Rmax - R[indice(i,j,largura)])/(Rmax - Rmin);

            /*soma a cor desse pixel na curva, considerando a coordenada d*/
            /*qual a representa��o mais adequada da curva? */
            /*duas sugest�es:
              - vetor de double(subdividindo um inteiro em um numero finito de casas decimais(umas dez);
              - curva real, guardando a informa��o de intensidade de cada coordenada, numa lista ou vetor guardando em cada elemento a coordenada(dist�ncia exata) e a intensidade de cor .*/
            
            /*usando o vetor de double*/
            
            /*soma o valor correspondente na curva*/
            curva[di+CRAIO]+=rij;
            /*isso aqui tah cagado*/
            /*os valores positivos tao errados, eu acho*/
            /*supondo que esteja certo...*/
          }
      }
      
  /* Normaliza a curva observada, para comparar com a ideal: */
  normaliza_curva(curva,n);

  /* na verdarde, c � uma medida de discrep�ncia, isto �, qto maior, pior � o casamento com a curva ideal. */
  /* para contornar este problema, a fun��o retorna o valor negativo de c */
  c=0;
  for (i=0; i<n; i++)
    {
      double e;
      e = curva[i]-ci[i];
      c+= e*e;
      /* Deve-se fazer uma c�pia da janela com o peso de cada pixel. */
      /* Por exemplo: um pixel na borda da janela deve valer pouco na compara��o */

    }
  return -c;
}

/* Esta Fun��o � a que compara duas janelas propriamente */
/* � uma adapta��o da fun��o calculacurvacasamento */
/* jp � a janela de pesos;
   janelas[teta] cont�m a janela ideal para o dado teta (o vetor janelas � criado pela fun��o cria_janelas_ideais */
void calcula_melhor_teta(
  int *R,              /* Imagem de entrada */ 
  int largura,         /* Largura da imagem */
  int v, int h,        /* Indices do pixel a analisar */
  double *jp,          /* Matriz de pesos dos pixels */
  double **janelas,    /* Lista de janelas */
  int NJ,              /* N�mero de janelas */
  double *melhor_prod, /* Melhor produto escalar dentre todas as janelas */
  int *melhor_i        /* Indice da janela que deu melhor casamento. Este valor indica o �ngulo ideal de casamento */
)
{
  int i;
  *melhor_i = -1;
  *melhor_prod = 0.0;
  for (i = 0; i < NJ; i++)
    { double prod = calcula_mascara_casamento(R, largura, v, h, jp, janelas[i]);
      if (fabs(prod) > fabs(*melhor_prod))
        { *melhor_prod = prod; *melhor_i = i; }
    }
}

double calcula_mascara_casamento(
  int *R,             /* Imagem de entrada */ 
  int largura,        /* Largura da imagem */
  int v, int h,       /* Indices do pixel a analisar */
  double *jp,         /* Matriz de pesos dos pixels */
  double *janela      /* Mascara */
)
{
  int i, j, a, b;
  double prod_escalar = 0.0;
  int n=2*JRAIO+1;

  /*estes dois 'for' percorrem a janela quadrada*/
  for (i=v-JRAIO; i<=v+JRAIO; i++)
    for (j=h-JRAIO; j<=h+JRAIO; j++)
      {
        a=i-v;
        b=j-h;
	if (a*a + b*b <= JRAIO*JRAIO) /*dentro do c�rculo*/
          {
	    /* Acumula produto escalar da amostra com a janela ideal */
	    int p = n*(a+JRAIO)+(b+JRAIO);

            double rij = R[indice(i,j,largura)];
            double tij = janela[p];
            double wij = jp[p];
            prod_escalar += wij*rij*tij;
          }
      }
  /* prod_escalar � positivo quando a amostra � igual � janela, negativo
    quando � o oposto (negativo) da janela, zero quando n�o tem nada a ver. */
  return prod_escalar;
}


/*H e V sao, respectivamente, largura e altura(em pixels) da imagem*/
int *bucketing(lista_t *lista, int H, int V)
{
  int largura, altura, tetai=-90, tetaf=89;
  int i, j, l,c;
  no_t *p;

  largura=(tetaf-tetai)/5 +1;
  altura=arredonda(sqrt((H*H)+(V*V))/4) +1;

  int tam = largura*altura+1;
  int *m = (int *)notnull(calloc(tam=(largura*altura+1), sizeof(int)), "no mem");
  for(i=0;i<altura;i++)
    for(j=0;j<largura;j++)
      m[indice(i,j,largura)]=0;
  m[0]=--tam;

  p = lista->primno;
  for (i=0; i<lista->tam; i++)
    {
      double cteta = cos(p->teta/180.0*M_PI);
      double steta = sin(p->teta/180.0*M_PI);
      /*******************************************************************/
      /* o valor de "l" pode ser negativo. verificar incoer�ncias.
	 ser� que n�o est� tudo errado???
	 acho que n�o */
      /*******************************************************************/
      l= abs(arredonda(distancia(p->v,p->h,cteta,steta)/4));
      c= (arredonda(p->teta/5) - tetai/5);
      m[indice(l,c,largura)+1]++;
      p=p->prox;
    }
  /*os picos dessa tabela m correspondem �s retas da imagem.*/

  return m;
}


lista_t *extraipicos(int n, int *m, int largura, int altura)
{
  lista_t *lista;
  no_t *no;
  int i, tam, tetai=-90;

  lista=novalista(n);
  tam=m[0];
  for (i=0; i<tam; i++)
    {
      no = novono();
      no->c=m[i+1];
      no->teta=tetai+5*coluna(i,largura);
      no->v=4*linha(i,largura);
      insereno(no, lista);
    }

  return lista;
}


/* a imagem est� em R;
   ld e cd s�o, respectivamente, a linha e a coluna do ponto aproximado 
   pmi e pmj guardam os valores de i e j exatos */
void detectacruzamentos(int *R, int largura, int altura, int ld, int cd, double *pmi, double *pmj)
{
  int i=0, j=0, mi, mj;/*contadores*/
  double *ci;
  double pi180=180*M_PI;
  double cmaxv = -1e30, cmaxh = -1e30, c=0, k, kv=0, kh=0;
  double cteta = cos(-90/pi180);
  double steta = sin(-90/pi180);


  ci = curvaideal(0, SIGMA, R, largura);

  cmaxv = cmaxh = -1e30;
  /*Para cada pixel da janela ao redor do ponto dado*/
  for (i=ld-JRAIO; i<ld+JRAIO; i++)
    {
      fprintf(stderr, ".");
      for (j=cd-JRAIO; j<cd+JRAIO; j++)
        {
	  c=calculacurvacasamento(R,largura,i,j,cteta,steta,ci);
	  if (c > cmaxv) { cmaxv = c; mj = j; }/* melhor casamento vertical */
	  c=calculacurvacasamento(R,largura,i,j,steta,-cteta,ci);
	  if (c > cmaxh) { cmaxh = c; mi = i; }/* melhor casamento horizontal */
	}
    }
  /*aqui encontrou o melhor pixel (mi,mj) com precis�o de 1 pixel */

  cmaxv = cmaxh = -1e30;
  /* agora deve variar a curva ideal 0.5 pixel a menos e a mais. */
  for (k=-0.8; k<=0.8; k+=0.1)
    {
      if (ci) free(ci);
      ci = curvaideal(k,SIGMA,R,largura);/* curva ideal deslocada de k */
      c=calculacurvacasamento(R,largura,mi,mj,cteta,steta,ci);
      if (c > cmaxv) { cmaxv = c; kv = k; }/* melhor casamento vertical */
      c=calculacurvacasamento(R,largura,mi,mj,steta,-cteta,ci);
      if (c > cmaxh) { cmaxh = c; kh = k; }/* melhor casamento horizontal */
    }

  printf("kh: %g.  kv: %g.\n", kh, kv);

  (*pmi) = (double)mi+kh;
  (*pmj) = (double)mj+kv;
  
  if (ci) free(ci);
}

/* nlinhas � o n�mero de linhas que ele vai tentar identificar, isto �, pega os nlinhas tro�os mais parecidos com uma linha*/
/* *R, largura e altura s�o os par�metros da imagem */
double *detectalinhas2(int nlinhas, int *R, int largura, int altura)
{
  int i=0, j=0;/*contadores*/
  double teta=0, maior=0; /*maior=maior nota atribuida a um pixel. utilizado para normalizar a matriz e criar uma imagem */
  double *ci, *mm;
  /* Par�metro {L} varia em {0..Lmax-1}: */
  int Lmax = (int)ceil(sqrt(largura*largura + altura*altura));
  /* Par�metro {teta} varia em {0..Tmax-1}: */
  int Tmax = 180;
  /* Numero de baldes na transformada: */
  int n=Lmax*Tmax; /* tamanho da matriz que acumula as notas de cada pixel, jah na forma teta X L */
  int *G;/* guarda uma foto da matriz com a acumula��o das notas. Eixo vertical = teta ; horizontal = L */
  lista_t *lista;
  no_t *no;
  double *M;/* guarda a matriz com a acumula��o das notas. Eixo vertical = teta ; horizontal = L */
  int *MLINHAS;/* guarda uma copia da imagem, mas com a cor de cada pixel relacionada a seu grau de aparencia com uma linha */

  ci = curvaideal(0, SIGMA, R, largura);
  
  M = rn_alloc(n);
  G = (int *)notnull(calloc(n,sizeof(int)), "no mem");
  MLINHAS = (int *)notnull(calloc(largura*altura,sizeof(int)), "no mem");
  for (i=0; i<n; i++) { M[i]=0; }

  /*Para cada pixel da imagem grande*/
  for (i=JRAIO; i<(altura-JRAIO); i++)
    { fprintf(stderr, ".");
      for (j=JRAIO; j<(largura-JRAIO); j++)
        { 
	  double cmax = -1e30, Gv;
	  double mteta, mcteta, msteta;
          int iteta, iL;

	  /* este la�o procura aproximadamente o melhor �ngulo. Esta busca ser� refinada */
          for (teta=-90; teta<90; teta+=5)
            { 
	      double cteta = cos(teta/180.0*M_PI);
              double steta = sin(teta/180.0*M_PI);
              double c;
              c=calculacurvacasamento(R,largura,i,j,cteta,steta,ci);
              if (c > cmax)
		{
		  cmax = c;
		  mteta=teta;
		  mcteta=cteta;
		  msteta=steta;
		}
            }
          /* cmax � a melhor nota para todos os angulos teta neste pixel */
          /* cmax = 0 se o casamento foi perfeito. */
          /* cmax << 0 se o casamento foi ruim. */
          /* Gv = 1.0 se casamento perfeito, 0.0 se p�ssimo. */
          Gv = (15 - fabs(cmax))/15.0 - 0.393;
          if (Gv < 0.0) { Gv = 0.0; }
          if (Gv > 1.0) { Gv = 1.0; }
	  MLINHAS[largura*i + j] = (int)(Gv*255);
	  /*deve pegar uma matriz (o bucket) e somar as notas de cada pixel*/
	  /* essas distancias devem ser do ponto (i,j) � reta de �ngulo teta que passa por (0,0) */
          if (i == 100)
            { fprintf(stderr, "pixel = [%3d,%3d] mteta = %3.0f  Gv = %5.3f", i,j,mteta,Gv); }

	  /* Este la�o refina a busca do teta */
          for (teta=mteta-10; teta<=mteta+15; teta++)
            { 
	      double cteta = cos(teta/180.0*M_PI);
              double steta = sin(teta/180.0*M_PI);
    
              iteta = (int)(teta*(Tmax/180) + 0.5 + Tmax) % Tmax;
	      /* um valor de teta indica um �ndice da matriz, portanto, aumenta todos os valores de teta pra que a contagem se inicie em zero. */
              iL = (int)((Lmax + distancia(i,j,cteta,steta))/2 + 0.5);
              if (i == 100)
                { fprintf(stderr, " iteta = %d iL = %d\n", iteta,iL); }
              if ((iL < 0) || (iL > Lmax))
                { fprintf(stderr, "** bad iL [%d,%d] = %d **\n", i,j,iL); }
              else
                { M[Lmax*iteta + iL] += Gv; }
            }
	}
    }

  salvaimagem(5, MLINHAS, NULL, NULL, "det.pgm", largura, altura, 255);
  free(MLINHAS);

  /*o m�ximo valor de maior � de 0.972401*/
  /* Escreve matriz: */
  for (i = 0; i < n; i++) 
    { 
      double Mi = M[i];
      if (Mi>maior) { maior = Mi; } 
    }
  fprintf(stderr, "maior = %f", maior);
  for (i = 0; i < n; i++) 
    { 
      int Gi = (int)(255*M[i]/maior + 0.5);
      if (Gi > 255) { Gi = 255; }
      G[i] = Gi;
    }

  /*esta imagem � uma foto da transformada de Hough */
  salvaimagem(5, G, NULL, NULL, "f8.pgm", Lmax, 180, 255);
  free(G);


  /* pega os nlinhas maiores valores */

  lista=novalista(nlinhas);
  for (i=0; i<n; i++)
    {
      if ( (M[i]>lista->cmin) || (lista->tam<lista->tammax) )
	{
	  no = novono();
	  no->c=M[i];/* a lista � ordenada por valores decrescentes de c */
	  no->v=-90+linha(i,Lmax);/*teta*/
	  no->h=coluna(i,Lmax);/*L*/
	  insereno(no, lista);
	}
    }

  mm=rn_alloc(2*nlinhas +1);
  mm[0] = 2*nlinhas;
  no=lista->primno;
  for (i=1; i<=lista->tam; i++)
    {
      mm[i]=(double)no->v; /*teta*/
      mm[i+nlinhas]=(double)no->h; /* L */
      no=no->prox;
    }
 liberalista(lista);

  return mm;
}

#define DEBUG_JANELAS 1

/* Fun��o de Detector Local de Linhas */
/* atribui uma nota a cada pixel, guardando o resultado em G */
/* G � uma c�pia de R, na qual a intensidade de cada pixel significa o quanto este pixel se parece com uma linha */
void detecta_linhas_usando_pesos(
  int *R,          /* Imagem real (um canal de cor) */
  int largura,     /*  */
  int altura,      /*  */
  int *G,          /* Imagem com valor maximo do prod. escalar (saida). */
  int *T,          /* Imagem com o teta correspondente (saida). */
  double *Tindice    /* cont�m os tetas T, s� que n�o adaptada para sair numa imagem */
)
{
  int i=0, j=0;/*contadores*/
  double *ci=NULL, *jp=NULL;/* curva ideal e janela de pesos */
  double **janelas=NULL;
  double sigma_pesos = ((double)JRAIO)/2.5;
  double max_prod; /* Valor m�ximo do prod. escalar */

  /* chutando fator sigma */
  ci = curvaideal(0, SIGMA, R, largura);
  jp = janela_de_pesos(sigma_pesos,R,largura);
  janelas = cria_janelas_ideais(NTETAS,ci,jp);
  max_prod = calcula_max_prod_janelas(NTETAS,janelas,jp);
  
  printf("max produto escalar = %10.4f\n", max_prod);

  if (DEBUG_JANELAS)
    { grava_janelas_ideais(NTETAS, janelas, jp); }

  /*Para cada pixel da imagem grande*/
  for (i=0; i<altura; i++)
    { 
      fprintf(stderr, ".");
      for (j=0; j<largura; j++)
        { 
	  double melhor_prod; int melhor_i;
          double Gv, Tv;
          if ((i<JRAIO) || (i>=altura-JRAIO) || (j<JRAIO) || (j>=largura-JRAIO))
            { /* Pixel muito prximo aa borda, ignore: */
              melhor_prod = 0.0;
              melhor_i = 0;
            }
          else
            { 
              calcula_melhor_teta(
                R, largura, i, j, jp, janelas, NTETAS, 
                &melhor_prod, &melhor_i
              );
              /* melhor_prod � o melhor prod. esc. p/ todos os angulos teta neste pixel */
              /* melhor_prod >> 0 se o casamento foi forte. */
              /* melhor_prod << 0 se o casamento foi forte mas oposto. */
              /* melhor_prod = 0 se nao tem nada a ver. */
              /* Gv = 1.0 se casamento perfeito, 0.0 se p�ssimo. */
            }
          Gv = 0.5*fabs(melhor_prod)/max_prod + 0.5;
          if (Gv < 0.0) { Gv = 0.0; }
          if (Gv > 1.0) { Gv = 1.0; }
          Tv = ((double)melhor_i)/((double)NTETAS - 1);
          G[indice(i,j,largura)] = (int)(255 * Gv + 0.5);
          /* �ndice da m�scara de melhor casamento: */
          if (Tindice) { Tindice[indice(i,j,largura)] = (double)melhor_i; }
          /* Tindice deve j� vir alocado. Se ele vier NULL ent�o despreza esta matriz */
          /* este fato � usado para rodar o programa s� com o detector local, sem extrair as linhas */
          T[indice(i,j,largura)] = (int)(255 * Tv + 0.5);
          /*este if eh soh para testes */
/* 	  if ((i > 75) && (i < 80) && (j > 45) && (j < 80)) */
/* 	    { fprintf(stderr, "(%f)", melhor_prod); } */
        }
    }

  libera_janelas_ideais(NTETAS, janelas);
  if (janelas) { free(janelas); }

}


/* usa uma gaussiana unidimensional */
void detectalinhas(int *R, int largura, int altura, int *G)
{
  int i=0, j=0;/*contadores*/
  double teta=0;
  double *ci;

  /* Na figura grande(f0a.ppm) o ponto (306,173) � uma linha*/
  /* (129,144) � reta na vertical(angulo 87 graus) na figura f0.ppm */
  /* chutando fator o */
  ci = curvaideal(0, SIGMA, R, largura);

  /*Para cada pixel da imagem grande*/
  for (i=JRAIO; i<(altura-JRAIO); i++)
    { 
      fprintf(stderr, ".");
      for (j=JRAIO; j<(largura-JRAIO); j++)
        { 
	  double cmax = -1e30, Gv;
          for (teta=-90; teta<90; teta+=5)
            {
	      double cteta = cos(teta/180.0*M_PI);
              double steta = sin(teta/180.0*M_PI);
              double c;
              c=calculacurvacasamento(R,largura,i,j,cteta,steta,ci);
              if (c > cmax) { cmax = c; }
            }
	  /*este if eh soh para testes */
	  /*          if ((i > 100) && (i < 150) && (j > 100) && (j < 150))
		      { fprintf(stderr, "(%f)", cmax); }*/
          /* cmax � a melhor nota para todos os angulos teta neste pixel */
          /* cmax = 0 se o casamento foi perfeito. */
          /* cmax << 0 se o casamento foi ruim. */
          /* Gv = 1.0 se casamento perfeito, 0.0 se p�ssimo. */
          Gv = (15 - fabs(cmax))/15.0;
          if (Gv < 0.0) { Gv = 0.0; }
          if (Gv > 1.0) { Gv = 1.0; }
          G[indice(i,j,largura)] = (int)(255 * Gv + 0.5);
        }
    }

  if (ci) free(ci);

}


/* G � a imagem depois de aplicado o detector local;
   T � a mesma imagem com o melhor teta encontrado para cada pixel como intensidade do pixel;
   H � o resultado da transformada de Hough;
   largura_H e altura_H s�o referentes a matriz H
   
   O elemento HH[i,j] eh o escore da linha com distancia L e inclinacao
   teta, onde 
   
     teta =  180*((i+0.5)/altura_H) + TETAMIN
     L    = Lmax*(2*(j+0.5)/largura_H - 1)
     
   onde Lmax = sqrt(altura^2 + largura^2).
  */
int gera_hough(int *R, double *T, int largura, int altura, int **HH, int *largura_H, int *altura_H)
{
  int i=0, j=0; /* contadores */
  int *H=NULL;
  double *Htemp=NULL;
  double maxcor=0;

  double Lmax = sqrt(largura*largura + altura*altura);
  /* maior dist�ncia poss�vel */

  int largurah = (int)ceil(Lmax);
  /* arbitr�rio; nao vale a pena ser maior que 2*Lmax. */
  int alturah=NTETAS/2; 
  /* arbitr�rio; nao vale a pena ser mior que NTETAS. */

  int tam = largurah*alturah;

  /*  [valor - valorminimo] / [ (valormaximo-valorminimo)/(elementos-1) ] */

  if ( (!HH) || (!R) ||(!T) || (!largura_H) || (!altura_H) )
    { printf("Par�metros incorretos em gera_hough.\n"); return -1; }

  if (*HH) free(*HH);

  H = (int *)notnull(calloc(tam, sizeof(int)), "no mem");
  Htemp = rn_alloc(tam);

  *HH = H;
  *largura_H = largurah;
  *altura_H = alturah;

  for (i=0; i<tam; i++) { H[i]=0; Htemp[i]=0.0; }

  /* passa nas imagens de R e T e soma os valores de R nos buckets corretos de H */
  for (i=0; i<altura; i++)
    for (j=0; j<largura; j++)
      {
	int indr = indice(i,j,largura);
	double Tindice = T[indr];
        double teta = TETAMIN + 180.0*(Tindice/((double)NTETAS));
        double peso = ((double)R[indr] - 128);
	double dist = distancia(i,j,cos_graus(teta),sin_graus(teta));
	/* melhoria a fazer: estes senos e cossenos j� foram calculados em outro lugar */

        /* Calcula indices na transformada de Hough para (teta,dist): */
	int lin_H = (int)(((double)alturah)*(teta-TETAMIN)/180.0);
	int col_H = (int)(((double)largurah)*(dist+Lmax)/(Lmax+Lmax));
	/* pode n�o parecer, mas lin_H depende de NTETAS, pois alturah = NTETAS */

	int indh = indice(lin_H,col_H,largurah);
	/* criei estas vari�veis pra ficar mais f�cil de debugar */
        
	Htemp[indh] += peso;
	if (maxcor<Htemp[indh]) { maxcor = Htemp[indh]; }
	/* fa�o isso porque a acumula��o poderia estar gerando valores maiores que 255 */

	if (lin_H == 0)
	  {
	    printf("Peso = %f  teta = %f   pixel = [%d,%d]\n", peso, teta, i, j);
          }
      }

  for (i=0; i<tam; i++) { H[i] = (int)((Htemp[i]/maxcor)*255.0 -0.5); }
  if (Htemp) free(Htemp);

  /* Pronto! H cont�m o resultado de minha transformada de Hough adaptada.
     Ainda n�o t� bom... */

  return (i*j);
}


/* n � o n�mero de linhas a detectar;
   R � a imagem original;
   G guarda o quanto cada pixel � parecido com uma linha;
   T guarda o melhor �ngulo encontrado para uma linha com centro neste pixel;
   largura e altura da imagem.
*/
double *extrai_linhas_com_pesos(int n, int *R, int *G, int *T, int largura, int altura, char *pref_out)
{
  int i=0, j=0;/* contadores */
  no_t *no=NULL;
  double *m=NULL;
  double *Tindice=NULL; /* guarda os valores */

  lista_t *lista=NULL;

  int *H=NULL; /* matriz depois da transformada de Hough adaptada */
  int largura_H=0, altura_H=0;

  int Lmax = sqrt(largura*largura+altura*altura); /* maior dist�ncia poss�vel */
  /* isso t� feio. est� definido dentro de gera hough, e � usado aqui de novo, ent�o tem que repetir a f�rmula */

  Tindice = rn_alloc(largura*altura);
  /* Cada elemento desta matriz cont�m o �ndice da m�scara que melhor casou para este mesmo pixel na imagem original. */
  /* para se obter o �ngulo a partir deste �ndice:
       �ngulo = 180*indice/NTETAS + TETAMIN ]
  */

  detecta_linhas_usando_pesos(R,largura,altura,G,T,Tindice);

  gera_hough(G,Tindice,largura,altura, &H, &largura_H, &altura_H);
  /* H cont�m a transformada de Hough */
  /* altura_H e largura_H s�o relativos � matriz H */

  fprintf(stderr, "Matriz de acumula��o %i x %i \n", largura_H, altura_H);

  /* Essa matriz n�o t� boa n�o, mas pelo menos ela j� existe... */


  /**********************************************************

     DEBUGAR ESTA ....

  ********************************************************/

  /* o pr�ximo passo seria extrair os picos dessa matriz.
     verificar se a fun��o extraipicos serve */
  /* n�o . a fun�o extrai picos n�o serve... */

  lista=novalista(TAMMAX);

  /*Para cada elemento de H*/
  for (i=0; i<altura_H; i++) /* as linhas s�o os tetas */
    {
      for (j=0; j<largura_H; j++) /* as colunas s�o as dist�ncias */
        {
	  no=novono();
	  no->c = H[indice(i,j,largura_H)];/* double */ /* ordena pelo menor valor de c */
	  /* um valor grande de H[ij] significa muitos pixels parecidos com linha para �ngulos e dist�ncias semelhantes. */
	  no->v = i;/* v e h(int) guardam as coordenadas vertical(teta) e horizontal(distancia) do bucket */
	  no->h = j;
	  insereno(no, lista);
	}
    }

  salvaimagem(5, H, NULL, NULL, txtcat(pref_out, "-hough.pgm"), largura_H, altura_H, 255);
  /* foto da transformada de Hough */

  if (H) free(H);

  if (m) free(m);
  m=rn_alloc(3*n+1);
  m[0] = n; /* n�mero de elementos v�lidos */
  no=lista->primno;
  for(i=1; i<=n; i++)
    {
      /* Converte indices v,h da transformada de Hough para teta,dist: */
      double teta = TETAMIN + 180.0*((double)no->v + 0.5)/((double)altura_H);
      double dist = - Lmax + (Lmax+Lmax)*((double)no->h + 0.5)/((double)largura_H);
      double peso = no->c;
      m[i] = teta;
      m[i+n] = dist;
      m[i+2*n] = peso;
      no=no->prox;
      //      fprintf(stderr, "[%d,%d] teta = %f  dist = %f \n", no->v, no->h, m[i], m[i+n]);
    }

  liberalista(lista);
  
  return m;
}


int *extrailinhas(int n, int *R, int largura, int altura)
{
  int i=0, j=0;/*contadores*/
  double teta=0;
  no_t *no=NULL;
  lista_t *lista=NULL;
  int *m;
  double *ci;

  /* Na figura grande(f0a.ppm) o ponto (306,173) � uma linha*/
  /* (129,144) � reta na vertical(angulo 87 graus) na figura f0.ppm */
  /* chutando o fator sigma */
  ci = curvaideal(0, SIGMA, R, largura);

  lista=novalista(TAMMAX);

  /*Para cada pixel da imagem grande*/
  for (i=JRAIO; i<(altura-JRAIO); i++)
    {
      for (j=JRAIO; j<(largura-JRAIO); j++)
        {
          for (teta=-90; teta<90; teta+=5)
            { 
	      double cteta = cos(teta/180.0*M_PI);
              double steta = sin(teta/180.0*M_PI);
              no=novono();
              no->teta=teta;
              no->v=i;
              no->h=j;
              no->c=calculacurvacasamento(R,largura,i,j,cteta,steta,ci);
              insereno(no, lista);
            }
        }
    }
  printf("Inicia o Bucketing.\n");
  m=bucketing(lista, largura, altura);
  printf("Termina o Bucketing.\n");
  liberalista(lista);
  printf("Antes de extrair picos.\n");
  lista = extraipicos(n,m,largura,altura);
  printf("Termina extra��o de picos.\n");
  free(m);
  m=(int *)notnull(calloc(2*n, sizeof(int)), "no mem");
  no=lista->primno;
  for(i=0; i<2*n; i++)
    {
      m[i]=no->v;
      m[2*i]=no->teta;
      no=no->prox;
    }
  printf("Terminou(quase).\n");
  liberalista(lista);
  
  return m;
}


/******************************************************************

   TODO LIST :                    
     
     verificar se calcula_mascara_casamento quebra com teta < 0

******************************************************************/
