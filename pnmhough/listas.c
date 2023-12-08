#include <listas.h>
#include <affirm.h>

lista_t *novalista(int tammax)
{
  lista_t *lista;
  if ( (lista=(lista_t *)calloc(1, sizeof(lista_t))) == NULL )
    printf("Problemas ao alocar nova lista.\n");
  lista->tam=0;
  lista->tammax=tammax;
  lista->cmin=0;
  lista->primno=NULL;
  return lista;
}

no_t *novono()
{
  no_t *no;
  if ( (no = (no_t *)calloc(1, sizeof(no_t))) == NULL )
    printf("Problemas ao alocar novo nó.\n");
  no->c=no->teta=0;
  no->v=no->h=0;
  no->prox=no->ant=NULL;
  return no;
}

int liberano(no_t *no)
{
  if (no) { free(no); }
  return 1;
}

int libno(no_t *no)
{
  if (no)
    { libno(no->prox);
      free(no);
    }
  return 1;
}

int liberalista(lista_t *lista)
{
  if (lista)
    {
      ((lista->primno)->ant)->prox=NULL;/*o último recebe NULL*/
      libno(lista->primno);
      free(lista);
    }
  return 1;
}

/*insere o no na lista*/
/*o valor de no * resultante é inconsistente*/
int insereno(no_t *no, lista_t *lista)
{
  no_t *p;

  if (lista->tam==0)/*lista vazia*/
    {
      no->prox=no;
      no->ant=no;
      lista->primno=no;
      lista->cmin = ((lista->primno)->ant)->c;
      lista->tam++;/*aumenta já aqui*/
    }
  else/*lista não vazia*/
    {
      /*insere na posição correta*/       
      p=lista->primno;
      if (p->c < no->c)/*no começo da lista*/
        {
          no->prox=p;
          no->ant=p->ant;
          p->ant=no;
          (no->ant)->prox=no;
          lista->primno=no;
        }
      else
        {
          while( (p->prox!=lista->primno) && (no->c <= p->c) )
            p=p->prox;
          if (no->c > p->c)/*no meio da lista*/
            {
              no->ant=p->ant;
              p->ant=no;
              no->prox=p;
              (no->ant)->prox=no;
            }
          else/*no fim da lista*/
            {
              no->prox=p->prox;
              p->prox=no;
              no->ant=p;
              (lista->primno)->ant=no;
              lista->cmin=no->c;
            }
        }
      lista->tam++;/*aumenta aqui no final do else*/
    }
  /*se ultrapassar o tammax tira o último*/
  if(lista->tam > lista->tammax)
    {
      p=(lista->primno)->ant;
      (lista->primno)->ant=p->ant;
      (p->ant)->prox=lista->primno;
      lista->cmin=(p->ant)->c;
      liberano(p);
      lista->tam--;
    }

  return (lista->tam);
}
