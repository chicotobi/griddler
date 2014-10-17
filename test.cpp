/******************************************************************************
* FILE: omp_mm.c
* DESCRIPTION:
*   OpenMp Example - Matrix Multiply - C Version
*   Demonstrates a matrix multiply using OpenMP. Threads share row iterations
*   according to a predefined chunk size.
* AUTHOR: Blaise Barney
* LAST REVISED: 06/28/05
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define NRA 100                 /* number of rows in matrix A */
#define NCA 100                 /* number of columns in matrix A */
#define NCB 100                  /* number of columns in matrix B */

int main (int argc, char *argv[])
{
  int	tid, nthreads, i, j, k, l, chunk;
  double	a[NRA][NCA],           /* matrix A to be multiplied */
      b[NCA][NCB],           /* matrix B to be multiplied */
      c[NRA][NCB];           /* result matrix C */

  chunk = 5;                    /* set loop iteration chunk size */

  double d[NCA];


  /*** Initialize matrices ***/
  for (i=0; i<NRA; i++)
    for (j=0; j<NCA; j++)
      a[i][j]= i+j;
  for (i=0; i<NCA; i++)
    for (j=0; j<NCB; j++)
      b[i][j]= i*j;
  for (i=0; i<NRA; i++)
    for (j=0; j<NCB; j++)
      c[i][j]= 0;
  for (i=0; i<NCA; i++)
    d[i]= i;

  /*** Spawn a parallel region explicitly scoping all variables ***/
#pragma omp parallel shared(a,b,c,d,nthreads,chunk) private(i,j,k)
  {
    for(l=0;l<10000;l++)
#pragma omp for schedule (static, chunk)
      for (i=0; i<NRA; i++)
        for(j=0; j<NCB; j++)
          for (k=0; k<NCA; k++)
            c[i][j] += a[i][k] * d[k] * b[k][j];
  }


}
