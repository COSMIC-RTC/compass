#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cublas_v2.h"

#define IDX2C(i,j,ld) (((j)*(ld))+(i))

void modify (float *m, int ldm, int n, int p, int q, float alpha,
	     float beta)
{
  cublasSscal (n-p, alpha, &m[IDX2C(p,q,ldm)], ldm);
  cublasSscal (ldm-p, beta, &m[IDX2C(p,q,ldm)], 1);
}
#define M 6 //# of columns idx i
#define N 5 //# of rows idx j

int idx_test(float *a)
{
  int i, j;
  cublasStatus stat;
  float* devPtrA;
  //a = (float *)malloc (M * N * sizeof (*a));

  for (j = 0; j < N; j++) {
    for (i = 0; i < M; i++) {
      a[IDX2C(i,j,M)] = i * M + j + 1;
    }
  }

  for (j = 0; j < N; j++) {
    // for each row
    for (i = 0; i < M; i++) {
      // print each column
      printf ("%7.0f (%d,%d)", a[IDX2C(i,j,M)],i,j);
    }
    printf ("\n");
  }
  /*
      1 (0,0)      7 (1,0)     13 (2,0)     19 (3,0)     25 (4,0)     31 (5,0)
      2 (0,1)      8 (1,1)     14 (2,1)     20 (3,1)     26 (4,1)     32 (5,1)
      3 (0,2)      9 (1,2)     15 (2,2)     21 (3,2)     27 (4,2)     33 (5,2)
      4 (0,3)     10 (1,3)     16 (2,3)     22 (3,3)     28 (4,3)     34 (5,3)
      5 (0,4)     11 (1,4)     17 (2,4)     23 (3,4)     29 (4,4)     35 (5,4)
  */

  printf ("\n");

  cublasInit();
  stat = cublasAlloc (M*N, sizeof(*a), (void**)&devPtrA);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf ("device memory allocation failed");
    cublasShutdown();
    return EXIT_FAILURE;
  }
  stat = cublasSetMatrix (M, N, sizeof(*a), a, M, devPtrA, M);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf ("data download failed");
    cublasFree (devPtrA);
    cublasShutdown();
    return EXIT_FAILURE;
  }
  modify (devPtrA, M, N, 1, 2, 16.0f, 12.0f);
  stat = cublasGetMatrix (M, N, sizeof(*a), devPtrA, M, a, M);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf ("data upload failed");
    cublasFree (devPtrA);
    cublasShutdown();
    return EXIT_FAILURE;
  }
  cublasFree (devPtrA);
  cublasShutdown();
  for (j = 0; j < N; j++) {
    for (i = 0; i < M; i++) {
      printf ("%7.0f", a[IDX2C(i,j,M)]);
    }
    printf ("\n");
  }
  return EXIT_SUCCESS;
}

/* in yorick :
a=array(0.0f,6*5);idx_test,&a;
b=reform(a,[2,6,5]);
"";
b(,1);b(,2);b(,3);b(,4);b(,5);



 */
