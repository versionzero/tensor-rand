#include "common.h"
#include "debug.h"
#include "matrix.h"
#include "tensor.h"
#include "vector.h"

matrix_t* matrix_malloc(int n, int m) {
  matrix_t *pm = (matrix_t*) malloc(sizeof(matrix_t));
  pm->M = (int*) malloc(n*m*sizeof(int)+1);
  pm->n = n;
  pm->m = m;
  return pm;
}

void matrix_free(matrix_t *pm) {
  free(pm->M);
  free(pm);
}

void matrix_print(matrix_t *pt) {
  int k;
  for (int i = 0; i < pt->n; i++) {
    for (int j = 0; j < pt->m; j++) {
      k = i * pt->n + j;
      dbgprintf("%4d ", pt->M[k]);
    }
    dbgprintf("\n");
  }
  dbgprintf("\n");
}
