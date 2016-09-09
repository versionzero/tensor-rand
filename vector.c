#include "common.h"
#include "debug.h"
#include "matrix.h"
#include "tensor.h"
#include "vector.h"

double* vector_malloc(int n) {
  double *pv = malloc(n*sizeof(double));
  return pv;
}

void vector_init(double *pv, int n, int val) {
  for (int i = 0; i < n; i++) {
    pv[i] = val;
  }
}

void vector_free(double *pv) {
  free(pv);
}
