#include "common.h"
#include "debug.h"
#include "matrix.h"
#include "tensor.h"
#include "vector.h"

double rand2() {
  return rand() / (double) RAND_MAX;
}

tensor_ecrs_t* tensor_ecrs_malloc(int n, int m) {
  tensor_ecrs_t *pt = (tensor_ecrs_t*) malloc(sizeof(tensor_ecrs_t));
  pt->R    = (int*) malloc(n*sizeof(int)+1);
  pt->CK   = (int*) malloc(m*sizeof(int)+1);
  pt->V    = (int*) malloc(m*sizeof(int)+1);
  pt->n    = n;
  pt->nnz  = 0;
  return pt;
}

void tensor_ecrs_init(tensor_ecrs_t *pt) {
  for (int i = 0; i < pt->n + 1; i++) {
    pt->R[i] = 0;
  }
  for (int i = 0; i < pt->nnz; i++) {
    pt->CK[i] = 0;
    pt->V[i]  = 0;
  }
}

void tensor_ecrs_free(tensor_ecrs_t *pt) {
  free(pt->R);
  free(pt->CK);
  free(pt->V);
  free(pt);
}

void tensor_ecrs_compress(tensor_ecrs_t *pt, cord_t *pc) {
  // For each non-zero, use the row number to fill in the R
  // buckets. Two adjacent buckets in R represent a range in CK and V
  // that have values in that row.
  pt->nnz  = pc->nnz;
  for (int i = 0; i < pt->nnz; i++) {
    int j = pc->c[i].t+1;
    pt->R[j] = i+1;
  }
  // If there is a gap in the row structure, between two adjacent
  // buckets, r0 and r1, fill in the buckets between them with the
  // same value as the r0, so that their range is zero, i.e. there are
  // no values in that row.
  pt->n = pc->n+1;
  for (int i = 1; i < pt->n; i++) {
    if (pt->R[i] < pt->R[i-1]) {
      pt->R[i] = pt->R[i-1];
    }
  }
  // Encode the column and tube indecies in such a way that they can
  // be decoded later, but consume less space. While we're at it,
  // stash the value of the entry as well.
  for (int i = 0; i < pc->nnz; i++) {
    pt->CK[i] = pc->c[i].c * pc->n + pc->c[i].r;
    pt->V[i]  = pc->c[i].v;
  }
}

void tensor_ecrs_mult(matrix_t *pm, tensor_ecrs_t *pt, double *v, int n) {
  int c, t, k, r0;
  for (int r = 1; r < pt->n; r++) {
    r0 = r-1;
    for (int i = pt->R[r0]; i < pt->R[r]; i++) {
      c = pt->CK[i] / n;
      t = pt->CK[i] % n;
      pm->M[r0*n + c] += v[t] * pt->V[i];
    }
  }
}

int tensor_ecrs_max_row(tensor_ecrs_t *pt) {
  int ndiag = 0;
  for (int r = 1; r < pt->n; r++) {
    ndiag = fmax(ndiag, pt->R[r] - pt->R[r-1]);
  }
  return ndiag;
}

void tensor_ecrs_print(tensor_ecrs_t *pt) {
  for (int i = 0; i < pt->n; i++) {
    dbgprintf(" %4d", i);
  }
  dbgprintf("\n");
  for (int i = 0; i < pt->n; i++) {
    dbgprintf("%4d ", pt->R[i]);
  }
  dbgprintf("\n");
  for (int i = 0; i < pt->nnz; i++) {
    dbgprintf("%4d ", pt->CK[i]);
  }
  dbgprintf("\n");
  for (int i = 0; i < pt->nnz; i++) {
    dbgprintf("%4d ", pt->V[i]);
  }
  dbgprintf("\n");
}

tensor_eellpack_t* tensor_eellpack_malloc(int n, int ndiag) {
  tensor_eellpack_t *pt = (tensor_eellpack_t*) malloc(sizeof(tensor_ecrs_t));
  pt->A     = matrix_malloc(n, ndiag);
  pt->B     = matrix_malloc(n, ndiag);
  pt->n     = n;
  pt->ndiag = ndiag;
  pt->nnz   = 0;
  return pt;
}

void tensor_eellpack_init(tensor_eellpack_t *pt) {
  /*
c
c fill coef with zero elements and jcoef with row numbers.------------ 
c
      do 4 j=1,ndiag 
         do 41 i=1,nrow
            coef(i,j) = 0.0d0
            jcoef(i,j) = i
 41      continue
 4    continue
  */
  int i, j, n, ndiag;
  n     = pt->n;
  ndiag = pt->ndiag;
  for (i = 0; i < n; i++) {
    for (j = 0; j < ndiag; j++) {
      pt->A->M[i*n+j] = 0.0;
      pt->B->M[i*n+j] = i;
      //dbgprintf("%d\n", j);
    }
  }
}

void tensor_eellpack_print(tensor_eellpack_t *pt) {
  for (int i = 0; i < pt->ndiag; i++) {
    dbgprintf(" %4d", i);
  }
  dbgprintf("\n");
  matrix_print(pt->A);
  dbgprintf("\n");
  matrix_print(pt->B);
  dbgprintf("\n");
}

void tensor_eellpack_free(tensor_eellpack_t *pt) {
  matrix_free(pt->A);
  matrix_free(pt->B);
  free(pt);
}

void tensor_ecrs_to_eellpack(tensor_ecrs_t *ptecrs, tensor_eellpack_t *pteell) {
}

cord_t* cord_malloc(int n, int m) {
  cord_t *pc = (cord_t*) malloc(sizeof(cord_t));
  pc->c   = (c_t*) malloc(m*sizeof(c_t));
  pc->n   = n;
  pc->nnz = 0;
  return pc;
}

void cord_free(cord_t *pc) {
  free(pc->c);
  free(pc);
}

void cord_gen(cord_t *pc, double entpb) {
  pc->nnz = 0;
  int n = pc->n;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        if (rand2() < entpb) {
          pc->c[pc->nnz].r = i;
          pc->c[pc->nnz].c = j;
          pc->c[pc->nnz].t = k;
          pc->c[pc->nnz].v = 1;
          pc->nnz++;
        }
      }
    }
  }
}

int cord_sort_cmp_ecrs(const void *p1, const void *p2) {
  int result;
  c_t *pa = (c_t*) p1;
  c_t *pb = (c_t*) p2;
  if ((result = pa->t - pb->t) == 0) {
    if ((result = pa->r - pb->r) == 0) {
      result = pa->c - pb->c;
    }
  }
  return result;
}

void cord_sort_ecrs(cord_t *pc) {
  qsort(pc->c, pc->nnz, sizeof(c_t), cord_sort_cmp_ecrs);
}

void cord_print(cord_t *pc) {
  for (int i = 0; i < pc->nnz; i++) {
    dbgprintf("(%4d, %4d, %4d) = %4d\n",
              pc->c[i].t, pc->c[i].r,
              pc->c[i].c, pc->c[i].v);
  }
}

void cord_print_stats(cord_t *pc) {
  int n = pc->n;
  int m = n*n*n;
  dbgprintf("%d/%d = %lf\n", pc->nnz, m, pc->nnz/(double)m);
}
