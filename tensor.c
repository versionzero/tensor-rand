#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>

int debug = 1;

void dbgprintf(const char *format, ...) {
  if (debug == 1) {
    va_list arguments;
    va_start(arguments, format);
    vfprintf(stderr, format, arguments);
    va_end(arguments);
  }
}

double rand2() {
  return rand() / (double) RAND_MAX;
}

typedef struct {
  int *R, *CK, *V;
  int n, nnz;
} tensor_ecrs_t;

typedef struct {
  int *M;
  int n, m;
} matrix_t;

typedef struct {
  matrix_t *A, *B;
  int n, nnz;
} tensor_eellpack_t;

typedef struct {
  int r, c, t, v;
} c_t;
  
typedef struct {
  c_t *c;
  int nnz, n;
} cord_t;

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
  int length = 0;
  for (int r = 1; r < pt->n; r++) {
    length = fmax(length, pt->R[r] - pt->R[r-1]);
  }
  return length;
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

tensor_eellpack_t* tensor_eellpack_malloc(int n) {
  tensor_eellpack_t *pt = (tensor_eellpack_t*) malloc(sizeof(tensor_ecrs_t));
  pt->A   = matrix_malloc(n, n);
  pt->B   = matrix_malloc(n, n);
  pt->n   = n;
  pt->nnz = 0;
  return pt;
}

void tensor_eellpack_init(tensor_eellpack_t *pt) {
  
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


void run_ecrs(int n, double entpb) {
  // Since we pick generated coordicates based on entpb, there may be
  // cases where we get sligtly more entries than were requested. As a
  // result, we we set an upper bound of 2*entpb to adjust for
  // this. The alternative is to allowcate n*n*n entries, since we
  // know it will never be the case that that baoud will be exceded.
  int m = (int) ceil(entpb*2*n*n*n);
  cord_t *pc = cord_malloc(n, m);
  cord_gen(pc, entpb);
  cord_sort_ecrs(pc);
  cord_print(pc);
  dbgprintf("\n");
  cord_print_stats(pc);
  dbgprintf("\n");
  tensor_ecrs_t *pt = tensor_ecrs_malloc(n, m);
  tensor_ecrs_init(pt);
  tensor_ecrs_compress(pt, pc);
  tensor_ecrs_print(pt);
  dbgprintf("\n");
  double *pv = vector_malloc(n);
  vector_init(pv, n, 1);
  matrix_t *pm = matrix_malloc(n, n);
  tensor_ecrs_mult(pm, pt, pv, n);
  matrix_print(pm);
  dbgprintf("\n");
  int length = tensor_ecrs_max_row(pt);
  dbgprintf("length: %d\n", length);
  matrix_free(pm);
  tensor_ecrs_free(pt);
  cord_free(pc);  
}

cord_t *cord_test(cord_t *pc) {
  int T[12] = { 0, 0, 0, 0, 0, 0, 1, 1, 1,  1,  1,  1 };
  int R[12] = { 0, 0, 1, 1, 2, 2, 0, 0, 1,  1,  2,  2 };
  int C[12] = { 1, 3, 0, 2, 1, 2, 0, 2, 1,  2,  0,  3 };
  int V[12] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
  pc->nnz = 0;
  int n = pc->n;
  for (int i = 0; i < n; i++) {
    pc->c[pc->nnz].r = R[i];
    pc->c[pc->nnz].c = C[i];
    pc->c[pc->nnz].t = T[i];
    pc->c[pc->nnz].v = V[i];
    pc->nnz++;
  }
  return pc;
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    dbgprintf("usage: prog <n>");
    exit(1);
  }
  int seed = time(NULL);
  if (argc == 3) {
    seed = atoi(argv[2]);
  }
  dbgprintf("seed: %d\n", seed);
  srand(seed);
  int n = atoi(argv[1]);
  double entpb = 0.01;
  run_ecrs(n, entpb);
}
