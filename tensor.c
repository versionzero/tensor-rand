#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

double rand2() {
  return rand() / (double) RAND_MAX;
}

struct tensor_t {
  int *R, *CK, *V;
  int N, nnz;
};

struct matrix_t {
  int *pM;
  int N, M;
};

struct c_t {
  int r, c, t, v;
};
  
struct cord_t {
  struct c_t *c;
  int nnz, N;
};

struct matrix_t* matrix_malloc(int n, int m) {
  struct matrix_t *pm = (struct matrix_t*) malloc(sizeof(struct matrix_t));
  pm->pM = (int*) malloc(n*m*sizeof(int)+1);
  pm->N  = n;
  pm->M  = m;
  return pm;
}

void matrix_free(struct matrix_t *pm) {
  free(pm->pM);
  free(pm);
}

void matrix_print(struct matrix_t *pt) {
  int k;
  for (int i = 0; i < pt->N; i++) {
    for (int j = 0; j < pt->M; j++) {
      k = i * pt->N + j;
      printf("%4d ", pt->pM[k]);
    }
    printf("\n");
  }
  printf("\n");
}

double* vector_malloc(int n) {
  double *pv = malloc(n*sizeof(double));
  return pv;
}

void vector_init(double *pv, int n) {
  for (int i = 0; i < n; i++) {
    pv[i] = 1;
  }
}

struct tensor_t* tensor_malloc(int n, int m) {
  struct tensor_t *pt = (struct tensor_t*) malloc(sizeof(struct tensor_t));
  pt->R    = (int*) malloc(n*sizeof(int)+1);
  pt->CK   = (int*) malloc(m*sizeof(int)+1);
  pt->V    = (int*) malloc(m*sizeof(int)+1);
  pt->N    = n;
  pt->nnz  = 0;
  return pt;
}

void tensor_init(struct tensor_t *pt) {
  for (int i = 0; i < pt->N; i++) {
    pt->R[i] = 0;
  }
  for (int i = 0; i < pt->nnz; i++) {
    pt->CK[i] = 0;
    pt->V[i]  = 0;
  }
}

void tensor_free(struct tensor_t *pt) {
  free(pt->R);
  free(pt->CK);
  free(pt->V);
  free(pt);
}

void tensor_compress(struct tensor_t *pt, struct cord_t *pc) {
  // For each non-zero, use the row number to fill in the R
  // buckets. Two adjacent buckets in R represent a range in CK and V
  // that have values in that row.
  pt->nnz  = pc->nnz;
  for (int i = 0; i < pt->nnz; i++) {
    int j = pc->c[i].r+1;
    pt->R[j] = i+1;
  }
  // If there is a gap in the row structure, between two adjacent
  // buckets, r0 and r1, fill in the buckets between them with the
  // same value as the r0, so that their range is zero, i.e. there are
  // no values in that row.
  pt->N = pc->N+1;
  for (int i = 1; i < pt->N; i++) {
    if (pt->R[i] < pt->R[i-1]) {
      pt->R[i] = pt->R[i-1];
    }
  }
  // Encode the column and tube indecies in such a way that they can
  // be decoded later, but consume less space. While we're at it,
  // stash the value of the entry as well.
  for (int i = 0; i < pc->nnz; i++) {
    pt->CK[i] = pc->c[i].c * pc->N + pc->c[i].t;
    pt->V[i]  = pc->c[i].v;
  }  
}

void tensor_mult(struct matrix_t *pm, struct tensor_t *pt, double *v, int n) {
  int c, t, k, r0;  
  for (int r = 1; r < pt->N; r++) {
    r0 = r-1;
    for (int i = pt->R[r0]; i < pt->R[r]; i++) {
      c = pt->CK[i] / n;
      t = pt->CK[i] % n;
      k = pt->R[r];
      printf("(%4d, %4d, %4d)\n", k, c, t);
    }
  }
}

void tensor_print(struct tensor_t *pt) {
  for (int i = 0; i < pt->N; i++) {
    printf("%4d ", i);
  }
  printf("\n");
  for (int i = 0; i < pt->N; i++) {
    printf("%4d ", pt->R[i]);
  }
  printf("\n");
  for (int i = 0; i < pt->nnz; i++) {
    printf("%4d ", pt->CK[i]);
  }
  printf("\n");
  for (int i = 0; i < pt->nnz; i++) {
    printf("%4d ", pt->V[i]);
  }
  printf("\n");
}

struct cord_t* cord_malloc(int n, int m) {
  struct cord_t *pc = (struct cord_t*) malloc(sizeof(struct cord_t));
  pc->c   = (struct c_t*) malloc(m*sizeof(struct c_t));
  pc->N   = n;
  pc->nnz = 0;
  return pc;
}

void cord_free(struct cord_t *pc) {
  free(pc->c);
  free(pc);
}

void cord_gen(struct cord_t *pc, double entpb) {
  srand(time(NULL));
  pc->nnz = 0;
  int n = pc->N;
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

int cord_sort_cmp(const void *p1, const void *p2) {
  int result;
  struct c_t *pa = (struct c_t*) p1;
  struct c_t *pb = (struct c_t*) p2;
  if ((result = pa->r - pb->r) == 0) {
    if ((result = pa->c - pb->c) == 0) {
      result = pa->t - pb->t;
    }
  }
  return result;
}

void cord_sort(struct cord_t *pc) {
  qsort(pc->c, pc->nnz, sizeof(struct c_t), cord_sort_cmp);
}

void cord_print(struct cord_t *pc) {
  for (int i = 0; i < pc->nnz; i++) {
    printf("(%4d, %4d, %4d) = %4d\n",
           pc->c[i].r, pc->c[i].c,
           pc->c[i].t, pc->c[i].v);
  }
}

void cord_print_stats(struct cord_t *pc) {
  int n = pc->N;
  int m = n*n*n;
  printf("%d/%d = %lf\n", pc->nnz, m, pc->nnz/(double)m);
}


void run(int n, double entpb) {
  // Since we pick generated coordicates based on entpb, there may be
  // cases where we get sligtly more entries than were requested. As a
  // result, we we set an upper bound of 2*entpb to adjust for
  // this. The alternative is to allowcate n*n*n entries, since we
  // know it will never be the case that that baoud will be exceded.
  int m = (int) ceil(entpb*2*n*n*n);
  struct cord_t *pc = cord_malloc(n, m);
  cord_gen(pc, entpb);
  cord_sort(pc);
  cord_print(pc);
  printf("\n");
  cord_print_stats(pc);
  printf("\n");
  struct tensor_t *pt = tensor_malloc(n, m);
  tensor_init(pt);
  tensor_compress(pt, pc);
  tensor_print(pt);
  printf("\n");
  double *pv = vector_malloc(n);
  vector_init(pv, n);
  struct matrix_t *pm = matrix_malloc(n, n);
  tensor_mult(pm, pt, pv, n);
  matrix_print(pm);
  printf("\n");
  matrix_free(pm);
  tensor_free(pt);
  cord_free(pc);  
}

struct cord_t *cord_test(struct cord_t *pc) {
  int T[12] = { 0, 0, 0, 0, 0, 0, 1, 1, 1,  1,  1,  1 };
  int R[12] = { 0, 0, 1, 1, 2, 2, 0, 0, 1,  1,  2,  2 };
  int C[12] = { 1, 3, 0, 2, 1, 2, 0, 2, 1,  2,  0,  3 };
  int V[12] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
  pc->nnz = 0;
  int n = pc->N;
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
  if (argc < 1) {
    printf("usage: prog <n>");
  }
  int n = atoi(argv[1]);
  double entpb = 0.01;
  run(n, entpb);
}
