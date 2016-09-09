#include "common.h"
#include "debug.h"
#include "gundersen.h"
#include "matrix.h"
#include "tensor.h"
#include "vector.h"

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
  tensor_ecrs_t *ptecrs = tensor_ecrs_malloc(n, m);
  tensor_ecrs_init(ptecrs);
  tensor_ecrs_compress(ptecrs, pc);
  tensor_ecrs_print(ptecrs);
  dbgprintf("\n");
  double *pv = vector_malloc(n);
  vector_init(pv, n, 1);
  matrix_t *pm = matrix_malloc(n, n);
  tensor_ecrs_mult(pm, ptecrs, pv, n);
  matrix_print(pm);
  dbgprintf("\n");
  int ndiag = tensor_ecrs_max_row(ptecrs);
  dbgprintf("ndiag: %d\n", ndiag);
  matrix_free(pm);
  tensor_ecrs_free(ptecrs);
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
