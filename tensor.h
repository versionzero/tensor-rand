
#ifndef _tensor_
#define _tensor_

typedef struct {
  int *R, *CK, *V;
  int n, nnz;
} tensor_ecrs_t;

typedef struct {
  matrix_t *A, *B;
  int n, nnz, ndiag;
} tensor_eellpack_t;

typedef struct {
  int r, c, t, v;
} c_t;
  
typedef struct {
  c_t *c;
  int nnz, n;
} cord_t;

tensor_ecrs_t* tensor_ecrs_malloc(int n, int m);
void tensor_ecrs_init(tensor_ecrs_t *pt);
void tensor_ecrs_free(tensor_ecrs_t *pt);
void tensor_ecrs_compress(tensor_ecrs_t *pt, cord_t *pc);
void tensor_ecrs_mult(matrix_t *pm, tensor_ecrs_t *pt, double *v, int n);
int tensor_ecrs_max_row(tensor_ecrs_t *pt);
void tensor_ecrs_print(tensor_ecrs_t *pt);

tensor_eellpack_t* tensor_eellpack_malloc(int n, int ndiag);
void tensor_eellpack_init(tensor_eellpack_t *pt);
void tensor_eellpack_print(tensor_eellpack_t *pt);
void tensor_eellpack_free(tensor_eellpack_t *pt);

void tensor_ecrs_to_eellpack(tensor_ecrs_t *ptecrs, tensor_eellpack_t *pteell);

cord_t* cord_malloc(int n, int m);
void cord_free(cord_t *pc);
void cord_gen(cord_t *pc, double entpb);
int cord_sort_cmp_ecrs(const void *p1, const void *p2);
void cord_sort_ecrs(cord_t *pc);
void cord_print(cord_t *pc);
void cord_print_stats(cord_t *pc);

#endif
