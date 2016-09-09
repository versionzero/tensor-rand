
#ifndef _matrix_
#define _matrix_

typedef struct {
  int *M;
  int n, m;
} matrix_t;

matrix_t* matrix_malloc(int n, int m);
void matrix_free(matrix_t *pm);
void matrix_print(matrix_t *pt);

#endif
