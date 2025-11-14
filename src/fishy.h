#ifndef FISHY_H
#define FISHY_H

#define IDX(i, j, ldm) (i * ldm + j)

#include <stdio.h>

typedef double (*RHSFunc1D)(double x);
typedef double (*RHSFunc2D)(double x, double y);

typedef struct {
  double upper;
  double diag;
  double lower;
} TridiagMat;

typedef struct {
  TridiagMat upper;
  TridiagMat diag;
  TridiagMat lower;
} BlockTridiagMat;

void solve_tridiag_jacobi(TridiagMat A_h, double *rhs_values, int n, double *sol);
void solve_tridiag_gs(TridiagMat A_h, double *rhs_values, int n, double *sol);
void solve_blocktridiag_gs(BlockTridiagMat A_h, double *rhs_values, int n, double *sol);

void solve_poisson1d(double a, double b, int n, RHSFunc1D f_rhs, double *sol) {
  double h = (b - a) / (n + 1);
  double *rhs_values = malloc(n * sizeof(double));

  for (int i = 0; i < n; i++) {
    double x = a + (i + 1) * h;
    rhs_values[i] = f_rhs(x);
  }

  double inv_h_sq = 1 / (h*h);
  TridiagMat A_h = {
    .upper = -1 * inv_h_sq,
    .diag = 2 * inv_h_sq,
    .lower = -1 * inv_h_sq,
  };

  solve_tridiag_gs(A_h, rhs_values, n, sol);

  free(rhs_values);
}

void solve_poisson2d(double a, double b, double c, double d, int n, RHSFunc2D f_rhs, double *sol) {
  double h1 = (b - a) / (n + 1);
  double h2 = (d - c) / (n + 1);
  double *rhs_values = malloc(n * n * sizeof(double));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double x = a + (i + 1) * h1;
      double y = c + (j + 1) * h2;
      rhs_values[IDX(i, j, n)] = f_rhs(x, y);
    }
  }

  double inv_h1_sq = 1 / (h1*h1);
  double inv_h2_sq = 1 / (h2*h2);
  TridiagMat B_tilde = {
    .upper = -1 * inv_h1_sq,
    .diag = 2 * inv_h1_sq + 2 * inv_h2_sq,
    .lower = -1 * inv_h1_sq,
  };

  TridiagMat I_h = {
    .upper = 0.0,
    .diag = -1 * inv_h2_sq,
    .lower = 0.0,
  };

  BlockTridiagMat A_h = {
    .upper = I_h,
    .diag = B_tilde,
    .lower = I_h,
  };

  solve_blocktridiag_gs(A_h, rhs_values, n, sol);

  free(rhs_values);
}

void solve_tridiag_jacobi(TridiagMat A_h, double *rhs_values, int n, double *sol) {
  const int MAX_ITER = 20;
  double *x_new = malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) sol[i] = 0.0;

  double inv_diag = 1 / A_h.diag;

  for (int iter = 0; iter < MAX_ITER; iter++) {
    // x_(k+1) = D^-1 (b - (L + U)x_k)
    for (int i = 0; i < n; i++) {
      double left = (i == 0) ? 0.0 : sol[i-1];
      double right = (i == n - 1) ? 0.0 : sol[i+1];
      x_new[i] = inv_diag * (rhs_values[i] - A_h.lower * left - A_h.upper * right);
    }
    for (int i = 0; i < n; i++) sol[i] = x_new[i];
  }
  free(x_new);
}

void solve_tridiag_gs(TridiagMat A_h, double *rhs_values, int n, double *sol) {
  const int MAX_ITER = 20;
  for (int i = 0; i < n; i++) sol[i] = 0.0;

  double inv_diag = 1 / A_h.diag;

  for (int iter = 0; iter < MAX_ITER; iter++) {
    // x_(k+1) = (L + D)^-1 (b - U x_k)
    for (int i = 0; i < n; i++) {
      double left = (i == 0) ? 0.0 : sol[i-1];
      double right = (i == n - 1) ? 0.0 : sol[i+1];
      sol[i] = inv_diag * (rhs_values[i] - A_h.lower * left - A_h.upper * right);
    }
  }
}

void solve_blocktridiag_gs(BlockTridiagMat A_h, double *rhs_values, int n, double *sol) {
  const int MAX_ITER = 5000;
  for (int k = 0; k < n * n; k++) sol[k] = 0.0;

  double inv_diag = 1 / A_h.diag.diag;

  for (int iter = 0; iter < MAX_ITER; iter++) {
    // Gauss-Seidel over the 2D grid (lexicographic: j outer, i inner)
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
        int idx = IDX(i, j, n);
        double left  = (i == 0)      ? 0.0 : sol[IDX(i - 1, j, n)];
        double right = (i == n - 1)  ? 0.0 : sol[IDX(i + 1, j, n)];
        double down  = (j == 0)      ? 0.0 : sol[IDX(i, j - 1, n)];
        double up    = (j == n - 1)  ? 0.0 : sol[IDX(i, j + 1, n)];

        sol[idx] = inv_diag * (
          rhs_values[idx]
          - A_h.diag.lower * left
          - A_h.diag.upper * right
          - A_h.lower.diag * down
          - A_h.upper.diag * up
        );
      }
    }
  }
}

#endif // !FISHY_H
