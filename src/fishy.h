#ifndef FISHY_H
#define FISHY_H

#include <stdio.h>

typedef double (*RHSFunc)(double x);

typedef struct {
  double upper;
  double diag;
  double lower;
} TridiagMat;

void solve_tridiag_jacobi(TridiagMat A_h, double *rhs_values, int n, double *sol);
void solve_tridiag_gs(TridiagMat A_h, double *rhs_values, int n, double *sol);

void solve_poisson1d(double a, double b, int n, RHSFunc f_rhs, double *sol) {
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

#endif // !FISHY_H
