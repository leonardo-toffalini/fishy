#ifndef FISHY_H
#define FISHY_H

#define IDX(i, j, ldm) ((i) * (ldm) + (j))

typedef float (*RHSFunc1D)(float x);
typedef float (*RHSFunc2D)(float x, float y);

typedef struct {
  float upper;
  float diag;
  float lower;
} TridiagMat;

typedef struct {
  TridiagMat upper;
  TridiagMat diag;
  TridiagMat lower;
} BlockTridiagMat;

void solve_tridiag_jacobi(TridiagMat A_h, float *rhs_values, int n, float *sol);
void solve_tridiag_gs(TridiagMat A_h, float *rhs_values, int n, float *sol);
void solve_blocktridiag_gs(BlockTridiagMat A_h, float *rhs_values, int n, float *sol);

void solve_poisson1d(float a, float b, float alpha, float beta, int n, RHSFunc1D f_rhs, float *sol) {
  float h = (b - a) / (n + 1);
  float *rhs_values = malloc((n + 2) * sizeof(float));
  rhs_values[0] = alpha;
  rhs_values[n+1] = beta;

  for (int i = 1; i < n + 1; i++) {
    float x = a + i * h;
    rhs_values[i] = f_rhs(x);
  }

  float inv_h_sq = 1 / (h*h);
  TridiagMat A_h = {
    .upper = -1 * inv_h_sq,
    .diag = 2 * inv_h_sq,
    .lower = -1 * inv_h_sq,
  };

  solve_tridiag_gs(A_h, rhs_values, n, sol);

  free(rhs_values);
}

void solve_poisson2d(float a, float b, float c, float d, int n, RHSFunc2D f_rhs, float *sol) {
  float h1 = (b - a) / (n + 1);
  float h2 = (d - c) / (n + 1);
  float *rhs_values = malloc(n * n * sizeof(float));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      float x = a + (i + 1) * h1;
      float y = c + (j + 1) * h2;
      rhs_values[IDX(i, j, n)] = f_rhs(x, y);
    }
  }

  float inv_h1_sq = 1 / (h1*h1);
  float inv_h2_sq = 1 / (h2*h2);
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

void solve_tridiag_jacobi(TridiagMat A_h, float *rhs_values, int n, float *sol) {
  const int MAX_ITER = 20;
  float *x_new = malloc(n * sizeof(float));
  for (int i = 0; i < n; i++) sol[i] = 0.0f;

  float inv_diag = 1 / A_h.diag;

  for (int iter = 0; iter < MAX_ITER; iter++) {
    // x_(k+1) = D^-1 (b - (L + U)x_k)
    for (int i = 0; i < n; i++) {
      float left = (i == 0) ? 0.0f : sol[i-1];
      float right = (i == n - 1) ? 0.0f : sol[i+1];
      x_new[i] = inv_diag * (rhs_values[i] - A_h.lower * left - A_h.upper * right);
    }
    for (int i = 0; i < n; i++) sol[i] = x_new[i];
  }
  free(x_new);
}

void solve_poisson2d_9pt(float a, float b, float c, float d, int n, RHSFunc2D f_rhs, float *sol) {
  float h1 = (b - a) / (n + 1);
  float h2 = (d - c) / (n + 1);
  float *rhs_values = malloc(n * n * sizeof(float));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      float x = a + (i + 1) * h1;
      float y = c + (j + 1) * h2;
      rhs_values[IDX(i, j, n)] = f_rhs(x, y);
    }
  }

  float inv_h1_sq = 1 / (h1*h1);
  float inv_h2_sq = 1 / (h2*h2);
  TridiagMat B_tilde = {
    .upper = 4.0f   / 6.0f * inv_h1_sq,
    .diag  = -20.0f / 6.0f * inv_h1_sq,
    .lower = 4.0f   / 6.0f * inv_h1_sq,
  };

  TridiagMat I_h = {
    .upper = 1.0f / 6.0f * inv_h2_sq,
    .diag  = 4.0f / 6.0f * inv_h2_sq,
    .lower = 1.0f / 6.0f * inv_h2_sq,
  };

  BlockTridiagMat A_h = {
    .upper = I_h,
    .diag = B_tilde,
    .lower = I_h,
  };

  solve_blocktridiag_gs(A_h, rhs_values, n, sol);

  free(rhs_values);
}

void solve_tridiag_gs(TridiagMat A_h, float *rhs_values, int n, float *sol) {
  const int MAX_ITER = 500;
  for (int i = 0; i < n+2; i++) sol[i] = 0.0;
  sol[0] = rhs_values[0];
  sol[n+1] = rhs_values[n+1];

  float inv_diag = 1 / A_h.diag;

  for (int iter = 0; iter < MAX_ITER; iter++) {
    // x_(k+1) = (L + D)^-1 (b - U x_k)
    for (int i = 1; i < n+1; i++) {
      float left = sol[i-1];
      float right = sol[i+1];
      sol[i] = inv_diag * (rhs_values[i] - A_h.lower * left - A_h.upper * right);
    }
  }
}

void solve_blocktridiag_gs(BlockTridiagMat A_h, float *rhs_values, int n, float *sol) {
  const int MAX_ITER = 500;
  for (int k = 0; k < n * n; k++) sol[k] = 0.0f;

  float inv_diag = 1 / A_h.diag.diag;

  for (int iter = 0; iter < MAX_ITER; iter++) {
    // Gauss-Seidel over the 2D grid (lexicographic: j outer, i inner)
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        int idx = IDX(i, j, n);
        float left     = (i == 0)      ? 0.0f : sol[IDX(i - 1, j, n)];
        float right    = (i == n - 1)  ? 0.0f : sol[IDX(i + 1, j, n)];
        float down     = (j == 0)      ? 0.0f : sol[IDX(i, j - 1, n)];
        float up       = (j == n - 1)  ? 0.0f : sol[IDX(i, j + 1, n)];
        float topleft  = (i == 0 || j == 0)         ? 0.0f : sol[IDX(i-1, j-1, n)];
        float topright = (i == 0 || j == n - 1)     ? 0.0f : sol[IDX(i-1, j+1, n)];
        float botleft  = (i == n - 1 || j == 0)     ? 0.0f : sol[IDX(i+1, j-1, n)];
        float botright = (i == n - 1 || j == n - 1) ? 0.0f : sol[IDX(i+1, j+1, n)];

        sol[idx] = inv_diag * (
          rhs_values[idx]
          - A_h.diag.lower * left
          - A_h.diag.upper * right
          - A_h.lower.diag * down
          - A_h.upper.diag * up

          - A_h.upper.lower * topleft
          - A_h.upper.upper * topright
          - A_h.lower.lower * botleft
          - A_h.lower.upper * botright
        );
      }
    }
  }
}

#endif // !FISHY_H
