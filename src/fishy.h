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

void solve_tridiag_gs(TridiagMat A_h, float *rhs_values, int n, float *sol);
void solve_blocktridiag_gs(BlockTridiagMat A_h, float *rhs_values, int n, float *sol);
void solve_helmholtz2d(float a, float b, float c, float d, float lam, int n, RHSFunc2D f_rhs, float *sol);

BlockTridiagMat block_tridiag_from_kernel(float *kernel) {
  // kernel must be an array of size 9, with a row major matrix layour
  // e.g. [[00, 01, 02], [10, 11, 12], [20, 21, 22]]

  TridiagMat Mid = {
    .upper = kernel[IDX(1, 2, 3)],
    .diag  = kernel[IDX(1, 1, 3)],
    .lower = kernel[IDX(1, 0, 3)],
  };

  TridiagMat Up = {
    .upper = kernel[IDX(0, 2, 3)],
    .diag  = kernel[IDX(0, 1, 3)],
    .lower = kernel[IDX(0, 0, 3)],
  };

  TridiagMat Down = {
    .upper = kernel[IDX(2, 2, 3)],
    .diag  = kernel[IDX(2, 1, 3)],
    .lower = kernel[IDX(2, 0, 3)],
  };

  BlockTridiagMat A_h = {
    .upper = Up,
    .diag  = Mid,
    .lower = Down,
  };

  return A_h;
}

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
    .upper = -1.0f * inv_h_sq,
    .diag  =  2.0f * inv_h_sq,
    .lower = -1.0f * inv_h_sq,
  };

  solve_tridiag_gs(A_h, rhs_values, n, sol);

  free(rhs_values);
}

void solve_poisson2d(float a, float b, float c, float d, int n, RHSFunc2D f_rhs, float *sol) {
  // solve an equation of the form: -Delta u = f
  solve_helmholtz2d(a, b, c, d, 0.0f, n, f_rhs, sol);
}

void solve_helmholtz2d(float a, float b, float c, float d, float lam, int n, RHSFunc2D f_rhs, float *sol) {
  // solve an equation of the form: -Delta u + lamda * u = f
  float h1 = (b - a) / (n + 1);
  float h2 = (d - c) / (n + 1);
  float *rhs_values = malloc((n + 2) * (n + 2) * sizeof(float));
  for (int i = 0; i < n + 2; i++) {
    rhs_values[IDX(i,   0, n+2)] = 0.0f; // first col
    rhs_values[IDX(i, n+1, n+2)] = 0.0f; // last col
    rhs_values[IDX(0,   i, n+2)] = 0.0f; // first row
    rhs_values[IDX(n+1, i, n+2)] = 0.0f; // last row
  }

  for (int i = 1; i < n + 1; i++) {
    for (int j = 1; j < n + 1; j++) {
      float x = a + i * h1;
      float y = c + j * h2;
      rhs_values[IDX(i, j, n+2)] = f_rhs(x, y);
    }
  }

  float inv_h1_sq = 1 / (h1*h1);
  float inv_h2_sq = 1 / (h2*h2);
  float kernel[9] = {
                 0.0f,                         -1.0f * inv_h2_sq,             0.0f,
    -1.0f * inv_h1_sq, lam + 2.0f * inv_h1_sq + 2.0f * inv_h2_sq, -1.0f * inv_h1_sq,
                 0.0f,                         -1.0f * inv_h2_sq,             0.0f,
  };

  BlockTridiagMat A_h = block_tridiag_from_kernel(kernel);

  solve_blocktridiag_gs(A_h, rhs_values, n, sol);

  free(rhs_values);
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
  float kernel[9] = {
    -1.0f/6.0f * inv_h2_sq, -4.0f/6.0f  * inv_h2_sq, -1.0f/6.0f * inv_h2_sq,
    -4.0f/6.0f * inv_h1_sq,  20.0f/6.0f * inv_h1_sq, -4.0f/6.0f * inv_h1_sq,
    -1.0f/6.0f * inv_h2_sq, -4.0f/6.0f  * inv_h2_sq, -1.0f/6.0f * inv_h2_sq,
  };

  BlockTridiagMat A_h = block_tridiag_from_kernel(kernel);

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
  for (int i = 1; i < n + 1; i++)
    for (int j = 1; j < n + 1; j++)
      sol[IDX(i, j, n+2)] = 0.0f;

  float inv_diag = 1 / A_h.diag.diag;

  for (int iter = 0; iter < MAX_ITER; iter++) {
    // Gauss-Seidel over the 2D grid (lexicographic: j outer, i inner)
    for (int i = 1; i < n+1; i++) {
      for (int j = 1; j < n+1; j++) {
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
