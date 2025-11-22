#ifndef FISHY_H
#define FISHY_H

#define IDX(i, j, ldm) ((i) * (ldm) + (j))

typedef enum {
  FIVE_POINT_STENCIL = 0,
  NINE_POINT_STENCIL,
} StencilType;

typedef float (*RHSFunc1D)(float x);
typedef float (*RHSFunc2D)(float x, float y);
typedef float (*BCFunc2D)(float x, float y);

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

typedef struct {
  float a, b, c, d, lam;
  int n;
  RHSFunc2D f_rhs;
  BCFunc2D bc_func;
  float *sol;
  StencilType stencil;
} FishyParams;

void solve_tridiag_gs(TridiagMat A_h, float *rhs_values, int n, float *sol);
void solve_blocktridiag_gs(BlockTridiagMat A_h, float *rhs_values, int n, float *sol);
void solve_helmholtz2d(FishyParams params);

float max_norm_error(float *ys, float *exact_ys, int n) {
  float err = 0.0f;
  for (int i = 0; i < n; i++)
    err = fmax(fabs(exact_ys[i] - ys[i]), err);

  return err;
}

void get_diff(float *ys, float *exact_ys, int n, float *diff) {
  for (int i = 0; i < n; i++)
    diff[i] = exact_ys[i] - ys[i];
}

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

void get_kernel_for_stencil(StencilType stencil, float h1, float h2, float lam, float *kernel) {
  float inv_h1_sq = 1 / (h1*h1);
  float inv_h2_sq = 1 / (h2*h2);

  switch (stencil) {
    case FIVE_POINT_STENCIL:
        kernel[0] =              0.0f; kernel[1] =                         -1.0f * inv_h2_sq; kernel[2] =              0.0f;
        kernel[3] = -1.0f * inv_h1_sq; kernel[4] = lam + 2.0f * inv_h1_sq + 2.0f * inv_h2_sq; kernel[5] = -1.0f * inv_h1_sq;
        kernel[6] =              0.0f; kernel[7] =                         -1.0f * inv_h2_sq; kernel[8] =              0.0f;
      break;
    case NINE_POINT_STENCIL:
    default:
      kernel[0] = -1.0f/6.0f * inv_h2_sq; kernel[1] = -4.0f/6.0f  * inv_h2_sq; kernel[2] = -1.0f/6.0f * inv_h2_sq;
      kernel[3] = -4.0f/6.0f * inv_h1_sq; kernel[4] =  20.0f/6.0f * inv_h1_sq; kernel[5] = -4.0f/6.0f * inv_h1_sq;
      kernel[6] = -1.0f/6.0f * inv_h2_sq; kernel[7] = -4.0f/6.0f  * inv_h2_sq; kernel[8] = -1.0f/6.0f * inv_h2_sq;
      break;
  };
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

void solve_poisson2d(float a, float b, float c, float d, int n, RHSFunc2D f_rhs, BCFunc2D bc_func, float *sol, StencilType stencil) {
  // solve an equation of the form: -Delta u = f
  solve_helmholtz2d((FishyParams){a, b, c, d, 0.0f, n, f_rhs, bc_func, sol, stencil});
}

void solve_helmholtz2d(FishyParams params) {
  // solve an equation of the form: -Delta u + lamda * u = f
  float a = params.a; float b = params.b; float c = params.c; float d = params.d;
  int n = params.n;
  float h1 = (b - a) / (n + 1);
  float h2 = (d - c) / (n + 1);
  float *rhs_values = malloc((n + 2) * (n + 2) * sizeof(float));
  for (int i = 0; i < n + 2; i++) {
    float x = a + i * h1;
    float y = c + i * h2;
    rhs_values[IDX(i,   0, n+2)] = params.bc_func(x, c); // y=c (first col/bottom)
    rhs_values[IDX(i, n+1, n+2)] = params.bc_func(x, d); // y=d (last col/top)
    rhs_values[IDX(0,   i, n+2)] = params.bc_func(a, y); // x=a (first row/left)
    rhs_values[IDX(n+1, i, n+2)] = params.bc_func(b, y); // x=b (last row/right)
  }

  for (int i = 1; i < n + 1; i++) {
    for (int j = 1; j < n + 1; j++) {
      float x = a + i * h1;
      float y = c + j * h2;
      rhs_values[IDX(i, j, n+2)] = params.f_rhs(x, y);
    }
  }

  float kernel[9];
  get_kernel_for_stencil(params.stencil, h1, h2, params.lam, kernel);
  BlockTridiagMat A_h = block_tridiag_from_kernel(kernel);

  solve_blocktridiag_gs(A_h, rhs_values, n, params.sol);

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
  const int MAX_ITER = 5000;
  for (int k = 0; k < (n + 2) * (n + 2); k++) sol[k] = 0.0f;
  for (int i = 0; i < n + 2; i++) {
    sol[IDX(i,   0, n+2)] = rhs_values[IDX(i,   0, n+2)];
    sol[IDX(i, n+1, n+2)] = rhs_values[IDX(i, n+1, n+2)];
    sol[IDX(0,   i, n+2)] = rhs_values[IDX(0,   i, n+2)];
    sol[IDX(n+1, i, n+2)] = rhs_values[IDX(n+1, i, n+2)];
  }

  float inv_diag = 1 / A_h.diag.diag;

  for (int iter = 0; iter < MAX_ITER; iter++) {
    for (int i = 1; i < n + 1; i++) {
      for (int j = 1; j < n + 1; j++) {
        int idx = IDX(i, j, n + 2);
        float left     = sol[IDX(i - 1, j, n+2)];
        float right    = sol[IDX(i + 1, j, n+2)];
        float down     = sol[IDX(i, j - 1, n+2)];
        float up       = sol[IDX(i, j + 1, n+2)];
        float topleft  = sol[IDX(i-1, j-1, n+2)];
        float topright = sol[IDX(i-1, j+1, n+2)];
        float botleft  = sol[IDX(i+1, j-1, n+2)];
        float botright = sol[IDX(i+1, j+1, n+2)];

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
