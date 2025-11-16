#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fishy.h"
#include "plot.h"

float f1(float x) {
  return 1.0f;
}

float f2(float x) {
  return -exp(x);
}

int main(void) {
  // float a = 0.0f, b = 1.0f;
  // float alpha = 1.0f, beta = 2.0f;
  // int n = 20;
  // float sol[n+2];
  // float xs[n+2];
  //
  // solve_poisson1d(a, b, alpha, beta, n, f2, sol);
  //
  // float h = (b - a) / (n + 1);
  // for (int i = 0; i < n + 2; i++)
  //   xs[i] = a + (i + 1) * h;

  // plot(xs, sol, n + 2, ".-", ORANGE);

  int n = 30, m = 40;

  float ys[n * m];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      ys[IDX(i, j, m)] = (float)(i + j);
    }
  }

  plot_surface(ys, n, m);

  return 0;
}
