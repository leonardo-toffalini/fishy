#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fishy.h"
#include "plot.h"

double f1(double x) {
  return 1;
}

double f2(double x) {
  return -exp(x);
}

int main(void) {
  double a = 0.0, b = 1;
  double alpha = 1.0, beta = 2.0;
  int n = 20;
  double sol[n+2];
  double xs[n+2];

  solve_poisson1d(a, b, alpha, beta, n, f2, sol);

  double h = (b - a) / (n + 1);
  for (int i = 0; i < n + 2; i++)
    xs[i] = a + (i + 1) * h;

  // plot(xs, sol, n + 2, ".-", ORANGE);

  double ys[20 * 30];
  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 30; j++) {
      ys[IDX(i, j, 30)] = (double)(i + j);
    }
  }

  imshow(ys, 20, 30);

  return 0;
}
