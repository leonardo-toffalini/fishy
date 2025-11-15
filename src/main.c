#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fishy.h"
#include "plot.h"
#include "utils.h"

double f1(double x) {
  return 1;
}

double f2(double x) {
  return -exp(x);
}

int main(void) {
  double a = 0.0, b = 1;
  double alpha = 1.0, beta = 2.0;
  int n = 100;
  double sol[n+2];
  double xs[n+2];

  solve_poisson1d(a, b, alpha, beta, n, f2, sol);
  // write_doubles("data.txt", sol, n);
  for (int i = 0; i < n + 2; i++)
    printf("%lf\n", sol[i]);

  double h = (b - a) / (n + 1);
  for (int i = 0; i < n + 2; i++)
    xs[i] = a + (i + 1) * h;

  plot_scatter(xs, sol, n);

  return 0;
}
