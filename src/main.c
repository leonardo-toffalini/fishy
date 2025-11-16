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

float f3(float x, float y) {
  return sin(PI * x) * sin(2 * PI * y);
}

float f4(float x, float y) {
  return 1.0f;
}

float f5(float x, float y) {
  return sin(PI * x) * sin(PI * y);
}

void demo_possion1d() {
  float a = 0.0f, b = 1.0f;
  float alpha = 1.0f, beta = 2.0f;
  int n = 20;
  float sol[n+2];
  float xs[n+2];

  solve_poisson1d(a, b, alpha, beta, n, f2, sol);

  float h = (b - a) / (n + 1);
  for (int i = 0; i < n + 2; i++)
    xs[i] = a + (i + 1) * h;

  plot(xs, sol, n + 2, ".-", BLUE);
}

void demo_possion2d() {
  float a = 0.0f, b = 1.0f, c = 0.0f, d = 4.0f;
  int n = 400, m = 400;
  float ys[n * m];

  solve_poisson2d(a, b, c, d, n, f5, ys);
  plot_surface(ys, n, m);
}

int main(void) {
  demo_possion2d();

  return 0;
}
