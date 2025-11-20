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

float f6(float x, float y) {
  return 2 * PI * PI * f5(x, y);
}

float f7(float x, float y) {
    return sin(PI * x) * cos(PI * y);
}

float g1(float x, float y) {
  return 0.0f;
}

float g2(float x, float y) {
  return (1 / (2 * PI * PI)) * sin(PI * x) * cos(PI * y);
}

float exact1(float x, float y) {
  return (1 / (2 * PI * PI)) * sin(PI * x) * cos(PI * y);
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
  float a = 0.0f, b = 1.0f, c = 0.0f, d = 1.0f;
  int n = 200, m = 200;
  float ys[(n + 2) * (m + 2)];

  float exact_ys[(n+2) * (m + 2)];
  float h1 = (b - a) / (n + 1);
  float h2 = (d - c) / (n + 1);
  for (int i = 0; i < n + 2; i++) {
    for (int j = 0; j < n + 2; j++) {
      float x = a + i * h1;
      float y = c + j * h2;
      exact_ys[IDX(i, j, n+2)] = exact1(x, y);
    }
  }

  solve_poisson2d(a, b, c, d, n, f7, g2, ys, NINE_POINT_STENCIL);

  float diff[(n+2) * (m+2)];
  get_diff(ys, exact_ys, (n+2) * (m+2), diff);
  float err = max_norm_error(ys, exact_ys, (n+2) * (n+2));
  printf("max norm error = %f\n", err);

  plot_surface(ys, n+2, m+2);
  // plot_surface(exact_ys, n+2, m+2);
  // plot_surface(diff, n+2, m+2);
}

int main(void) {
  demo_possion2d();

  return 0;
}
