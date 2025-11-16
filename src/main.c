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
  float a = 0.0f, b = 1.0f;
  float alpha = 1.0f, beta = 2.0f;
  int n = 20;
  float sol[n+2];
  float xs[n+2];

  solve_poisson1d(a, b, alpha, beta, n, f2, sol);

  float h = (b - a) / (n + 1);
  for (int i = 0; i < n + 2; i++)
    xs[i] = a + (i + 1) * h;

  // plot(xs, sol, n + 2, ".-", ORANGE);

  float ys[20 * 30];
  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 30; j++) {
      ys[IDX(i, j, 30)] = (float)(i + j);
    }
  }

  imshow(ys, 20, 30);

  return 0;
}
