#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fishy.h"
#include "../include/raylib.h"

double f(double x) {
  return -sin(3.141592653589793 * x);
}

float translate(float a, float b, float x, float pad) {
  return (x - a)/(b - a) * (800.0f - 2 * pad) + pad;
}

double max(double *x, int n) {
  double res = x[0];
  for (int i = 0; i < n; i++)
    res = x[i] > res ? x[i] : res;

  return res;
}

double min(double *x, int n) {
  double res = x[0];
  for (int i = 0; i < n; i++)
    res = x[i] < res ? x[i] : res;

  return res;
}

void plot(double *xs, double *ys, int n) {
  float a = min(xs, n), b = max(xs, n);
  float c = min(ys, n), d = max(ys, n);
  float pad = 20.0f;

  InitWindow(800, 600, "Raylib Plot");
  while (!WindowShouldClose()) {
    for (int i = 0; i < n; i++)
      DrawCircle(
        translate(a, b, xs[i], pad),
        translate(c, d, ys[i], pad),
        1, RAYWHITE
      );

    EndDrawing();
  }
  CloseWindow();
}

int main(void) {
  double a = 0.0, b = 1.0;
  int n = 100;
  double sol[n];
  double xs[n];

  solve_poisson1d(a, b, n, f, sol);
  for (int i = 0; i < n; i++)
    printf("%lf ", sol[i]);
  printf("\n");

  double h = (b - a) / (n + 1);
  for (int i = 0; i < n; i++)
    xs[i] = a + (i + 1) * h;

  plot(xs, sol, n);

  return 0;
}
