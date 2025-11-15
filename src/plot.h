#ifndef PLOT_H
#define PLOT_H

#include "../include/raylib.h"

#define HEIGHT 600
#define WIDTH 800

float translateX(float a, float b, float x, float pad) {
  return (x - a)/(b - a) * (WIDTH - 2 * pad) + pad;
}

float translateY(float c, float d, float y, float pad) {
  float t = (y - c)/(d - c) * (HEIGHT - 2 * pad) + pad;
  return HEIGHT - t;
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

void plot_scatter(double *xs, double *ys, int n) {
  float a = min(xs, n), b = max(xs, n);
  float c = min(ys, n), d = max(ys, n);
  int num_ticks = 10;
  float xticks[num_ticks];
  float yticks[num_ticks];
  float pad = 20.0f;
  char buf[100];

  float xstride = (b - a) / num_ticks;
  float ystride = (d - c) / num_ticks;
  for (int i = 0; i < num_ticks; i ++) {
    xticks[i] = a + (i + 1) * xstride;
    yticks[i] = c + (i + 1) * ystride;
  }

  InitWindow(WIDTH, HEIGHT, "Raylib Scatter");
  while (!WindowShouldClose()) {
    BeginDrawing();
    ClearBackground(BLACK);
    DrawLine(0, HEIGHT-pad, WIDTH, HEIGHT-pad, RAYWHITE);
    DrawLine(pad, 0, pad, HEIGHT, RAYWHITE);
    for (int i = 0; i < n; i++) {
      float screenx = translateX(a, b, xs[i], pad);
      float screeny = translateY(c, d, ys[i], pad);
      DrawCircle(
        screenx,
        screeny,
        1, BLUE
      );

    }

    for (int i = 0; i < num_ticks; i++) {
      float screenx = translateX(a, b, xticks[i], pad);
      float screeny = translateY(c, d, yticks[i], pad);
      sprintf(buf, "%.2f", xticks[i]);
      DrawText(buf, screenx - 8, HEIGHT-pad + 6, 12, RAYWHITE);
      DrawLine(screenx, HEIGHT-pad + 5, screenx, HEIGHT-pad - 5, RAYWHITE);

      sprintf(buf, "%.2f", yticks[i]);
      DrawText(buf, pad + 8, screeny - 6, 12, RAYWHITE);
      DrawLine(pad - 5, screeny, pad + 5, screeny, RAYWHITE);

    }

    EndDrawing();
  }
  CloseWindow();
}

#endif // !PLOT_H
