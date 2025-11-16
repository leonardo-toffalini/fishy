#ifndef PLOT_H
#define PLOT_H

#include "../include/raylib.h"
#include "colormap.h"
#include <string.h>

#define HEIGHT 900
#define WIDTH 1600

float translateX(float a, float b, float x, float pad) {
  return (x - a)/(b - a) * (WIDTH - 2 * pad) + pad;
}

float translateY(float c, float d, float y, float pad) {
  float t = (y - c)/(d - c) * (HEIGHT - 2 * pad) + pad;
  return HEIGHT - t;
}

float max(float *x, int n) {
  float res = x[0];
  for (int i = 0; i < n; i++)
    res = x[i] > res ? x[i] : res;

  return res;
}

float min(float *x, int n) {
  float res = x[0];
  for (int i = 0; i < n; i++)
    res = x[i] < res ? x[i] : res;

  return res;
}

void plot(float *xs, float *ys, int n, char *line_style, Color col) {
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

  SetConfigFlags(FLAG_MSAA_4X_HINT);

  InitWindow(WIDTH, HEIGHT, "Raylib Scatter");
  while (!WindowShouldClose()) {
    BeginDrawing();
    ClearBackground(BLACK);
    DrawLine(0, HEIGHT-pad, WIDTH, HEIGHT-pad, RAYWHITE);
    DrawLine(pad, 0, pad, HEIGHT, RAYWHITE);
    for (int i = 0; i < n; i++) {
      float screenx1 = translateX(a, b, xs[i], pad);
      float screeny1 = translateY(c, d, ys[i], pad);

      if ((strcmp(line_style, "-") == 0 || strcmp(line_style, ".-") == 0) && i + 1 < n) {
        float screenx2 = translateX(a, b, xs[i+1], pad);
        float screeny2 = translateY(c, d, ys[i+1], pad);

        DrawLine(screenx1, screeny1, screenx2, screeny2, col);
      }

      if (strcmp(line_style, ".") == 0 || strcmp(line_style, ".-") == 0) {
        DrawCircle(
          screenx1,
          screeny1,
          2, col
        );
      }

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

void imshow(float *ys, int n, int m) {
  float ymin = min(ys, n * m), ymax = max(ys, n * m);
  float pad = 40.0f;
  float dx = (WIDTH - 2 * pad) / (float)m;
  float dy = (HEIGHT - 2 * pad) / (float)n;
  printf("dx = %f\n", dx);
  printf("dy = %f\n", dy);
  char buf[100];

  SetConfigFlags(FLAG_MSAA_4X_HINT);

  InitWindow(WIDTH, HEIGHT, "Raylib Heatmap");
  while (!WindowShouldClose()) {
    BeginDrawing();
    ClearBackground((Color){224, 217, 199});
    DrawLine(0, HEIGHT-pad, WIDTH, HEIGHT-pad, BLACK);
    DrawLine(pad, 0, pad, HEIGHT, BLACK);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        float screenx = translateX(0, m, j, pad);
        float screeny = translateY(0, n, i, pad);
        Color c = get_color(ys[IDX(i, j, m)], (float)ymin, (float)ymax);
        DrawRectangle(screenx, screeny - dy, dx, dy, c);
        DrawRectangleLines(screenx, screeny - dy, dx, dy, BLACK);
        // DrawCircle(screenx, screeny, 5, c);
      }
    }

    for (int j = 1; j < m + 1; j++) {
      float screenx = translateX(0, m, j, pad);
      sprintf(buf, "%d", j);
      DrawText(buf, screenx - 8, HEIGHT-pad + 6, 12, BLACK);
      DrawLine(screenx, HEIGHT-pad + 5, screenx, HEIGHT-pad - 5, BLACK);
    }
    for (int i = 1; i < n + 1; i++) {
      float screeny = translateY(0, n, i, pad);
      sprintf(buf, "%d", i);
      DrawText(buf, pad - 22, screeny - 6, 12, BLACK);
      DrawLine(pad - 5, screeny, pad + 5, screeny, BLACK);
    }

    EndDrawing();
  }
  CloseWindow();
}

void plot_surface(float *ys, int n, int m) {

}

#endif // !PLOT_H
