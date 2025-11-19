#ifndef PLOT_H
#define PLOT_H

#include "../include/raylib.h"
#include "colormap.h"
#include <string.h>

#define HEIGHT 720
#define WIDTH 1280
#define BG_COLOR ((Color){224, 217, 199, 255})

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
    ClearBackground(BG_COLOR);
    DrawLineEx((Vector2){0, HEIGHT-pad}, (Vector2){WIDTH, HEIGHT-pad}, 2, BLACK);
    DrawLineEx((Vector2){pad, 0}, (Vector2){pad, HEIGHT}, 2, BLACK);
    for (int i = 0; i < n; i++) {
      float screenx1 = translateX(a, b, xs[i], pad);
      float screeny1 = translateY(c, d, ys[i], pad);

      if ((strcmp(line_style, "-") == 0 || strcmp(line_style, ".-") == 0) && i + 1 < n) {
        float screenx2 = translateX(a, b, xs[i+1], pad);
        float screeny2 = translateY(c, d, ys[i+1], pad);

        DrawLineEx((Vector2){screenx1, screeny1}, (Vector2){screenx2, screeny2}, 2, col);
      }

      if (strcmp(line_style, ".") == 0 || strcmp(line_style, ".-") == 0) {
        DrawCircle(
          screenx1,
          screeny1,
          4, col
        );
      }

    }

    for (int i = 0; i < num_ticks; i++) {
      float screenx = translateX(a, b, xticks[i], pad);
      float screeny = translateY(c, d, yticks[i], pad);
      sprintf(buf, "%.2f", xticks[i]);
      DrawText(buf, screenx - 8, HEIGHT-pad + 6, 12, BLACK);
      DrawLine(screenx, HEIGHT-pad + 5, screenx, HEIGHT-pad - 5, BLACK);

      sprintf(buf, "%.2f", yticks[i]);
      DrawText(buf, pad + 8, screeny - 6, 12, BLACK);
      DrawLine(pad - 5, screeny, pad + 5, screeny, BLACK);

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
    ClearBackground(BG_COLOR);
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

void normalize_float_array(float *ys, int n, int m) {
  float ymin = min(ys, n * m), ymax = max(ys, n * m);

  // avoid divide by 0
  if (ymin == ymax) {
    for (int i = 0; i < n * m; i++)
      ys[i] /= ys[i];
    return;
  }

  for (int i = 0; i < n * m; i++)
    ys[i] = (ys[i] - ymin) / (ymax - ymin);
}

void plot_surface(float *ys, int n, int m) {
  normalize_float_array(ys, n, m);

  SetConfigFlags(FLAG_MSAA_4X_HINT);

  InitWindow(WIDTH, HEIGHT, "Raylib Surface");

  Camera camera = {0};
  camera.position = (Vector3){18.0f, 21.0f, 18.0f};
  camera.target = (Vector3){0.0f, 0.0f, 0.0f};
  camera.up = (Vector3){0.0f, 1.0f, 0.0f};
  camera.fovy = 45.0f;
  camera.projection = CAMERA_PERSPECTIVE;

  Image img = {
    .data = ys,
    .width = m,
    .height = n,
    .format = PIXELFORMAT_UNCOMPRESSED_R32,
    .mipmaps = 1
  };

  Shader color_shader = LoadShader(NULL, "src/shaders/viridis.frag");

  float mesh_scale = fmin(16.0f / n, 16.0f / m);
  const float REF_TEXTURE_SIZE = 200.0f;

  float max_dim = fmax((float)n, (float)m);
  float fixed_tex_scale = REF_TEXTURE_SIZE / max_dim;

  Texture2D texture = LoadTextureFromImage(img);
  Mesh mesh = GenMeshHeightmap(img, (Vector3){n * mesh_scale, 24, m * mesh_scale});
  Model model = LoadModelFromMesh(mesh);
  for (int i = 0; i < model.materialCount; i++) {
    model.materials[i].shader = color_shader;
  }
  model.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = texture;
  Vector3 mapPosition = {-8.0f, 0.0f, -8.0f};

  CameraMode camera_mode = CAMERA_ORBITAL;

  DisableCursor();
  SetTargetFPS(60);
  while (!WindowShouldClose()) {
    if (IsKeyPressed(KEY_ENTER)) {
      camera.target = (Vector3){0.0f, 0.0f, 0.0f};
      if (camera_mode == CAMERA_ORBITAL) camera_mode = CAMERA_FREE;
      else camera_mode = CAMERA_ORBITAL;
      camera.target = (Vector3){0.0f, 0.0f, 0.0f};
    }
    UpdateCamera(&camera, camera_mode);

    BeginDrawing();
    ClearBackground(BG_COLOR);
    BeginMode3D(camera);
    DrawModel(model, mapPosition, 1.0f, WHITE);
    // DrawModelWires(model, mapPosition, 1.0f, WHITE);
    DrawGrid(20, 1.0f);
    EndMode3D();

    BeginShaderMode(color_shader);
    DrawTextureEx(texture, (Vector2){10, 10}, 0.0f, fixed_tex_scale, WHITE);
    EndShaderMode();

    EndDrawing();
  }

  UnloadTexture(texture);
  UnloadShader(color_shader);
  CloseWindow();

}

#endif // !PLOT_H
