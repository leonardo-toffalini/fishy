Color hotcold[256] ;

typedef struct {
  float r, g, b;
} FloatColor;

Color get_color(double v, double vmin, double vmax) {
  FloatColor c = {1.0, 1.0, 1.0};
  double dv = vmax - vmin;
  v = fmin(fmax(v, vmin), vmax);

  if (v < (vmin + 0.25 * dv)) {
     c.r = 0;
     c.g = 4 * (v - vmin) / dv;
  } else if (v < (vmin + 0.5 * dv)) {
     c.r = 0;
     c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
  } else if (v < (vmin + 0.75 * dv)) {
     c.r = 4 * (v - vmin - 0.5 * dv) / dv;
     c.b = 0;
  } else {
     c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
     c.b = 0;
  }

  Color res = {0};
  res.r = 255.0 * c.r;
  res.g = 255.0 * c.g;
  res.b = 255.0 * c.b;
  res.a = 255;

  return res;
}

