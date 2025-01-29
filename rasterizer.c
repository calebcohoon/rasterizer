#include <conio.h>
#include <dos.h>
#include <math.h>
#include <stdlib.h>

#define INF 32767 // Infinity for 16 bit compilers
#define SW 320    // Screen width
#define SH 200    // Screen height
#define VIEWPORT_SIZE 1
#define PROJ_PLANE_Z 1.0f

#define SWAP_VECTOR2(a, b)                                                     \
  do {                                                                         \
    struct vector2 temp = (a);                                                 \
    (a) = (b);                                                                 \
    (b) = temp;                                                                \
  } while (0)

struct vector2 {
  int x;
  int y;
  float h;
};

struct vector3 {
  float x;
  float y;
  float z;
};

// Build an optimized palette for the main colors being used
void init_palette() {
  int color, shade;
  unsigned char index;
  float intensity;
  unsigned char base_colors[1][3] = {
      {0, 63, 0} // Green
  };

  for (color = 0; color < 1; color++) {
    for (shade = 0; shade < 256; shade++) {
      index = color * 256 + shade;
      intensity = shade / 255.0f;

      outp(0x3C8, index);

      // For the very last shade, just make it white
      // so we have a color for the background
      if (color == 0 && shade == 255) {
        outp(0x3C9, 63);
        outp(0x3C9, 63);
        outp(0x3C9, 63);
      } else {
        outp(0x3C9, (unsigned char)(base_colors[color][0] * intensity));
        outp(0x3C9, (unsigned char)(base_colors[color][1] * intensity));
        outp(0x3C9, (unsigned char)(base_colors[color][2] * intensity));
      }
    }
  }
}

unsigned char shade_color(unsigned char color, float intensity) {
  unsigned char shade;
  if (intensity > 1.0f)
    intensity = 1.0f;
  if (intensity < 0.0f)
    intensity = 0.0f;

  shade = (unsigned char)(intensity * 255.0f);
  return color * 256 + shade;
}

void set_mode(unsigned char mode) {
  union REGS regs;

  regs.w.ax = mode;

  int386(0x10, &regs, &regs);
}

void set_pixel(int x, int y, char color) {
  unsigned char *screen = (unsigned char *)0xA0000;
  int ax, ay;

  ax = SW / 2 + x;
  ay = SH / 2 - y;

  if (ax < 0 || ax >= SW || ay < 0 || ay >= SH) {
    return;
  }

  screen[ay * SW + ax] = color;
}

struct vector3 vec3_add(struct vector3 *a, struct vector3 *b) {
  struct vector3 result;

  result.x = a->x + b->x;
  result.y = a->y + b->y;
  result.z = a->z + b->z;

  return result;
}

struct vector3 vec3_sub(struct vector3 *a, struct vector3 *b) {
  struct vector3 result;

  result.x = b->x - a->x;
  result.y = b->y - a->y;
  result.z = b->z - a->z;

  return result;
}

struct vector3 vec3_mul(struct vector3 *v, float s) {
  struct vector3 result;

  result.x = v->x * s;
  result.y = v->y * s;
  result.z = v->z * s;

  return result;
}

struct vector3 vec3_div(struct vector3 *v, float d) {
  struct vector3 result;

  result.x = v->x / d;
  result.y = v->y / d;
  result.z = v->z / d;

  return result;
}

float vec3_len(struct vector3 *v) {
  return sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
}

float vec3_dot(struct vector3 *a, struct vector3 *b) {
  return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}

void interpolate(float i0, float d0, float i1, float d1, float *arry,
                 int *o_len) {
  float a, d, i;
  int idx = 0;

  if (i0 == i1) {
    arry[0] = d0;
    *o_len = 1;
    return;
  }

  a = (d1 - d0) / (float)(i1 - i0);
  d = d0;

  for (i = i0; i <= i1; i++) {
    arry[idx] = d;
    d = d + a;
    idx++;
  }

  *o_len = idx;
}

void concat_sides(float *x01, int x01_len, float *x12, int x12_len,
                  float *x012) {
  int i;

  for (i = 0; i <= x01_len - 1; i++) {
    x012[i] = x01[i];
  }

  for (i = 0; i <= x12_len; i++) {
    x012[i + (x01_len - 1)] = x12[i];
  }
}

void draw_line(struct vector2 *p0, struct vector2 *p1, unsigned char color) {
  struct vector2 p0_local = *p0;
  struct vector2 p1_local = *p1;
  int x, y, len;
  float slope_values[SW];
  float a;

  if (abs(p1_local.x - p0_local.x) > abs(p1_local.y - p0_local.y)) {
    if (p0_local.x > p1_local.x) {
      SWAP_VECTOR2(p0_local, p1_local);
    }

    interpolate(p0_local.x, p0_local.y, p1_local.x, p1_local.y, &slope_values,
                &len);

    for (x = p0_local.x; x <= p1_local.x; x++) {
      set_pixel(x, slope_values[x - p0_local.x], color);
    }
  } else {
    if (p0_local.y > p1_local.y) {
      SWAP_VECTOR2(p0_local, p1_local);
    }

    interpolate(p0_local.y, p0_local.x, p1_local.y, p1_local.x, &slope_values,
                &len);

    for (y = p0_local.y; y <= p1_local.y; y++) {
      set_pixel(slope_values[y - p0_local.y], y, color);
    }
  }
}

void draw_wireframe_triangle(struct vector2 *p0, struct vector2 *p1,
                             struct vector2 *p2, unsigned char color) {
  draw_line(p0, p1, color);
  draw_line(p1, p2, color);
  draw_line(p2, p0, color);
}

void draw_filled_triangle(struct vector2 *p0, struct vector2 *p1,
                          struct vector2 *p2, unsigned char color) {
  struct vector2 p0_local = *p0;
  struct vector2 p1_local = *p1;
  struct vector2 p2_local = *p2;
  int m, x, y;
  float *x_left, *x_right;
  float x01[SH], x12[SH], x02[SH], x012[SH];
  int x01_len, x12_len, x02_len;

  if (p1_local.y < p0_local.y) {
    SWAP_VECTOR2(p1_local, p0_local);
  }

  if (p2_local.y < p0_local.y) {
    SWAP_VECTOR2(p2_local, p0_local);
  }

  if (p2_local.y < p1_local.y) {
    SWAP_VECTOR2(p2_local, p1_local);
  }

  interpolate(p0_local.y, p0_local.x, p1_local.y, p1_local.x, &x01, &x01_len);
  interpolate(p1_local.y, p1_local.x, p2_local.y, p2_local.x, &x12, &x12_len);
  interpolate(p0_local.y, p0_local.x, p2_local.y, p2_local.x, &x02, &x02_len);

  concat_sides(&x01, x01_len, &x12, x12_len, &x012);

  m = x02_len / 2;

  if (x02[m] < x012[m]) {
    x_left = x02;
    x_right = x012;
  } else {
    x_left = x012;
    x_right = x02;
  }

  for (y = p0_local.y; y <= p2_local.y; y++) {
    for (x = x_left[y - p0_local.y]; x <= x_right[y - p0_local.y]; x++) {
      set_pixel(x, y, color);
    }
  }
}

void draw_shaded_triangle(struct vector2 *p0, struct vector2 *p1,
                          struct vector2 *p2, unsigned char color) {
  struct vector2 p0_local = *p0;
  struct vector2 p1_local = *p1;
  struct vector2 p2_local = *p2;
  int m, x, y;
  float *x_left, *x_right, *h_left, *h_right;
  float x01[SH], x12[SH], x02[SH], x012[SH];
  float h01[SH], h12[SH], h02[SH], h012[SH];
  int x01_len, x12_len, x02_len;
  int h01_len, h12_len, h02_len;

  if (p1_local.y < p0_local.y) {
    SWAP_VECTOR2(p1_local, p0_local);
  }

  if (p2_local.y < p0_local.y) {
    SWAP_VECTOR2(p2_local, p0_local);
  }

  if (p2_local.y < p1_local.y) {
    SWAP_VECTOR2(p2_local, p1_local);
  }

  interpolate(p0_local.y, p0_local.x, p1_local.y, p1_local.x, &x01, &x01_len);
  interpolate(p0_local.y, p0_local.h, p1_local.y, p1_local.h, &h01, &h01_len);

  interpolate(p1_local.y, p1_local.x, p2_local.y, p2_local.x, &x12, &x12_len);
  interpolate(p1_local.y, p1_local.h, p2_local.y, p2_local.h, &h12, &h12_len);

  interpolate(p0_local.y, p0_local.x, p2_local.y, p2_local.x, &x02, &x02_len);
  interpolate(p0_local.y, p0_local.h, p2_local.y, p2_local.h, &h02, &h02_len);

  concat_sides(&x01, x01_len, &x12, x12_len, &x012);
  concat_sides(&h01, h01_len, &h12, h12_len, &h012);

  m = x02_len / 2;

  if (x02[m] < x012[m]) {
    x_left = x02;
    h_left = h02;
    x_right = x012;
    h_right = h012;
  } else {
    x_left = x012;
    h_left = h012;
    x_right = x02;
    h_right = h02;
  }

  for (y = p0_local.y; y <= p2_local.y; y++) {
    int xl = x_left[y - p0_local.y];
    int xr = x_right[y - p0_local.y];
    float h_segment[SW];
    int h_segment_len;

    interpolate(xl, h_left[y - p0_local.y], xr, h_right[y - p0_local.y],
                &h_segment, &h_segment_len);

    for (x = xl; x <= xr; x++) {
      set_pixel(x, y, shade_color(color, h_segment[x - xl]));
    }
  }
}

struct vector2 viewport_to_canvas(float x, float y) {
  struct vector2 result;

  result.x = x * (SW / VIEWPORT_SIZE);
  result.y = y * (SH / VIEWPORT_SIZE);

  return result;
}

struct vector2 project_vertex(struct vector3 *v) {
  return viewport_to_canvas(v->x * PROJ_PLANE_Z / v->z,
                            v->y * PROJ_PLANE_Z / v->z);
}

int main(void) {
  struct vector2 p0 = {-50, -62, 0.3f};
  struct vector2 p1 = {50, 12, 0.1f};
  struct vector2 p2 = {5, 62, 1.0f};

  struct vector3 vAf = {-0.1f, 0.1f, 2};
  struct vector3 vBf = {0.1f, 0.1f, 2};
  struct vector3 vCf = {0.1f, -0.1f, 2};
  struct vector3 vDf = {-0.1f, -0.1f, 2};

  struct vector3 vAb = {-0.1f, 0.1f, 1};
  struct vector3 vBb = {0.1f, 0.1f, 1};
  struct vector3 vCb = {0.1f, -0.1f, 1};
  struct vector3 vDb = {-0.1f, -0.1f, 1};

  struct vector2 pAf = project_vertex(&vAf);
  struct vector2 pBf = project_vertex(&vBf);
  struct vector2 pCf = project_vertex(&vCf);
  struct vector2 pDf = project_vertex(&vDf);

  struct vector2 pAb = project_vertex(&vAb);
  struct vector2 pBb = project_vertex(&vBb);
  struct vector2 pCb = project_vertex(&vCb);
  struct vector2 pDb = project_vertex(&vDb);

  set_mode(0x13);

  init_palette();

  draw_shaded_triangle(&p0, &p1, &p2, 0);
  draw_wireframe_triangle(&p0, &p1, &p2, shade_color(0, 255));

  draw_line(&pAf, &pBf, shade_color(0, 255));
  draw_line(&pBf, &pCf, shade_color(0, 255));
  draw_line(&pCf, &pDf, shade_color(0, 255));
  draw_line(&pDf, &pAf, shade_color(0, 255));

  draw_line(&pAb, &pBb, shade_color(0, 255));
  draw_line(&pBb, &pCb, shade_color(0, 255));
  draw_line(&pCb, &pDb, shade_color(0, 255));
  draw_line(&pDb, &pAb, shade_color(0, 255));

  draw_line(&pAf, &pAb, shade_color(0, 255));
  draw_line(&pBf, &pBb, shade_color(0, 255));
  draw_line(&pCf, &pCb, shade_color(0, 255));
  draw_line(&pDf, &pDb, shade_color(0, 255));

  getch();

  set_mode(0x03);

  return 0;
}
