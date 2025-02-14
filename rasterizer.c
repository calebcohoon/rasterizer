#include <conio.h>
#include <dos.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define INF 32767 // Infinity for 16 bit compilers
#define SW 320    // Screen width
#define SH 200    // Screen height
#define VIEWPORT_SIZE 1
#define CLIPPING_PLANE_COUNT 5
#define PROJ_PLANE_Z 1.0f
#define PI 3.14159f
#define DEPTH_EPSILON 0.00001f
#define CHAR_WIDTH 8
#define CHAR_HEIGHT 8
#define FIRST_CHAR 65
#define CHAR_COUNT 26
#define AMBIENT_LIGHT 0
#define POINT_LIGHT 1
#define DIRECTIONAL_LIGHT 2

#define SWAP_VECTOR2(a, b)                                                     \
  do {                                                                         \
    struct vector2 temp = (a);                                                 \
    (a) = (b);                                                                 \
    (b) = temp;                                                                \
  } while (0)

struct mat4x4 {
  float a, b, c, d;
  float e, f, g, h;
  float i, j, k, l;
  float m, n, o, p;
};

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

struct vector4 {
  float x;
  float y;
  float z;
  float w;
};

struct triangle {
  int vertex_index[3];
};

struct model {
  char *name;
  float radius;
  struct vector3 center;
  int vert_count;
  int tri_count;
  struct vector3 *vertices;
  struct triangle *triangles;
};

struct instance {
  struct model *model;
  unsigned char color;
  float scale;
  struct vector3 position;
  struct mat4x4 orientation;
  struct mat4x4 transform;
};

struct plane {
  struct vector3 normal;
  float distance;
};

struct camera {
  struct vector3 position;
  struct mat4x4 orientation;
  struct plane clipping_planes[5];
};

struct light {
  int type;
  float intensity;
  struct vector3 vector;
};

unsigned char back_buffer[SW * SH];
float depth_buffer[SW * SH];
const unsigned char font_data[CHAR_COUNT][8] = {
    /* A */ {0x00, 0x66, 0x66, 0x7E, 0x66, 0x66, 0x3C, 0x18},
    /* B */ {0x00, 0x7C, 0x66, 0x66, 0x7C, 0x66, 0x66, 0x7C},
    /* C */ {0x00, 0x3C, 0x66, 0x60, 0x60, 0x60, 0x66, 0x3C},
    /* D */ {0x00, 0x78, 0x6C, 0x66, 0x66, 0x66, 0x6C, 0x78},
    /* E */ {0x00, 0x7E, 0x60, 0x60, 0x78, 0x60, 0x60, 0x7E},
    /* F */ {0x00, 0x60, 0x60, 0x60, 0x78, 0x60, 0x60, 0x7E},
    /* G */ {0x00, 0x3C, 0x66, 0x66, 0x6E, 0x60, 0x66, 0x3C},
    /* H */ {0x00, 0x66, 0x66, 0x66, 0x7E, 0x66, 0x66, 0x66},
    /* I */ {0x00, 0x3C, 0x18, 0x18, 0x18, 0x18, 0x18, 0x3C},
    /* J */ {0x00, 0x38, 0x6C, 0x6C, 0x0C, 0x0C, 0x0C, 0x1E},
    /* K */ {0x00, 0x66, 0x6C, 0x78, 0x70, 0x78, 0x6C, 0x66},
    /* L */ {0x00, 0x7E, 0x60, 0x60, 0x60, 0x60, 0x60, 0x60},
    /* M */ {0x00, 0x63, 0x63, 0x63, 0x6B, 0x7F, 0x77, 0x63},
    /* N */ {0x00, 0x66, 0x66, 0x6E, 0x7E, 0x7E, 0x76, 0x66},
    /* O */ {0x00, 0x3C, 0x66, 0x66, 0x66, 0x66, 0x66, 0x3C},
    /* P */ {0x00, 0x60, 0x60, 0x60, 0x7C, 0x66, 0x66, 0x7C},
    /* Q */ {0x00, 0x36, 0x6C, 0x6A, 0x66, 0x66, 0x66, 0x3C},
    /* R */ {0x00, 0x66, 0x6C, 0x78, 0x7C, 0x66, 0x66, 0x7C},
    /* S */ {0x00, 0x3C, 0x66, 0x06, 0x3C, 0x60, 0x66, 0x3C},
    /* T */ {0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7E},
    /* U */ {0x00, 0x3C, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66},
    /* V */ {0x00, 0x18, 0x3C, 0x66, 0x66, 0x66, 0x66, 0x66},
    /* W */ {0x00, 0x63, 0x77, 0x7F, 0x6B, 0x63, 0x63, 0x63},
    /* X */ {0x00, 0x66, 0x66, 0x3C, 0x18, 0x3C, 0x66, 0x66},
    /* Y */ {0x00, 0x18, 0x18, 0x18, 0x3C, 0x66, 0x66, 0x66},
    /* Z */ {0x00, 0x7E, 0x60, 0x30, 0x18, 0x0C, 0x06, 0x7E}};

// Build an optimized palette for the main colors being used
void init_palette() {
  int color, shade;
  unsigned char index;
  float intensity;
  unsigned char base_colors[8][3] = {{63, 0, 0}, // Red
                                     {0, 63, 0}, // Green
                                     {0, 0, 63}, // Blue
                                     {63, 63, 0},  {63, 0, 63}, {0, 63, 63},
                                     {32, 63, 63}, {32, 63, 63}};

  for (color = 0; color < 8; color++) {
    for (shade = 0; shade < 32; shade++) {
      index = color * 32 + shade;
      intensity = shade / 31.0f;

      outp(0x3C8, index);

      // Add a single shade of white in the the last color's shade
      if (color == 7 && shade == 31) {
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

void init_depth_buffer(void) {
  memset(depth_buffer, 0, sizeof(float) * SW * SH);
}

unsigned char shade_color(unsigned char color, float intensity) {
  unsigned char shade;
  if (intensity > 1.0f)
    intensity = 1.0f;
  if (intensity < 0.0f)
    intensity = 0.0f;

  shade = (unsigned char)(intensity * 31.0f);
  return color * 32 + shade;
}

void set_mode(unsigned char mode) {
  union REGS regs;

  regs.w.ax = mode;

  int386(0x10, &regs, &regs);
}

void set_pixel(int x, int y, char color) {
  int ax, ay;

  ax = SW / 2 + x;
  ay = SH / 2 - y;

  if (ax < 0 || ax >= SW || ay < 0 || ay >= SH) {
    return;
  }

  back_buffer[ay * SW + ax] = color;
}

void draw_char(int x, int y, char c, unsigned char color) {
  int row, col;
  int char_index;
  unsigned char row_data;

  if (c < 'A' || c > 'Z') {
    return;
  }

  char_index = c - FIRST_CHAR;

  for (row = 0; row < CHAR_HEIGHT; row++) {
    row_data = font_data[char_index][row];

    for (col = 0; col < CHAR_WIDTH; col++) {
      if (row_data & (0x80 >> col)) {
        set_pixel(x + col, y + row, color);
      }
    }
  }
}

void draw_text(int x, int y, const char *text, unsigned char color) {
  int pos_x = x;

  while (*text) {
    draw_char(pos_x, y, *text, color);
    pos_x += CHAR_WIDTH + 1;
    text++;
  }
}

int is_closer_pixel(int x, int y, float inv_z) {
  int ax, ay, offset;

  ax = SW / 2 + x;
  ay = SH / 2 - y;

  if (ax < 0 || ax >= SW || ay < 0 || ay >= SH) {
    return 0;
  }

  if (depth_buffer[ay * SW + ax] < (inv_z - DEPTH_EPSILON)) {
    depth_buffer[ay * SW + ax] = inv_z;
    return 1;
  }

  return 0;
}

void clear_screen(unsigned char color) { memset(back_buffer, color, SW * SH); }

void present_buffer(void) {
  unsigned char *screen = (unsigned char *)0xA0000;

  memcpy(screen, back_buffer, SW * SH);
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

struct vector3 vec3_cross(struct vector3 *a, struct vector3 *b) {
  struct vector3 result;

  result.x = (a->y * b->z) - (a->z * b->y);
  result.y = (a->z * b->x) - (a->x * b->z);
  result.z = (a->x * b->y) - (a->y * b->x);

  return result;
}

struct vector3 vec3_make(float x, float y, float z) {
  struct vector3 result;

  result.x = x;
  result.y = y;
  result.z = z;

  return result;
}

struct vector4 vec4_make(float x, float y, float z, float w) {
  struct vector4 result;

  result.x = x;
  result.y = y;
  result.z = z;
  result.w = w;

  return result;
}

struct vector4 mat4x4_mul_vector4(const struct mat4x4 *matrix,
                                  const struct vector4 *vector) {
  struct vector4 result;

  result.x = matrix->a * vector->x + matrix->b * vector->y +
             matrix->c * vector->z + matrix->d * vector->w;

  result.y = matrix->e * vector->x + matrix->f * vector->y +
             matrix->g * vector->z + matrix->h * vector->w;

  result.z = matrix->i * vector->x + matrix->j * vector->y +
             matrix->k * vector->z + matrix->l * vector->w;

  result.w = matrix->m * vector->x + matrix->n * vector->y +
             matrix->o * vector->z + matrix->p * vector->w;

  return result;
}

struct mat4x4 mat4x4_identity(void) {
  struct mat4x4 result;

  result.a = 1.0f;
  result.f = 1.0f;
  result.k = 1.0f;
  result.p = 1.0f;

  result.b = 0.0f;
  result.c = 0.0f;
  result.d = 0.0f;
  result.e = 0.0f;
  result.g = 0.0f;
  result.h = 0.0f;
  result.i = 0.0f;
  result.j = 0.0f;
  result.l = 0.0f;
  result.m = 0.0f;
  result.n = 0.0f;
  result.o = 0.0f;

  return result;
}

struct mat4x4 mat4x4_transpose(struct mat4x4 *m) {
  struct mat4x4 result;

  result.a = m->a;
  result.b = m->e;
  result.c = m->i;
  result.d = m->m;
  result.e = m->b;
  result.f = m->f;
  result.g = m->j;
  result.h = m->n;
  result.i = m->c;
  result.j = m->g;
  result.k = m->k;
  result.l = m->o;
  result.m = m->d;
  result.n = m->h;
  result.o = m->l;
  result.p = m->p;

  return result;
}

struct mat4x4 mat4x4_mul_mat4x4(struct mat4x4 *m1, struct mat4x4 *m2) {
  struct mat4x4 result;

  result.a = m1->a * m2->a + m1->b * m2->e + m1->c * m2->i + m1->d * m2->m;
  result.b = m1->a * m2->b + m1->b * m2->f + m1->c * m2->j + m1->d * m2->n;
  result.c = m1->a * m2->c + m1->b * m2->g + m1->c * m2->k + m1->d * m2->o;
  result.d = m1->a * m2->d + m1->b * m2->h + m1->c * m2->l + m1->d * m2->p;

  result.e = m1->e * m2->a + m1->f * m2->e + m1->g * m2->i + m1->h * m2->m;
  result.f = m1->e * m2->b + m1->f * m2->f + m1->g * m2->j + m1->h * m2->n;
  result.g = m1->e * m2->c + m1->f * m2->g + m1->g * m2->k + m1->h * m2->o;
  result.h = m1->e * m2->d + m1->f * m2->h + m1->g * m2->l + m1->h * m2->p;

  result.i = m1->i * m2->a + m1->j * m2->e + m1->k * m2->i + m1->l * m2->m;
  result.j = m1->i * m2->b + m1->j * m2->f + m1->k * m2->j + m1->l * m2->n;
  result.k = m1->i * m2->c + m1->j * m2->g + m1->k * m2->k + m1->l * m2->o;
  result.l = m1->i * m2->d + m1->j * m2->h + m1->k * m2->l + m1->l * m2->p;

  result.m = m1->m * m2->a + m1->n * m2->e + m1->o * m2->i + m1->p * m2->m;
  result.n = m1->m * m2->b + m1->n * m2->f + m1->o * m2->j + m1->p * m2->n;
  result.o = m1->m * m2->c + m1->n * m2->g + m1->o * m2->k + m1->p * m2->o;
  result.p = m1->m * m2->d + m1->n * m2->h + m1->o * m2->l + m1->p * m2->p;

  return result;
}

struct mat4x4 mat4x4_rotate_y(float degrees) {
  float cos_r = cos(degrees * PI / 180.0f);
  float sin_r = sin(degrees * PI / 180.0f);
  struct mat4x4 result;

  result.a = cos_r;
  result.b = 0;
  result.c = -sin_r;
  result.d = 0;

  result.e = 0;
  result.f = 1;
  result.g = 0;
  result.h = 0;

  result.i = sin_r;
  result.j = 0;
  result.k = cos_r;
  result.l = 0;

  result.m = 0;
  result.n = 0;
  result.o = 0;
  result.p = 1;

  return result;
}

struct mat4x4 mat4x4_translate_vector3(struct vector3 *v) {
  struct mat4x4 result;

  result.a = 1;
  result.b = 0;
  result.c = 0;
  result.d = v->x;

  result.e = 0;
  result.f = 1;
  result.g = 0;
  result.h = v->y;

  result.i = 0;
  result.j = 0;
  result.k = 1;
  result.l = v->z;

  result.m = 0;
  result.n = 0;
  result.o = 0;
  result.p = 1;

  return result;
}

struct mat4x4 mat4x4_scale(float scale) {
  struct mat4x4 result;

  result.a = scale;
  result.b = 0;
  result.c = 0;
  result.d = 0;

  result.e = 0;
  result.f = scale;
  result.g = 0;
  result.h = 0;

  result.i = 0;
  result.j = 0;
  result.k = scale;
  result.l = 0;

  result.m = 0;
  result.n = 0;
  result.o = 0;
  result.p = 1;

  return result;
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

  for (i = i0; i <= i1 && idx < SW; i++) {
    arry[idx] = d;
    d = d + a;
    idx++;
  }

  *o_len = idx;
}

void concat_sides(float *x01, int x01_len, float *x12, int x12_len,
                  float *x012) {
  int i;

  // Ignore the last value in the x01 array just like in the book
  for (i = 0; i <= x01_len - 1; i++) {
    x012[i] = x01[i];
  }

  for (i = 0; i <= x12_len; i++) {
    x012[i + (x01_len - 1)] = x12[i];
  }
}

struct vector3 compute_triangle_normal(struct vector3 *v0, struct vector3 *v1,
                                       struct vector3 *v2) {
  struct vector3 v0v1 = vec3_sub(v0, v1);
  struct vector3 v0v2 = vec3_sub(v0, v2);

  return vec3_cross(&v0v1, &v0v2);
}

float compute_illumination(struct vector3 *vertex, struct vector3 *normal,
                           struct camera *camera, struct light *lights,
                           int light_len) {
  float illumination = 0.0f, cos_alpha;
  struct vector3 vl;
  int l;

  for (l = 0; l < light_len; l++) {
    struct light *curr_light = &lights[l];

    if (curr_light->type == AMBIENT_LIGHT) {
      illumination += curr_light->intensity;
      continue;
    }

    if (curr_light->type == DIRECTIONAL_LIGHT) {
      struct vector4 light_vector = vec4_make(
          curr_light->vector.x, curr_light->vector.y, curr_light->vector.z, 1);
      struct mat4x4 camera_matrix = mat4x4_transpose(&camera->orientation);
      struct vector4 rotated_light =
          mat4x4_mul_vector4(&camera_matrix, &light_vector);

      vl = vec3_make(rotated_light.x, rotated_light.y, rotated_light.z);
    } else if (curr_light->type == POINT_LIGHT) {
      struct vector4 light_vector = vec4_make(
          curr_light->vector.x, curr_light->vector.y, curr_light->vector.z, 1);
      struct vector3 camera_pos_inv = vec3_mul(&camera->position, -1);
      struct mat4x4 camera_translation =
          mat4x4_translate_vector3(&camera_pos_inv);
      struct mat4x4 camera_transposed = mat4x4_transpose(&camera->orientation);
      struct mat4x4 camera_matrix =
          mat4x4_mul_mat4x4(&camera_transposed, &camera_translation);
      struct vector4 transformed_light =
          mat4x4_mul_vector4(&camera_matrix, &light_vector);
      struct vector3 transformed_light_v3 = vec3_make(
          transformed_light.x, transformed_light.y, transformed_light.z);
      struct vector3 vertex_inv = vec3_mul(vertex, -1);
      vl = vec3_add(&vertex_inv, &transformed_light_v3);
    }

    cos_alpha = vec3_dot(normal, &vl) / (vec3_len(&vl) * vec3_len(normal));
    if (cos_alpha > 0) {
      illumination += cos_alpha * curr_light->intensity;
    }
  }

  return illumination;
}

void draw_filled_triangle(struct vector2 *p0, struct vector2 *p1,
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
  float intensity = 0;

  // Add bounds checking for y-coordinates
  if (abs(p0_local.y) >= SH || abs(p1_local.y) >= SH || abs(p2_local.y) >= SH) {
    return;
  }

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
      if (is_closer_pixel(x, y, h_segment[x - xl])) {
        set_pixel(x, y, color);
      }
    }
  }
}

void render_triangle(struct triangle *triangle,
                     struct vector3 *transformed_verts,
                     struct vector2 *proj_verts, unsigned char color,
                     struct camera *camera, struct light *lights,
                     int light_len) {
  struct vector3 v0 = transformed_verts[triangle->vertex_index[0]];
  struct vector3 v1 = transformed_verts[triangle->vertex_index[1]];
  struct vector3 v2 = transformed_verts[triangle->vertex_index[2]];
  struct vector3 normal = compute_triangle_normal(&v0, &v1, &v2);
  struct vector3 vert_to_camera = vec3_mul(&v0, -1);
  struct vector3 center, normal_inv;
  float intensity;

  // Cull back face geometry
  if (vec3_dot(&normal, &vert_to_camera) <= 0) {
    return;
  }

  // Compute lighting intensity
  normal_inv = compute_triangle_normal(&v0, &v1, &v2);
  center = vec3_make((v0.x + v1.x + v2.x) / 3.0f, (v0.y + v1.y + v2.y) / 3.0f,
                     (v0.z + v1.z + v2.z) / 3.0f);
  intensity =
      compute_illumination(&center, &normal_inv, camera, lights, light_len);

  draw_filled_triangle(&proj_verts[triangle->vertex_index[0]],
                       &proj_verts[triangle->vertex_index[1]],
                       &proj_verts[triangle->vertex_index[2]],
                       shade_color(color, intensity));
}

struct vector2 viewport_to_canvas(float x, float y) {
  struct vector2 result;
  float aspect_ratio = 1.2f; // Account for having non-square pixels

  result.x = x * SW / VIEWPORT_SIZE;
  result.y = (y * aspect_ratio) * SH / VIEWPORT_SIZE;

  return result;
}

struct vector2 project_vertex(struct vector4 *v) {
  // Prevent a divide by zero
  if (v->z < 0.1f || v->z == 0) {
    v->z = 0.1f;
  }

  return viewport_to_canvas(v->x * PROJ_PLANE_Z / v->z,
                            v->y * PROJ_PLANE_Z / v->z);
}

void setup_instance_transform(struct instance *i) {
  struct mat4x4 scale = mat4x4_scale(i->scale);
  struct mat4x4 translation = mat4x4_translate_vector3(&i->position);
  struct mat4x4 composed = mat4x4_mul_mat4x4(&i->orientation, &scale);

  i->transform = mat4x4_mul_mat4x4(&translation, &composed);
}

int is_clipped(struct plane *planes, struct model *model, float scale,
               struct mat4x4 *transform) {
  int p;
  float radius = model->radius * scale;
  struct vector4 center =
      vec4_make(model->center.x, model->center.y, model->center.z, 1);
  center = mat4x4_mul_vector4(transform, &center);

  for (p = 0; p < CLIPPING_PLANE_COUNT; p++) {
    struct vector3 center_v3 = vec3_make(center.x, center.y, center.z);
    float distance =
        vec3_dot(&planes[p].normal, &center_v3) + planes[p].distance;
    if (distance < -radius) {
      return 1;
    }
  }

  return 0;
}

void render_model(struct model *model, unsigned char color,
                  struct mat4x4 *matrix, struct camera *camera,
                  struct light *lights, int light_len) {
  int v, t;
  struct vector3 *transformed_verts;
  struct vector2 *projected_verts;

  /* Allocate memory for transformed and projected vertices */
  transformed_verts =
      (struct vector3 *)malloc(model->vert_count * sizeof(struct vector3));
  projected_verts =
      (struct vector2 *)malloc(model->vert_count * sizeof(struct vector2));

  if (!transformed_verts || !projected_verts) {
    /* Handle allocation failure */
    if (transformed_verts)
      free(transformed_verts);
    if (projected_verts)
      free(projected_verts);
    return;
  }

  for (v = 0; v < model->vert_count; v++) {
    struct vector3 vertex = model->vertices[v];
    struct vector4 vertex_h = vec4_make(vertex.x, vertex.y, vertex.z, 1);
    struct vector4 vertex_mult = mat4x4_mul_vector4(matrix, &vertex_h);

    transformed_verts[v] =
        vec3_make(vertex_mult.x, vertex_mult.y, vertex_mult.z);
    projected_verts[v] = project_vertex(&vertex_mult);

    /* Repurpose H with the transformed Z later used with z-buffer */
    projected_verts[v].h = 1.0f / vertex_mult.z;
  }

  for (t = 0; t < model->tri_count; t++) {
    render_triangle(&model->triangles[t], transformed_verts, projected_verts,
                    color, camera, lights, light_len);
  }

  /* Free allocated memory */
  free(transformed_verts);
  free(projected_verts);
}

void render_scene(struct camera *camera, struct instance *instances,
                  int instance_len, struct light *lights, int light_len) {
  int i;
  struct mat4x4 cameraMatrix, m1, m2;
  struct vector3 neg_camera_pos;

  neg_camera_pos = vec3_mul(&camera->position, -1);
  m1 = mat4x4_transpose(&camera->orientation);
  m2 = mat4x4_translate_vector3(&neg_camera_pos);
  cameraMatrix = mat4x4_mul_mat4x4(&m1, &m2);

  for (i = 0; i < instance_len; i++) {
    struct mat4x4 transform =
        mat4x4_mul_mat4x4(&cameraMatrix, &instances[i].transform);

    int status = is_clipped(camera->clipping_planes, instances[i].model,
                            instances[i].scale, &transform);

    if (status == 0) {
      render_model(instances[i].model, instances[i].color, &transform, camera,
                   lights, light_len);
    }
  }
}

struct model *generate_sphere(int divs) {
  struct model *sphere;
  struct vector3 *vertices;
  struct triangle *triangles;
  float delta_angle;
  int d, i;
  int vertex_count, triangle_count;
  int i0, i1, i2;
  float y, radius;

  vertex_count = divs * (divs + 1);
  triangle_count = divs * divs * 2;

  sphere = (struct model *)malloc(sizeof(struct model));
  vertices = (struct vector3 *)malloc(vertex_count * sizeof(struct vector3));
  triangles =
      (struct triangle *)malloc(triangle_count * sizeof(struct triangle));

  if (!sphere || !vertices || !triangles) {
    if (vertices)
      free(vertices);
    if (triangles)
      free(triangles);
    if (sphere)
      free(sphere);
    return NULL;
  }

  delta_angle = 2.0f * PI / (float)divs;

  for (d = 0; d < divs + 1; d++) {
    y = (2.0f / (float)divs) * ((float)d - (float)divs / 2.0f);
    radius = sqrt(1.0f - y * y);

    for (i = 0; i < divs; i++) {
      int vertex_index = d * divs + i;
      vertices[vertex_index].x = radius * cos(i * delta_angle);
      vertices[vertex_index].y = y;
      vertices[vertex_index].z = radius * sin(i * delta_angle);
    }
  }

  for (d = 0; d < divs; d++) {
    for (i = 0; i < divs; i++) {
      int triangle_index = (d * divs + i) * 2;

      i0 = d * divs + i;
      i1 = (d + 1) * divs + ((i + 1) % divs);
      i2 = divs * d + ((i + 1) % divs);

      triangles[triangle_index].vertex_index[0] = i0;
      triangles[triangle_index].vertex_index[1] = i1;
      triangles[triangle_index].vertex_index[2] = i2;

      triangles[triangle_index + 1].vertex_index[0] = i0;
      triangles[triangle_index + 1].vertex_index[1] = i0 + divs;
      triangles[triangle_index + 1].vertex_index[2] = i1;
    }
  }

  sphere->name = "sphere";
  sphere->vertices = vertices;
  sphere->triangles = triangles;
  sphere->vert_count = vertex_count;
  sphere->tri_count = triangle_count;
  sphere->center = vec3_make(0.0f, 0.0f, 0.0f);
  sphere->radius = 1.0f;

  return sphere;
}

struct model *generate_plane(int divs, float size) {
  struct model *plane;
  struct vector3 *vertices;
  struct triangle *triangles;
  int vertex_count, triangle_count;
  int x, z, i;
  float half_size, step;

  vertex_count = divs * divs;
  triangle_count = (divs - 1) * (divs - 1) * 2;

  plane = (struct model *)malloc(sizeof(struct model));
  vertices = (struct vector3 *)malloc(vertex_count * sizeof(struct vector3));
  triangles =
      (struct triangle *)malloc(triangle_count * sizeof(struct triangle));

  if (!plane || !vertices || !triangles) {
    if (vertices)
      free(vertices);
    if (triangles)
      free(triangles);
    if (plane)
      free(plane);
    return NULL;
  }

  half_size = size * 0.5f;
  step = size / (float)(divs - 1);

  for (z = 0; z < divs; z++) {
    for (x = 0; x < divs; x++) {
      int vertex_index = z * divs + x;
      vertices[vertex_index].x = -half_size + x * step;
      vertices[vertex_index].y = 0.0f;
      vertices[vertex_index].z = -half_size + z * step;
    }
  }

  i = 0;
  for (z = 0; z < divs - 1; z++) {
    for (x = 0; x < divs - 1; x++) {
      int top_left = z * divs + x;
      int top_right = top_left + 1;
      int bottom_left = (z + 1) * divs + x;
      int bottom_right = bottom_left + 1;

      triangles[i].vertex_index[0] = top_left;
      triangles[i].vertex_index[1] = bottom_left;
      triangles[i].vertex_index[2] = top_right;
      i++;

      triangles[i].vertex_index[0] = bottom_left;
      triangles[i].vertex_index[1] = bottom_right;
      triangles[i].vertex_index[2] = top_right;
      i++;
    }
  }

  plane->name = "plane";
  plane->vertices = vertices;
  plane->triangles = triangles;
  plane->vert_count = vertex_count;
  plane->tri_count = triangle_count;
  plane->center = vec3_make(0.0f, 0.0f, 0.0f);
  plane->radius = (float)sqrt(2 * (size * size)) * 0.5f;

  return plane;
}

void free_model(struct model *model) {
  if (model) {
    if (model->vertices)
      free(model->vertices);
    if (model->triangles)
      free(model->triangles);
    free(model);
  }
}

int main(void) {
  float sqrt_2 = 1.0f / sqrt(2);
  struct camera camera;
  struct model *sphere = generate_sphere(13);
  struct model *plane = generate_plane(20, 30.0f);
  struct instance instances[4];
  struct light lights[3];

  lights[0].type = AMBIENT_LIGHT;
  lights[0].intensity = 0.3f;

  lights[1].type = DIRECTIONAL_LIGHT;
  lights[1].intensity = 0.2f;
  lights[1].vector = vec3_make(-1, 0, 4);

  lights[2].type = POINT_LIGHT;
  lights[2].intensity = 0.6f;
  lights[2].vector = vec3_make(-2, 2, -2);

  instances[0].model = sphere;
  instances[0].scale = 1;
  instances[0].position.x = 0;
  instances[0].position.y = -1;
  instances[0].position.z = 3;
  instances[0].orientation = mat4x4_identity();
  instances[0].color = 0;
  setup_instance_transform(&instances[0]);

  instances[1].model = sphere;
  instances[1].scale = 1;
  instances[1].position.x = -2;
  instances[1].position.y = 0;
  instances[1].position.z = 4;
  instances[1].orientation = mat4x4_identity();
  instances[1].color = 1;
  setup_instance_transform(&instances[1]);

  instances[2].model = sphere;
  instances[2].scale = 1;
  instances[2].position.x = 2;
  instances[2].position.y = 0;
  instances[2].position.z = 4;
  instances[2].orientation = mat4x4_identity();
  instances[2].color = 2;
  setup_instance_transform(&instances[2]);

  instances[3].model = plane;
  instances[3].scale = 1.0f;
  instances[3].position = vec3_make(0.0f, -1.0f, 10.0f);
  instances[3].orientation = mat4x4_identity();
  instances[3].color = 3;
  setup_instance_transform(&instances[3]);

  // Position the camera
  camera.position.x = 0;
  camera.position.y = 0;
  camera.position.z = -0.5f;
  camera.orientation = mat4x4_identity();

  // Near plane
  camera.clipping_planes[0].normal = vec3_make(0, 0, 1);
  camera.clipping_planes[0].distance = -1;

  // Left plane
  camera.clipping_planes[1].normal = vec3_make(sqrt_2, 0, sqrt_2);
  camera.clipping_planes[1].distance = 0;

  // Right plane
  camera.clipping_planes[2].normal = vec3_make(-sqrt_2, 0, sqrt_2);
  camera.clipping_planes[2].distance = 0;

  // Top plane
  camera.clipping_planes[3].normal = vec3_make(0, -sqrt_2, sqrt_2);
  camera.clipping_planes[3].distance = 0;

  // Bottom plane
  camera.clipping_planes[4].normal = vec3_make(0, sqrt_2, sqrt_2);
  camera.clipping_planes[4].distance = 0;

  set_mode(0x13);

  init_palette();

  init_depth_buffer();

  clear_screen(shade_color(7, 31) /* White */);

  render_scene(&camera, instances, 4, lights, 3);

  present_buffer();

  getch();

  set_mode(0x03);

  free_model(sphere);
  free_model(plane);

  return 0;
}
