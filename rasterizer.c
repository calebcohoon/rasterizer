#include <conio.h>
#include <dos.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define INF 32767 // Infinity for 16 bit compilers
#define SW 320    // Screen width
#define SH 200    // Screen height
#define VIEWPORT_SIZE 1
#define PROJ_PLANE_Z 1.0f
#define PI 3.14159f

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
  unsigned char color;
};

struct model {
  char *name;
  int vert_count;
  int tri_count;
  struct vector3 *vertices;
  struct triangle *triangles;
};

struct instance {
  struct model *model;
  float scale;
  struct vector3 position;
  struct mat4x4 orientation;
  struct mat4x4 transform;
};

struct camera {
  struct vector3 position;
  struct mat4x4 orientation;
};

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
  unsigned char *screen = (unsigned char *)0xA0000;
  int ax, ay;

  ax = SW / 2 + x;
  ay = SH / 2 - y;

  if (ax < 0 || ax >= SW || ay < 0 || ay >= SH) {
    return;
  }

  screen[ay * SW + ax] = color;
}

void clear_screen(unsigned char color) {
  unsigned char *screen = (unsigned char *)0xA0000;

  memset(screen, color, SW * SH);
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

  for (i = i0; i <= i1; i++) {
    arry[idx] = d;
    d = d + a;
    idx++;
  }

  *o_len = idx - 1;
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

void render_triangle(struct triangle *triangle, struct vector2 *proj_verts) {
  draw_wireframe_triangle(&proj_verts[triangle->vertex_index[0]],
                          &proj_verts[triangle->vertex_index[1]],
                          &proj_verts[triangle->vertex_index[2]],
                          shade_color(triangle->color, 31));
}

struct vector2 viewport_to_canvas(float x, float y) {
  struct vector2 result;

  result.x = x * SW / VIEWPORT_SIZE;
  result.y = y * SH / VIEWPORT_SIZE;

  return result;
}

struct vector2 project_vertex(struct vector4 *v) {
  return viewport_to_canvas(v->x * PROJ_PLANE_Z / v->z,
                            v->y * PROJ_PLANE_Z / v->z);
}

void setup_instance_transform(struct instance *i) {
  struct mat4x4 scale = mat4x4_scale(i->scale);
  struct mat4x4 translation = mat4x4_translate_vector3(&i->position);
  struct mat4x4 composed = mat4x4_mul_mat4x4(&i->orientation, &scale);

  i->transform = mat4x4_mul_mat4x4(&translation, &composed);
}

void render_model(struct model *model, struct mat4x4 *matrix) {
  int v, t;
  struct vector2 proj_verts[8];

  if (model->vert_count > 8) {
    return;
  }

  for (v = 0; v < model->vert_count; v++) {
    struct vector3 vertex = model->vertices[v];
    struct vector4 vertex_h = vec4_make(vertex.x, vertex.y, vertex.z, 1);
    struct vector4 vertex_mult = mat4x4_mul_vector4(matrix, &vertex_h);

    proj_verts[v] = project_vertex(&vertex_mult);
  }

  for (t = 0; t < model->tri_count; t++) {
    render_triangle(&model->triangles[t], proj_verts);
  }
}

void render_scene(struct camera *camera, struct instance *instances, int len) {
  int i;
  struct mat4x4 cameraMatrix, m1, m2;
  struct vector3 neg_camera_pos;

  neg_camera_pos = vec3_mul(&camera->position, -1);
  m1 = mat4x4_transpose(&camera->orientation);
  m2 = mat4x4_translate_vector3(&neg_camera_pos);
  cameraMatrix = mat4x4_mul_mat4x4(&m1, &m2);

  for (i = 0; i < len; i++) {
    struct mat4x4 transform =
        mat4x4_mul_mat4x4(&cameraMatrix, &instances[i].transform);

    render_model(instances[i].model, &transform);
  }
}

int main(void) {
  // For the triangle based cube
  int i;
  struct camera camera;
  struct vector3 vertices[8] = {{1, 1, 1},    {-1, 1, 1}, {-1, -1, 1},
                                {1, -1, 1},   {1, 1, -1}, {-1, 1, -1},
                                {-1, -1, -1}, {1, -1, -1}};
  struct triangle triangles[12] = {
      {0, 1, 2, 0}, {0, 2, 3, 0}, {4, 0, 3, 1}, {4, 3, 7, 1},
      {5, 4, 7, 2}, {5, 7, 6, 2}, {1, 5, 6, 3}, {1, 6, 2, 3},
      {4, 5, 1, 4}, {4, 1, 0, 4}, {2, 6, 7, 5}, {2, 7, 3, 5},
  };

  // For drawing the instances of the cube
  struct model the_cube;
  struct instance cube_instances[2];

  // Setup the cube instances
  the_cube.name = "cool cube";
  the_cube.vert_count = 8;
  the_cube.tri_count = 12;
  the_cube.vertices = vertices;
  the_cube.triangles = triangles;

  cube_instances[0].model = &the_cube;
  cube_instances[0].scale = 0.75f;
  cube_instances[0].position.x = -1.5f;
  cube_instances[0].position.y = 0;
  cube_instances[0].position.z = 7;
  cube_instances[0].orientation = mat4x4_identity();
  setup_instance_transform(&cube_instances[0]);

  cube_instances[1].model = &the_cube;
  cube_instances[1].scale = 1;
  cube_instances[1].position.x = 1.25f;
  cube_instances[1].position.y = 2.5f;
  cube_instances[1].position.z = 7.5f;
  cube_instances[1].orientation = mat4x4_rotate_y(195);
  setup_instance_transform(&cube_instances[1]);

  // Position the camera
  camera.position.x = -3;
  camera.position.y = 1;
  camera.position.z = 2;
  camera.orientation = mat4x4_rotate_y(-30);

  set_mode(0x13);

  // Clear screen to white
  clear_screen(shade_color(7, 31));

  init_palette();

  // Render instance of the cube models
  render_scene(&camera, cube_instances, 2);

  getch();

  set_mode(0x03);

  return 0;
}
