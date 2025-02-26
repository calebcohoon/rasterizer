#include <conio.h>
#include <dos.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Constants & Definitions */
#define INF 32767              // Infinity for 16 bit compilers
#define SW 320                 // Screen width
#define SH 200                 // Screen height
#define VIEWPORT_SIZE 1        // Size of the viewport in world units
#define PROJ_PLANE_Z 1.0f      // Distance to projection plane
#define PI 3.14159f            // Pi constant
#define DEPTH_EPSILON 0.00001f // Small value for depth comparisons

/* Text rendering */
#define CHAR_WIDTH 8
#define CHAR_HEIGHT 8
#define FIRST_CHAR 65
#define CHAR_COUNT 26

/* Light types */
#define AMBIENT_LIGHT 0
#define POINT_LIGHT 1
#define DIRECTIONAL_LIGHT 2

/* Frustum clipping */
#define CLIPPING_PLANE_COUNT 5

/* Utility macros */
#define SWAP_VECTOR2(a, b)                                                     \
  do {                                                                         \
    vector2_t temp = (a);                                                      \
    (a) = (b);                                                                 \
    (b) = temp;                                                                \
  } while (0)

/* Core Data Structures */
typedef struct {
  float a, b, c, d;
  float e, f, g, h;
  float i, j, k, l;
  float m, n, o, p;
} mat4x4_t;

typedef struct {
  int x;
  int y;
  float h; // Used for Z-buffer calculations
} vector2_t;

typedef struct {
  float x;
  float y;
  float z;
} vector3_t;

typedef struct {
  float x;
  float y;
  float z;
  float w;
} vector4_t;

/* 3D Geometry */
typedef struct {
  int vertex_index[3];
} triangle_t;

typedef struct {
  char *name;
  float radius; // Bounding sphere radius for culling
  vector3_t center;
  int vert_count;
  int tri_count;
  vector3_t *vertices;
  triangle_t *triangles;
} model_t;

typedef struct {
  model_t *model;
  unsigned char color;
  float scale;
  vector3_t position;
  mat4x4_t orientation;
  mat4x4_t transform;
} instance_t;

/* Scene Elements */
typedef struct {
  vector3_t normal;
  float distance;
} plane_t;

typedef struct {
  vector3_t position;
  mat4x4_t orientation;
  plane_t clipping_planes[CLIPPING_PLANE_COUNT];
} camera_t;

typedef struct {
  int type;
  float intensity;
  vector3_t vector; // Position or direction depending on type
} light_t;

/* Global Variables */
unsigned char back_buffer[SW * SH];
float depth_buffer[SW * SH];

/* Graphics System Function Prototypes */
void init_palette(void);
void init_depth_buffer(void);
unsigned char shade_color(unsigned char color, float intensity);
void set_mode(unsigned char mode);
void set_pixel(int x, int y, char color);
void clear_screen(unsigned char color);
void present_buffer(void);
int is_closer_pixel(int x, int y, float inv_z);

/* Vector & Matrix Operations Function Prototypes */
vector3_t vec3_add(vector3_t *a, vector3_t *b);
vector3_t vec3_sub(vector3_t *a, vector3_t *b);
vector3_t vec3_mul(vector3_t *v, float s);
vector3_t vec3_div(vector3_t *v, float d);
float vec3_len(vector3_t *v);
float vec3_dot(vector3_t *a, vector3_t *b);
vector3_t vec3_cross(vector3_t *a, vector3_t *b);
vector3_t vec3_make(float x, float y, float z);
vector4_t vec4_make(float x, float y, float z, float w);
mat4x4_t mat4x4_identity(void);
mat4x4_t mat4x4_transpose(mat4x4_t *m);
mat4x4_t mat4x4_mul_mat4x4(mat4x4_t *m1, mat4x4_t *m2);
vector4_t mat4x4_mul_vector4(const mat4x4_t *matrix, const vector4_t *vector);
mat4x4_t mat4x4_rotate_y(float degrees);
mat4x4_t mat4x4_translate_vector3(vector3_t *v);
mat4x4_t mat4x4_scale(float scale);

/* Triangle Rasterization Function Prototypes */
void interpolate(float i0, float d0, float i1, float d1, float *arry,
                 int *o_len);
void concat_sides(float *x01, int x01_len, float *x12, int x12_len,
                  float *x012);
void draw_filled_triangle(vector2_t *p0, vector2_t *p1, vector2_t *p2,
                          unsigned char color);

/* Lighting & Shading Function Prototypes */
vector3_t compute_triangle_normal(vector3_t *v0, vector3_t *v1, vector3_t *v2);
float compute_illumination(vector3_t *vertex, vector3_t *normal,
                           camera_t *camera, light_t *lights, int light_len);

/* 3D Rendering Pipeline Function Prototypes */
vector2_t viewport_to_canvas(float x, float y);
vector2_t project_vertex(vector4_t *v);
void setup_instance_transform(instance_t *i);
int is_clipped(plane_t *planes, model_t *model, float scale,
               mat4x4_t *transform);
void render_triangle(triangle_t *triangle, vector3_t *transformed_verts,
                     vector2_t *proj_verts, unsigned char color,
                     camera_t *camera, light_t *lights, int light_len);
void render_model(model_t *model, unsigned char color, mat4x4_t *matrix,
                  camera_t *camera, light_t *lights, int light_len);
void render_scene(camera_t *camera, instance_t *instances, int instance_len,
                  light_t *lights, int light_len);

/* Model Generation Function Prototypes */
model_t *generate_sphere(int divs);
model_t *generate_plane(int divs, float size);
void free_model(model_t *model);

/* Implementations */
void init_palette(void) {
  int color, shade;
  unsigned char index;
  float intensity;
  unsigned char base_colors[8][3] = {
      {63, 0, 0},   // Red
      {0, 63, 0},   // Green
      {0, 0, 63},   // Blue
      {63, 63, 0},  // Yellow
      {63, 0, 63},  // Magenta
      {0, 63, 63},  // Cyan
      {32, 63, 63}, // Teal
      {32, 63, 63}  // Teal (duplicate, last reserved for white)
  };

  // Build an optimized palette for the main colors being used
  for (color = 0; color < 8; color++) {
    for (shade = 0; shade < 32; shade++) {
      index = color * 32 + shade;
      intensity = shade / 31.0f;

      outp(0x3C8, index);

      // Add a single shade of white in the last color's shade
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

vector3_t vec3_add(vector3_t *a, vector3_t *b) {
  vector3_t result;

  result.x = a->x + b->x;
  result.y = a->y + b->y;
  result.z = a->z + b->z;

  return result;
}

vector3_t vec3_sub(vector3_t *a, vector3_t *b) {
  vector3_t result;

  result.x = b->x - a->x;
  result.y = b->y - a->y;
  result.z = b->z - a->z;

  return result;
}

vector3_t vec3_mul(vector3_t *v, float s) {
  vector3_t result;

  result.x = v->x * s;
  result.y = v->y * s;
  result.z = v->z * s;

  return result;
}

vector3_t vec3_div(vector3_t *v, float d) {
  vector3_t result;
  float recip = 1.0f / d; // Calculate reciprocal once for efficiency

  result.x = v->x * recip;
  result.y = v->y * recip;
  result.z = v->z * recip;

  return result;
}

float vec3_len(vector3_t *v) {
  return sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
}

float vec3_dot(vector3_t *a, vector3_t *b) {
  return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}

vector3_t vec3_cross(vector3_t *a, vector3_t *b) {
  vector3_t result;

  result.x = (a->y * b->z) - (a->z * b->y);
  result.y = (a->z * b->x) - (a->x * b->z);
  result.z = (a->x * b->y) - (a->y * b->x);

  return result;
}

vector3_t vec3_make(float x, float y, float z) {
  vector3_t result;

  result.x = x;
  result.y = y;
  result.z = z;

  return result;
}

vector4_t vec4_make(float x, float y, float z, float w) {
  vector4_t result;

  result.x = x;
  result.y = y;
  result.z = z;
  result.w = w;

  return result;
}

vector4_t mat4x4_mul_vector4(const mat4x4_t *matrix, const vector4_t *vector) {
  vector4_t result;

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

mat4x4_t mat4x4_identity(void) {
  mat4x4_t result;

  // Initialize to zero
  memset(&result, 0, sizeof(mat4x4_t));

  // Set diagonal to 1
  result.a = 1.0f;
  result.f = 1.0f;
  result.k = 1.0f;
  result.p = 1.0f;

  return result;
}

mat4x4_t mat4x4_transpose(mat4x4_t *m) {
  mat4x4_t result;

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

mat4x4_t mat4x4_mul_mat4x4(mat4x4_t *m1, mat4x4_t *m2) {
  mat4x4_t result;

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

mat4x4_t mat4x4_rotate_y(float degrees) {
  float cos_r = cos(degrees * PI / 180.0f);
  float sin_r = sin(degrees * PI / 180.0f);
  mat4x4_t result;

  // Initialize to zero
  memset(&result, 0, sizeof(mat4x4_t));

  result.a = cos_r;
  result.c = -sin_r;
  result.f = 1;
  result.i = sin_r;
  result.k = cos_r;
  result.p = 1;

  return result;
}

mat4x4_t mat4x4_translate_vector3(vector3_t *v) {
  mat4x4_t result;

  // Initialize to identity matrix
  result = mat4x4_identity();

  // Set translation components
  result.d = v->x;
  result.h = v->y;
  result.l = v->z;

  return result;
}

mat4x4_t mat4x4_scale(float scale) {
  mat4x4_t result;

  // Initialize to zero
  memset(&result, 0, sizeof(mat4x4_t));

  // Set scaling components
  result.a = scale;
  result.f = scale;
  result.k = scale;
  result.p = 1.0f;

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

void draw_filled_triangle(vector2_t *p0, vector2_t *p1, vector2_t *p2,
                          unsigned char color) {
  vector2_t p0_local = *p0;
  vector2_t p1_local = *p1;
  vector2_t p2_local = *p2;
  int m, x, y;
  float *x_left, *x_right, *h_left, *h_right;
  float x01[SH], x12[SH], x02[SH], x012[SH];
  float h01[SH], h12[SH], h02[SH], h012[SH];
  int x01_len, x12_len, x02_len;
  int h01_len, h12_len, h02_len;

  // Add bounds checking for y-coordinates
  if (abs(p0_local.y) >= SH || abs(p1_local.y) >= SH || abs(p2_local.y) >= SH) {
    return;
  }

  // Sort vertices by y-coordinate (p0 lowest, p2 highest)
  if (p1_local.y < p0_local.y) {
    SWAP_VECTOR2(p1_local, p0_local);
  }

  if (p2_local.y < p0_local.y) {
    SWAP_VECTOR2(p2_local, p0_local);
  }

  if (p2_local.y < p1_local.y) {
    SWAP_VECTOR2(p2_local, p1_local);
  }

  // Interpolate x-coordinates and depth (h) for the three sides of the triangle
  interpolate(p0_local.y, p0_local.x, p1_local.y, p1_local.x, &x01, &x01_len);
  interpolate(p0_local.y, p0_local.h, p1_local.y, p1_local.h, &h01, &h01_len);

  interpolate(p1_local.y, p1_local.x, p2_local.y, p2_local.x, &x12, &x12_len);
  interpolate(p1_local.y, p1_local.h, p2_local.y, p2_local.h, &h12, &h12_len);

  interpolate(p0_local.y, p0_local.x, p2_local.y, p2_local.x, &x02, &x02_len);
  interpolate(p0_local.y, p0_local.h, p2_local.y, p2_local.h, &h02, &h02_len);

  // Concatenate the short sides (p0-p1 + p1-p2)
  concat_sides(&x01, x01_len, &x12, x12_len, &x012);
  concat_sides(&h01, h01_len, &h12, h12_len, &h012);

  // Determine which is the left and right side
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

  // Fill the triangle scanline by scanline
  for (y = p0_local.y; y <= p2_local.y; y++) {
    int xl = x_left[y - p0_local.y];
    int xr = x_right[y - p0_local.y];
    float h_segment[SW];
    int h_segment_len;

    // Interpolate depth values along the scanline
    interpolate(xl, h_left[y - p0_local.y], xr, h_right[y - p0_local.y],
                &h_segment, &h_segment_len);

    // Draw each pixel if it's closer than what's already there
    for (x = xl; x <= xr; x++) {
      if (is_closer_pixel(x, y, h_segment[x - xl])) {
        set_pixel(x, y, color);
      }
    }
  }
}

vector3_t compute_triangle_normal(vector3_t *v0, vector3_t *v1, vector3_t *v2) {
  vector3_t v0v1 = vec3_sub(v0, v1);
  vector3_t v0v2 = vec3_sub(v0, v2);

  return vec3_cross(&v0v1, &v0v2);
}

float compute_illumination(vector3_t *vertex, vector3_t *normal,
                           camera_t *camera, light_t *lights, int light_len) {
  float illumination = 0.0f, cos_alpha;
  vector3_t vl;
  int l;

  for (l = 0; l < light_len; l++) {
    light_t *curr_light = &lights[l];

    if (curr_light->type == AMBIENT_LIGHT) {
      illumination += curr_light->intensity;
      continue;
    }

    if (curr_light->type == DIRECTIONAL_LIGHT) {
      vector4_t light_vector = vec4_make(
          curr_light->vector.x, curr_light->vector.y, curr_light->vector.z, 1);
      mat4x4_t camera_matrix = mat4x4_transpose(&camera->orientation);
      vector4_t rotated_light =
          mat4x4_mul_vector4(&camera_matrix, &light_vector);

      vl = vec3_make(rotated_light.x, rotated_light.y, rotated_light.z);
    } else if (curr_light->type == POINT_LIGHT) {
      vector4_t light_vector = vec4_make(
          curr_light->vector.x, curr_light->vector.y, curr_light->vector.z, 1);
      vector3_t camera_pos_inv = vec3_mul(&camera->position, -1);
      mat4x4_t camera_translation = mat4x4_translate_vector3(&camera_pos_inv);
      mat4x4_t camera_transposed = mat4x4_transpose(&camera->orientation);
      mat4x4_t camera_matrix =
          mat4x4_mul_mat4x4(&camera_transposed, &camera_translation);
      vector4_t transformed_light =
          mat4x4_mul_vector4(&camera_matrix, &light_vector);
      vector3_t transformed_light_v3 = vec3_make(
          transformed_light.x, transformed_light.y, transformed_light.z);
      vector3_t vertex_inv = vec3_mul(vertex, -1);
      vl = vec3_add(&vertex_inv, &transformed_light_v3);
    }

    cos_alpha = vec3_dot(normal, &vl) / (vec3_len(&vl) * vec3_len(normal));
    if (cos_alpha > 0) {
      illumination += cos_alpha * curr_light->intensity;
    }
  }

  return illumination;
}

vector2_t viewport_to_canvas(float x, float y) {
  vector2_t result;
  float aspect_ratio = 1.2f; // Account for having non-square pixels

  result.x = x * SW / VIEWPORT_SIZE;
  result.y = (y * aspect_ratio) * SH / VIEWPORT_SIZE;

  return result;
}

vector2_t project_vertex(vector4_t *v) {
  vector2_t result;

  // Prevent a divide by zero
  if (v->z < 0.1f || v->z == 0) {
    v->z = 0.1f;
  }

  // Project 3D coordinates to 2D
  result = viewport_to_canvas(v->x * PROJ_PLANE_Z / v->z,
                              v->y * PROJ_PLANE_Z / v->z);

  // Store 1/z for depth buffer
  result.h = 1.0f / v->z;

  return result;
}

void setup_instance_transform(instance_t *i) {
  mat4x4_t scale = mat4x4_scale(i->scale);
  mat4x4_t translation = mat4x4_translate_vector3(&i->position);
  mat4x4_t composed = mat4x4_mul_mat4x4(&i->orientation, &scale);

  i->transform = mat4x4_mul_mat4x4(&translation, &composed);
}

int is_clipped(plane_t *planes, model_t *model, float scale,
               mat4x4_t *transform) {
  int p;
  float radius = model->radius * scale;
  vector4_t center =
      vec4_make(model->center.x, model->center.y, model->center.z, 1);
  center = mat4x4_mul_vector4(transform, &center);

  for (p = 0; p < CLIPPING_PLANE_COUNT; p++) {
    vector3_t center_v3 = vec3_make(center.x, center.y, center.z);
    float distance =
        vec3_dot(&planes[p].normal, &center_v3) + planes[p].distance;
    if (distance < -radius) {
      return 1;
    }
  }

  return 0;
}

void render_triangle(triangle_t *triangle, vector3_t *transformed_verts,
                     vector2_t *proj_verts, unsigned char color,
                     camera_t *camera, light_t *lights, int light_len) {
  vector3_t v0 = transformed_verts[triangle->vertex_index[0]];
  vector3_t v1 = transformed_verts[triangle->vertex_index[1]];
  vector3_t v2 = transformed_verts[triangle->vertex_index[2]];
  vector3_t normal = compute_triangle_normal(&v0, &v1, &v2);
  vector3_t vert_to_camera = vec3_mul(&v0, -1);
  vector3_t center, normal_inv;
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

void render_model(model_t *model, unsigned char color, mat4x4_t *matrix,
                  camera_t *camera, light_t *lights, int light_len) {
  int v, t;
  vector3_t *transformed_verts;
  vector2_t *projected_verts;

  /* Allocate memory for transformed and projected vertices */
  transformed_verts =
      (vector3_t *)malloc(model->vert_count * sizeof(vector3_t));
  projected_verts = (vector2_t *)malloc(model->vert_count * sizeof(vector2_t));

  if (!transformed_verts || !projected_verts) {
    /* Handle allocation failure */
    if (transformed_verts)
      free(transformed_verts);
    if (projected_verts)
      free(projected_verts);
    return;
  }

  // Transform and project all vertices
  for (v = 0; v < model->vert_count; v++) {
    vector3_t vertex = model->vertices[v];
    vector4_t vertex_h = vec4_make(vertex.x, vertex.y, vertex.z, 1);
    vector4_t vertex_mult = mat4x4_mul_vector4(matrix, &vertex_h);

    transformed_verts[v] =
        vec3_make(vertex_mult.x, vertex_mult.y, vertex_mult.z);
    projected_verts[v] = project_vertex(&vertex_mult);
  }

  // Render all triangles
  for (t = 0; t < model->tri_count; t++) {
    render_triangle(&model->triangles[t], transformed_verts, projected_verts,
                    color, camera, lights, light_len);
  }

  /* Free allocated memory */
  free(transformed_verts);
  free(projected_verts);
}

void render_scene(camera_t *camera, instance_t *instances, int instance_len,
                  light_t *lights, int light_len) {
  int i;
  mat4x4_t cameraMatrix, m1, m2;
  vector3_t neg_camera_pos;

  // Create camera transformation matrix
  neg_camera_pos = vec3_mul(&camera->position, -1);
  m1 = mat4x4_transpose(&camera->orientation);
  m2 = mat4x4_translate_vector3(&neg_camera_pos);
  cameraMatrix = mat4x4_mul_mat4x4(&m1, &m2);

  // Process each instance in the scene
  for (i = 0; i < instance_len; i++) {
    mat4x4_t transform =
        mat4x4_mul_mat4x4(&cameraMatrix, &instances[i].transform);

    // Check if the object is within the frustum
    int status = is_clipped(camera->clipping_planes, instances[i].model,
                            instances[i].scale, &transform);

    if (status == 0) {
      render_model(instances[i].model, instances[i].color, &transform, camera,
                   lights, light_len);
    }
  }
}

model_t *generate_sphere(int divs) {
  model_t *sphere;
  vector3_t *vertices;
  triangle_t *triangles;
  float delta_angle;
  int d, i;
  int vertex_count, triangle_count;
  int i0, i1, i2;
  float y, radius;

  vertex_count = divs * (divs + 1);
  triangle_count = divs * divs * 2;

  sphere = (model_t *)malloc(sizeof(model_t));
  vertices = (vector3_t *)malloc(vertex_count * sizeof(vector3_t));
  triangles = (triangle_t *)malloc(triangle_count * sizeof(triangle_t));

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

model_t *generate_plane(int divs, float size) {
  model_t *plane;
  vector3_t *vertices;
  triangle_t *triangles;
  int vertex_count, triangle_count;
  int x, z, i;
  float half_size, step;

  vertex_count = divs * divs;
  triangle_count = (divs - 1) * (divs - 1) * 2;

  plane = (model_t *)malloc(sizeof(model_t));
  vertices = (vector3_t *)malloc(vertex_count * sizeof(vector3_t));
  triangles = (triangle_t *)malloc(triangle_count * sizeof(triangle_t));

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

void free_model(model_t *model) {
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
  camera_t camera;
  model_t *sphere = generate_sphere(13);
  model_t *plane = generate_plane(20, 30.0f);
  instance_t instances[4];
  light_t lights[3];

  // Setup lights
  lights[0].type = AMBIENT_LIGHT;
  lights[0].intensity = 0.3f;

  lights[1].type = DIRECTIONAL_LIGHT;
  lights[1].intensity = 0.2f;
  lights[1].vector = vec3_make(-1, 0, 4);

  lights[2].type = POINT_LIGHT;
  lights[2].intensity = 0.6f;
  lights[2].vector = vec3_make(-2, 2, -2);

  // Setup scene objects
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

  // Setup view frustum clipping planes
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

  // Initialize graphics
  set_mode(0x13);
  init_palette();
  init_depth_buffer();
  clear_screen(shade_color(7, 31) /* White */);

  // Render the scene
  render_scene(&camera, instances, 4, lights, 3);

  // Present the frame and wait for a keypress
  present_buffer();
  getch();

  // Cleanup and exit
  set_mode(0x03);
  free_model(sphere);
  free_model(plane);

  return 0;
}
