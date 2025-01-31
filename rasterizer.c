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
  struct vector3 position;
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

struct vector2 project_vertex(struct vector3 *v) {
  return viewport_to_canvas(v->x * PROJ_PLANE_Z / v->z,
                            v->y * PROJ_PLANE_Z / v->z);
}

void render_object(struct vector3 *vertices, int vert_len,
                   struct triangle *triangles, int tri_len,
                   struct vector3 *position) {
  int v, t;
  struct vector2 proj_verts[8];

  if (vert_len > 8) {
    return;
  }

  for (v = 0; v < vert_len; v++) {
    struct vector3 translated = vec3_add(&vertices[v], position);
    proj_verts[v] = project_vertex(&translated);
  }

  for (t = 0; t < tri_len; t++) {
    render_triangle(&triangles[t], proj_verts);
  }
}

void render_instance(struct instance *instance) {
  render_object(instance->model->vertices, instance->model->vert_count,
                instance->model->triangles, instance->model->tri_count,
                &instance->position);
}

void render_scene(struct instance *instances, int len) {
  int i;

  for (i = 0; i < len; i++) {
    render_instance(&instances[i]);
  }
}

int main(void) {
  // For the shaded triangle
  struct vector2 p0 = {-50, -62, 0.3f};
  struct vector2 p1 = {50, 12, 0.1f};
  struct vector2 p2 = {5, 62, 1.0f};

  // For simple quad 3d cube
  struct vector3 vAf = {-2, -0.5f, 5};
  struct vector3 vBf = {-2, 0.5f, 5};
  struct vector3 vCf = {-1, 0.5f, 5};
  struct vector3 vDf = {-1, -0.5f, 5};
  struct vector3 vAb = {-2, -0.5f, 6};
  struct vector3 vBb = {-2, 0.5f, 6};
  struct vector3 vCb = {-1, 0.5f, 6};
  struct vector3 vDb = {-1, -0.5f, 6};
  struct vector2 pAf = project_vertex(&vAf);
  struct vector2 pBf = project_vertex(&vBf);
  struct vector2 pCf = project_vertex(&vCf);
  struct vector2 pDf = project_vertex(&vDf);
  struct vector2 pAb = project_vertex(&vAb);
  struct vector2 pBb = project_vertex(&vBb);
  struct vector2 pCb = project_vertex(&vCb);
  struct vector2 pDb = project_vertex(&vDb);

  // For the triangle based cube
  int i;
  struct vector3 position = {4, -1, 15};
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

  set_mode(0x13);

  init_palette();

  // The green shaded triangle
  draw_shaded_triangle(&p0, &p1, &p2, 1);
  draw_wireframe_triangle(&p0, &p1, &p2, shade_color(7, 31));

  // The simple quad cube
  draw_line(&pAf, &pBf, shade_color(2, 31));
  draw_line(&pBf, &pCf, shade_color(2, 31));
  draw_line(&pCf, &pDf, shade_color(2, 31));
  draw_line(&pDf, &pAf, shade_color(2, 31));
  draw_line(&pAb, &pBb, shade_color(0, 31));
  draw_line(&pBb, &pCb, shade_color(0, 31));
  draw_line(&pCb, &pDb, shade_color(0, 31));
  draw_line(&pDb, &pAb, shade_color(0, 31));
  draw_line(&pAf, &pAb, shade_color(1, 31));
  draw_line(&pBf, &pBb, shade_color(1, 31));
  draw_line(&pCf, &pCb, shade_color(1, 31));
  draw_line(&pDf, &pDb, shade_color(1, 31));

  // Render the model
  render_object(vertices, 8, triangles, 12, &position);

  // Render instance of the cube model
  the_cube.name = "cool cube";
  the_cube.vert_count = 8;
  the_cube.tri_count = 12;
  the_cube.vertices = vertices;
  the_cube.triangles = triangles;
  cube_instances[0].model = &the_cube;
  cube_instances[0].position.x = -1;
  cube_instances[0].position.y = -2;
  cube_instances[0].position.z = 8;
  cube_instances[1].model = &the_cube;
  cube_instances[1].position.x = 1;
  cube_instances[1].position.y = 2;
  cube_instances[1].position.z = 8;

  render_scene(cube_instances, 2);

  getch();

  set_mode(0x03);

  return 0;
}
