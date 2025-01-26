#include <dos.h>
#include <i86.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>

#define INF 32767 // Infinity for 16 bit compilers
#define SW 320    // Screen width
#define SH 200    // Screen height

#define SWAP_VECTOR2(a, b) do { \
    struct vector2 temp = (a);  \
    (a) = (b);                 \
    (b) = temp;                \
} while(0)

struct vector2
{
    int x;
    int y;
};

struct vector3
{
    int x;
    int y;
    int z;
};

void set_mode(unsigned char mode)
{
    union REGS regs;

    regs.h.ah = 0;
    regs.h.al = mode;

    int86(0x10, &regs, &regs);
}

void set_pixel(int x, int y, char color)
{
    char far *screen = MK_FP(0xA000, 0);
    int ax, ay;

    ax = SW / 2 + x;
    ay = SH / 2 - y;

    if (ax < 0 || ax >= SW || ay < 0 || ay >= SH)
    {
        return;
    }

    screen[ay * SW + ax] = color;
}

struct vector3 vec3_add(struct vector3 *a, struct vector3 *b)
{
    struct vector3 result;

    result.x = a->x + b->x;
    result.y = a->y + b->y;
    result.z = a->z + b->z;

    return result;
}

struct vector3 vec3_sub(struct vector3 *a, struct vector3 *b)
{
    struct vector3 result;

    result.x = b->x - a->x;
    result.y = b->y - a->y;
    result.z = b->z - a->z;

    return result;
}

struct vector3 vec3_mul(struct vector3 *v, float s)
{
    struct vector3 result;

    result.x = v->x * s;
    result.y = v->y * s;
    result.z = v->z * s;

    return result;
}

struct vector3 vec3_div(struct vector3 *v, float d)
{
    struct vector3 result;

    result.x = v->x / d;
    result.y = v->y / d;
    result.z = v->z / d;

    return result;
}

float vec3_len(struct vector3 *v)
{
    return sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
}

float vec3_dot(struct vector3 *a, struct vector3 *b)
{
    return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}

int *interpolate(int i0, int d0, int i1, int d1, int *size)
{
    int *values = 0;
    float a, d;
    int i, sizei1_i0, idx = 0;

    if (i0 == i1) {
        values = malloc(sizeof(int));

        if (values == 0) {
            return 0;
        }
        
        values[0] = d0;

        *size = 1;
        
        return values;
    }

    sizei1_i0 = abs(i0 - i1);
    values = malloc(sizei1_i0 * sizeof(int));

    if (values == 0) {
        return 0;
    }

    a = (d1 - d0) / (float)(i1 - i0);
    d = d0;

    for (i = i0; i <= i1; i++) {
        values[idx] = d;
        d = d + a;
        idx++;
    }

    *size = sizei1_i0;

    return values;
}

void draw_line(struct vector2 *p0, struct vector2 *p1, unsigned char color)
{
    struct vector2 p0_local = *p0;
    struct vector2 p1_local = *p1;
    float a;
    int x, y, size, *slope_values;

    if (abs(p1_local.x - p0_local.x) > abs(p1_local.y - p0_local.y)) {
        if (p0_local.x > p1_local.x) {
            SWAP_VECTOR2(p0_local, p1_local);
        }

        slope_values = interpolate(p0_local.x, p0_local.y, 
            p1_local.x, p1_local.y, &size);

        if (slope_values == 0) {
            return;
        }

        for (x = p0_local.x; x < p1_local.x; x++) {
            set_pixel(x, slope_values[x - p0_local.x], color);
        }
    } else {
        if (p0_local.y > p1_local.y) {
            SWAP_VECTOR2(p0_local, p1_local);
        }

        slope_values = interpolate(p0_local.y, p0_local.x, 
            p1_local.y, p1_local.x, &size);

        if (slope_values == 0) {
            return;
        }

        for (y = p0_local.y; y < p1_local.y; y++) {
            set_pixel(slope_values[y - p0_local.y], y, color);
        }
    }

    free(slope_values);
}

void draw_wireframe_triangle(struct vector2 *p0, struct vector2 *p1, struct vector2 *p2,
    unsigned char color) 
{
    draw_line(p0, p1, color);
    draw_line(p1, p2, color);
    draw_line(p2, p0, color);
}

int *concat_sides(int *x01, int size01, int *x12, int size12, int *size012) {
    int size = ((size01 - 1) + size12);
    int *concated = malloc(size * sizeof(int));
    int i;

    for (i = 0; i <= size01 - 1; i++) {
        concated[i] = x01[i];
    }

    for (i = 0; i <= size12; i++) {
        concated[i + (size01 - 1)] = x12[i];
    }

    *size012 = size;

    return concated;
}

void draw_filled_triangle(struct vector2 *p0, struct vector2 *p1, struct vector2 *p2,
    unsigned char color) {
    struct vector2 p0_local = *p0;
    struct vector2 p1_local = *p1;
    struct vector2 p2_local = *p2;
    int *x01, *x12, *x02, *x012;
    int size01, size12, size02, size012, m;
    int x, y, *x_left, *x_right;

    if (p1_local.y < p0_local.y) {
        SWAP_VECTOR2(p1_local, p0_local);
    }

    if (p2_local.y < p0_local.y) {
        SWAP_VECTOR2(p2_local, p0_local);
    }

    if (p2_local.y < p1_local.y) {
        SWAP_VECTOR2(p2_local, p1_local);
    }

    x01 = interpolate(p0_local.y, p0_local.x, p1_local.y, p1_local.x, &size01);
    x12 = interpolate(p1_local.y, p1_local.x, p2_local.y, p2_local.x, &size12);
    x02 = interpolate(p0_local.y, p0_local.x, p2_local.y, p2_local.x, &size02);
    x012 = concat_sides(x01, size01, x12, size12, &size012);
    m = size02 / 2;

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

    free(x01);
    free(x12);
    free(x02);
    free(x012);
}

int main(void)
{
    struct vector2 p0 = { -50, -62 };
    struct vector2 p1 = { 50, 12 };
    struct vector2 p2 = { 5, 62 };

    set_mode(0x13);

    draw_filled_triangle(&p0, &p1, &p2, 2);
    draw_wireframe_triangle(&p0, &p1, &p2, 15);

    getch();

    set_mode(0x03);

    return 0;
}
