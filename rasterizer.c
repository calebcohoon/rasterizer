#include <dos.h>
#include <i86.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>

#define INF 32767 // Infinity for 16 bit compilers
#define SW 320    // Screen width
#define SH 200    // Screen height

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

int *interpolate(int i0, int d0, int i1, int d1)
{
    int *values = 0;
    float a, d;
    int i, idx = 0;

    if (i0 == i1) {
        values = malloc(sizeof(int));
        
        values[0] = d0;
        
        return values;
    }

    values = malloc(abs(i1 - i0) * sizeof(int));
    a = (d1 - d0) / (float)(i1 - i0);
    d = d0;

    for (i = i0; i < i1; i++) {
        values[idx] = d;
        d = d + a;
        idx++;
    }

    return values;
}

void draw_line(struct vector2 *p0, struct vector2 *p1, unsigned char color)
{
    struct vector2 p0_local = *p0;
    struct vector2 p1_local = *p1;
    float a;
    int x, y, *slope_values;

    if (abs(p1_local.x - p0_local.x) > abs(p1_local.y - p0_local.y)) {
        if (p0_local.x > p1_local.x) {
            struct vector2 temp = p0_local;
            
            p0_local = p1_local;
            p1_local = temp;
        }

        slope_values = interpolate(p0_local.x, p0_local.y, p1_local.x, p1_local.y);

        for (x = p0_local.x; x < p1_local.x; x++) {
            set_pixel(x, slope_values[x - p0_local.x], color);
        }
    } else {
        if (p0_local.y > p1_local.y) {
            struct vector2 temp = p0_local;
            
            p0_local = p1_local;
            p1_local = temp;
        }

        slope_values = interpolate(p0_local.y, p0_local.x, p1_local.y, p1_local.x);

        for (y = p0_local.y; y < p1_local.y; y++) {
            set_pixel(slope_values[y - p0_local.y], y, color);
        }
    }

    free(slope_values);
}

int main(void)
{
    struct vector2 p0, p1, p2, p3;

    p0.x = -100;
    p0.y = -50;
    p1.x = 120;
    p1.y = 60;

    p2.x = -25;
    p2.y = -100;
    p3.x = 30;
    p3.y = 120;

    set_mode(0x13);

    draw_line(&p0, &p1, 15);
    draw_line(&p2, &p3, 15);

    getch();

    set_mode(0x03);

    return 0;
}
