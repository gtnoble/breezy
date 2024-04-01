#ifndef BREEZY_MATH
#define BREEZY_MATH

#include <stddef.h>
#include <stdlib.h>

#include "arena.h"

#define BREEZY_CONST_VECTOR(x, y, z) (const Vector3D) {x, y, z}
#define BREEZY_INCHES_TO_MM(inches) inches * 25.4

typedef double Vector3D[3];

static const Vector3D x_unit_vector = {1.0, 0.0, 0.0};
static const Vector3D y_unit_vector = {0.0, 1.0, 0.0};
static const Vector3D z_unit_vector = {0.0, 0.0, 1.0};
static const Vector3D zeros_vector = {0.0, 0.0, 0.0};

typedef Vector3D Matrix[3];

typedef struct {
    size_t num_vertices;
    Vector3D *vertices[];
} Points;

#define BREEZY_DOUBLE_ARRAY(arena, ...) \
    bzy_new_double_array((double []) {__VA_ARGS__}, sizeof((double []) {__VA_ARGS__}) / sizeof(double), arena)

typedef struct {
    size_t num_elements;
    double elements[];
} DoubleArray;

double bzy_degrees_to_radians(double degrees);

DoubleArray *bzy_new_double_array(const double elements[], size_t num_elements, Arena *arena);

Vector3D *bzy_vector_normalizen(const Vector3D vector, Arena *arena);
Vector3D *bzy_vector_addn(const Vector3D v1, const Vector3D v2, Arena *arena);
Vector3D *bzy_vector_subn(const Vector3D v1, const Vector3D v2, Arena *arena);
Vector3D *bzy_scalar_vector_multn(double scalar, const Vector3D vector, Arena *arena);
Vector3D *bzy_hole_clearance(const Vector3D hole_direction, const Vector3D plane_normal, double radius, Arena *arena);
double *bzy_flatten_vectors(Vector3D * const vectors[], size_t num_vectors, Arena *arena);


Points *bzy_vectors_to_points(const Vector3D vectors[], size_t num_vectors, Arena *arena);
Points *bzy_cartesian_product(
    const DoubleArray *x_values,
    const DoubleArray *y_values,
    const DoubleArray *z_values,
    Arena *arena
);
Points *bzy_translate_points(const Points *points, const Vector3D direction, Arena *arena);
Points *bzy_rotate_points(const Points *points, const Vector3D axis, double angle, Arena *arena);
Points *bzy_polygonalize(const Points *points, Arena *arena);
Vector3D *bzy_polygon_normal(Points *polygon, Arena *arena);

#endif
