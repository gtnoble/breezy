#define __POSIX_C_SOURCE 202403L

#include <math.h>
#include <stdalign.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <float.h>

#include "bzy_math.h"
#include "arena.h"

double bzy_degrees_to_radians(double degrees) {
    return (degrees / 180.0) * M_PI;
}

Vector3D *new_vector(Arena *arena) {
    return arena_allocate_aligned(sizeof(Vector3D), alignof(Vector3D), arena);
}

void bzy_vector_copy(const Vector3D vector, Vector3D copy) {
    assert(vector != NULL);
    assert(copy != NULL);

    for (int i = 0; i < 3; i++) {
        copy[i] = vector[i];
    }
}

double *bzy_flatten_vectors(Vector3D * const vectors[], size_t num_vectors, Arena *arena) {
    double *flattened = arena_allocate_aligned(
        sizeof(double) * 3 * num_vectors, alignof(double), arena
    );

    size_t flattened_index = 0;
    for (size_t i = 0; i < num_vectors; i++) {
        for (size_t j = 0; j < 3; j++) {
            flattened[flattened_index] = (*vectors[i])[j];
            flattened_index++;
        }
    }

    return flattened;
}

Vector3D *duplicate_vector(const Vector3D vector, Arena *arena) {
    assert(vector != NULL);
    assert(arena != NULL);

    Vector3D *duplicated_vector = new_vector(arena);
    bzy_vector_copy(vector, *duplicated_vector);
    return duplicated_vector;
}

bool bzy_vector_equal(const Vector3D v1, const Vector3D v2) {
    for (int i = 0; i < 3; i++) {
        if (v1[i] != v2[i]) {
            return false;
        }
    }
    return true;
}

void bzy_vector_add(const Vector3D v1, const Vector3D v2, Vector3D result) {
    for (int i = 0; i < 3; i++) {
        result[i] = v1[i] + v2[i];
    }
}

Vector3D *bzy_vector_addn(const Vector3D v1, const Vector3D v2, Arena *arena) {
    Vector3D *result = new_vector(arena);
    bzy_vector_add(v1, v2, *result);
    return result;
}

void bzy_vector_sub(const Vector3D v1, const Vector3D v2, Vector3D result) {
    for (int i = 0; i < 3; i++) {
        result[i] = v1[i] - v2[i];
    }
}

Vector3D *bzy_vector_subn(const Vector3D v1, const Vector3D v2, Arena *arena) {
    Vector3D *result = new_vector(arena);
    bzy_vector_sub(v1, v2, *result);
    return result;
}

void bzy_scalar_vector_mult(double scalar, const Vector3D vector, Vector3D result) {
    for (int i = 0; i < 3; i++) {
        result[i] = scalar * vector[i];
    }
}

Vector3D *bzy_scalar_vector_multn(double scalar, const Vector3D vector, Arena *arena) {
    Vector3D *result = new_vector(arena);
    bzy_scalar_vector_mult(scalar, vector, *result);
    return result;
}

void bzy_vector_negate(const Vector3D vector, Vector3D result) {
    bzy_scalar_vector_mult(-1.0, vector, result);
}

void bzy_hadamard_product(const Vector3D v1, const Vector3D v2, Vector3D result) {
    for (int i = 0; i < 3; i++) {
        result[i] = v1[i] * v2[2];
    }
}

Vector3D *bzy_hadamard_productn(const Vector3D v1, const Vector3D v2, Arena *arena) {
    Vector3D *result = new_vector(arena);
    bzy_hadamard_product(v1, v2, *result);
    return result;
}

double bzy_vector_length(const Vector3D vector) {
    double sum_squares = 0;
    for (int i = 0; i < 3; i++) {
        sum_squares += pow(vector[i], 2);
    }
    return sqrt(sum_squares);
}

void bzy_vector_normalize(const Vector3D vector, Vector3D result) {
    double length = bzy_vector_length(vector);
    assert(length > 0);
    for (int i = 0; i < 3; i++) {
        result[i] = vector[i] / length;
    }
}

Vector3D *bzy_vector_normalizen(const Vector3D vector, Arena *arena) {
    Vector3D *result = new_vector(arena);
    bzy_vector_normalize(vector, *result);
    return result;
}

double bzy_vector_dot(const Vector3D v1, const Vector3D v2) {
    double sum = 0;
    for (int i = 0; i < 3; i++) {
        sum += v1[i] * v2[i];
    }
    return sum;
}


void bzy_matrix_vector_mult(const Matrix matrix, const Vector3D vector, Vector3D result) {
    for (int i = 0; i < 3; i++) {
        result[i] = bzy_vector_dot(matrix[i], vector);
    }
}

void bzy_vector_cross(const Vector3D v1, const Vector3D v2, Vector3D result) {
    const Matrix matrix = {
        {0.0, -v1[2], v1[1]},
        {v1[2], 0.0, -v1[0]},
        {-v1[1], v1[0], 0.0}
    };
    bzy_matrix_vector_mult(
        matrix,
        v2, 
        result
    );
}

Vector3D *bzy_vector_crossn(const Vector3D v1, const Vector3D v2, Arena *arena) {
    Vector3D *result = new_vector(arena);
    bzy_vector_cross(v1, v2, *result);
    return result;
}

void bzy_clear_vector(Vector3D vector) {
    for (int i = 0; i < 3; i++) {
        vector[i] = 0.0;
    }
}

Vector3D *bzy_rotate_vector(
    const Vector3D axis_vector, 
    double angle, 
    const Vector3D vector, 
    Arena *arena
) {

    Vector3D normalized_axis_vector;
    bzy_vector_normalize(axis_vector, normalized_axis_vector);

    Arena scratch = make_arena(1000000);

    Vector3D *result = bzy_vector_addn(
        *bzy_vector_addn(
            *bzy_scalar_vector_multn(cos(angle), vector, &scratch),
            *bzy_scalar_vector_multn(
                sin(angle), 
                *bzy_vector_crossn(
                    normalized_axis_vector, 
                    vector, 
                    &scratch
                ), 
                &scratch
            ),
            &scratch
        ),
        *bzy_scalar_vector_multn(
            bzy_vector_dot(normalized_axis_vector, vector) * (1 - cos(angle)), 
            normalized_axis_vector, 
            &scratch
        ),
        arena
    );

    free_arena(&scratch);

    return result;
}

double bzy_point_distance(const Vector3D p1, const Vector3D p2) {
    Vector3D difference;
    bzy_vector_sub(p1, p2, difference);
    return bzy_vector_length(difference);
}

Points *new_points(size_t num_points, Arena *arena) {
    Points *points = arena_allocate_aligned(
        sizeof(Points) + sizeof(Vector3D *) * num_points, alignof(Points), arena
    );
    points->num_vertices = num_points;
    return points;
}

Points *bzy_vectors_to_points(const Vector3D vectors[], size_t num_vectors, Arena *arena) {
    Points *points = new_points(num_vectors, arena);
    for (size_t i = 0; i < num_vectors; i++) {
        points->vertices[i] = duplicate_vector(vectors[i], arena);
    }
    return points;
}

BzyDoubles *bzy_new_doubles(size_t num_elements, Arena *arena) {
    assert(arena != NULL);
    
    BzyDoubles *array = arena_allocate_aligned(
        sizeof(BzyDoubles), 
        alignof(BzyDoubles), 
        arena
    );

    array->length = num_elements;

    array->elements = arena_allocate_aligned(
        sizeof(double) * num_elements, 
        alignof(double), 
        arena
    );

    return array;
}

BzyDoubles bzy_make_doubles(double elements[], size_t num_elements) {
    BzyDoubles array = {
        .elements = elements,
        .length = num_elements
    };
    return array;
}

double *bzy_doubles_get(size_t index, BzyDoubles doubles) {
    assert(index < doubles.length);
    return &doubles.elements[index];
}

Points *bzy_cartesian_product(
    const BzyDoubles x_values,
    const BzyDoubles y_values,
    const BzyDoubles z_values,
    Arena *arena
) {
    assert(x_values.length > 0);
    assert(y_values.length > 0);
    assert(z_values.length > 0);
    assert(arena != NULL);

    size_t num_points = x_values.length * y_values.length * z_values.length;
    Points *product = new_points(
        num_points,
        arena
    );

    size_t index = 0;
    for (size_t i = 0; i < x_values.length; i++) {
        for (size_t j = 0; j < y_values.length; j++) {
            for (size_t k = 0; k < z_values.length; k++) {
                assert(index < num_points);
                product->vertices[index] = duplicate_vector(
                    (Vector3D) {
                        *bzy_doubles_get(i, x_values),
                        *bzy_doubles_get(j, y_values),
                        *bzy_doubles_get(k, z_values)
                    }, 
                    arena
                );
                index++;
            }
        }
    }

    return product;
}

Points *duplicate_points(Points *points, Arena *arena) {
    assert(points != NULL);
    assert(points->vertices != NULL);
    assert(arena != NULL);

    Points *duplicated_points = new_points(points->num_vertices, arena);
    for (size_t i = 0; i < points->num_vertices; i++) {
        duplicated_points->vertices[i] = duplicate_vector(*points->vertices[i], arena);
    }
    return duplicated_points;
}

Points *bzy_translate_points(const Points *points, const Vector3D direction, Arena *arena) {
    Points *translated_points = new_points(points->num_vertices, arena);
    
    for (size_t i = 0; i < points->num_vertices; i++) {
        translated_points->vertices[i] = 
            bzy_vector_addn(direction, *points->vertices[i], arena);
    }

    return translated_points;
}

Points *bzy_rotate_points(const Points *points, const Vector3D axis, double angle, Arena *arena) {
    Points *rotated_points = new_points(points->num_vertices, arena);

    for (size_t i = 0; i < points->num_vertices; i++) {
        rotated_points->vertices[i] = bzy_rotate_vector(
            axis, 
            angle, 
            *points->vertices[i], 
            arena
        );
    }

    return rotated_points;
}


Points *bzy_flip_axis(const Vector3D axis, const Points *vertices, Arena *arena) {
    Points *flipped_vertices = new_points(vertices->num_vertices, arena);

    Vector3D flipped_axis;
    bzy_vector_negate(axis, flipped_axis);

    for (size_t i = 0; i < vertices->num_vertices; i++) {
        flipped_vertices->vertices[i] = bzy_hadamard_productn(
            *vertices->vertices[i], 
            flipped_axis, 
            arena
        );
    }
    return flipped_vertices;
}

double bzy_path_length(Points *path) {
    double length = 0;
    for (size_t i = 0; i < path->num_vertices - 1; i++) {
        length += bzy_point_distance(*path->vertices[i], *path->vertices[i+1]);
    }
    return length;
}

bool increment_indices(size_t num_indices, size_t *indices) {

    bool increment_index = true;
    for (size_t i = 0; i < num_indices; i++) {
        indices[i] += increment_index;
        if (indices[i] == num_indices) {
            indices[i] = 0;
            increment_index = true;
        }
        else {
            increment_index = false;
        }
    }

    if (increment_index) {
        return false;
    }
    else {
        return true;
    }
}

//@todo this can be done in O(n) time with a boolean valued array
bool check_indices_unique(size_t num_indices, const size_t *indices) {
    for (size_t i = 0; i < num_indices; i++) {
        for (size_t j = 0; j < num_indices; j++) {
            if (i != j && indices[i] == indices[j]) {
                return false;
            }
        }
    }
    return true;
}

bool permute_indices(size_t num_indices, size_t *indices) {
    while (increment_indices(num_indices, indices)) {
        if (check_indices_unique(num_indices, indices)) {
            return true;
        }
    }
    return false;
}

void apply_permutation(const Points *points, Points *permuted_points, size_t permutation[]) {
    for (size_t i = 0; i < points->num_vertices; i++) {
        size_t source_index = permutation[i];
        assert(source_index < points->num_vertices);

        permuted_points->vertices[source_index] = points->vertices[i];
    }
}

Points *bzy_polygonalize(const Points *points, Arena *arena) {
    assert(points != NULL);
    assert(arena != NULL);

    Arena permutation_arena = make_arena(1000000);

    size_t *permutation = arena_allocate_zeros(
        sizeof(size_t) * points->num_vertices, 
        &permutation_arena
    );
    Points *current_shortest_path = new_points(
        points->num_vertices, 
        &permutation_arena
    );
    Points *candidate_shortest_path = new_points(
        points->num_vertices, 
        &permutation_arena
    );

    double current_shortest_distance = DBL_MAX;

    while (permute_indices(points->num_vertices, permutation)) {
        apply_permutation(
            points, 
            candidate_shortest_path, 
            permutation
        );
        double candidate_path_distance = bzy_path_length(candidate_shortest_path);
        assert(candidate_path_distance >= 0);

        if (candidate_path_distance <= current_shortest_distance) {
            Points *swap;
            swap = current_shortest_path;
            current_shortest_path = candidate_shortest_path;
            candidate_shortest_path = swap;
            current_shortest_distance = candidate_path_distance;
        }
    }

    Points *shortest_path = duplicate_points(current_shortest_path, arena);
    free_arena(&permutation_arena);

    return shortest_path;
}

Vector3D *bzy_polygon_normal(const Points *polygon, Arena *arena) {
    assert(polygon != NULL);
    assert(arena != NULL);

    Arena scratch = make_arena(1000000);

    Vector3D *v1 = bzy_vector_subn(*polygon->vertices[1], *polygon->vertices[0], &scratch);
    assert(bzy_vector_length(*v1) > 0);

    Vector3D *v2 = bzy_vector_subn(*polygon->vertices[2], *polygon->vertices[0], &scratch);
    assert(bzy_vector_length(*v2) > 0);

    Vector3D *normal = bzy_vector_normalizen(
        *bzy_vector_crossn(
            *v1,
            *v2,
            &scratch), 
        arena
    );
    assert(bzy_vector_length(*normal) > 0);

    free_arena(&scratch);

    return normal;
}

Vector3D *bzy_hole_clearance(const Vector3D hole_direction, const Vector3D plane_normal, double radius, Arena *arena) {
    Arena scratch = make_arena(10000000);

    assert(radius >= 0);
    assert(bzy_vector_length(hole_direction) > 0.0);
    assert(bzy_vector_length(plane_normal) > 0.0);

    Vector3D *unit_normal = bzy_vector_normalizen(plane_normal, &scratch);
    Vector3D *unit_direction = bzy_vector_normalizen(hole_direction, &scratch);

    double direction_normal_component = bzy_vector_dot(*unit_normal, *unit_direction);
    
    Vector3D *clearance = bzy_scalar_vector_multn(
        fabs( 
            sqrt(1 - pow(direction_normal_component, 2)) / direction_normal_component
        ) *
        radius, 
        *unit_direction, 
        arena
    );

    free_arena(&scratch);

    return clearance;
}



