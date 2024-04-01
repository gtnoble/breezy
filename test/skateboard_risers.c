#include <assert.h>

#include "arena.h"
#include "db.h"
#include "bzy_math.h"

const double k_mounting_hole_spacing_width = BREEZY_INCHES_TO_MM(1.625);
const double k_old_school_hole_spacing_length = BREEZY_INCHES_TO_MM(2.5);
const double k_new_school_hole_spacing_length = BREEZY_INCHES_TO_MM(2.125);
const double k_hole_diameter = BREEZY_INCHES_TO_MM(0.215);
const double k_hole_to_edge_distance = 1.5 + k_hole_diameter / 2.0;

const double k_base_width = k_mounting_hole_spacing_width + 2.0 * k_hole_to_edge_distance;
const double k_base_length = k_old_school_hole_spacing_length + 2.0 * k_hole_to_edge_distance;

Points *base_points_to_top_points(Points *base_points, double wedge_angle, double base_height, Arena *arena) {
    Arena scratch = make_arena(10000000);

    Points *top_points = bzy_translate_points(
        bzy_rotate_points(base_points, y_unit_vector, wedge_angle, &scratch), 
        *bzy_scalar_vector_multn(base_height, z_unit_vector, &scratch),
        arena
    );

    free_arena(&scratch);

    return top_points;
}

char *make_riser(double base_height, double wedge_angle, char *basename) {
    Arena scratch = make_arena(10000000);

    Points *base_vertices = bzy_cartesian_product(
        BREEZY_DOUBLE_ARRAY(&scratch, 0.0, k_base_length), 
        BREEZY_DOUBLE_ARRAY(&scratch, k_base_width / 2.0, -k_base_width / 2.0), 
        BREEZY_DOUBLE_ARRAY(&scratch, 0.0), 
        &scratch
    );
    Points *top_vertices = base_points_to_top_points(base_vertices, wedge_angle, base_height, &scratch);

    char *body = bzy_make_arb8("body.s", base_vertices, top_vertices);

    Points *base_hole_centers = bzy_translate_points(
        bzy_cartesian_product(
            BREEZY_DOUBLE_ARRAY(
                &scratch, 
                0.0, 
                k_old_school_hole_spacing_length, 
                k_new_school_hole_spacing_length
            ), 
            BREEZY_DOUBLE_ARRAY(
                &scratch,
                k_mounting_hole_spacing_width / 2.0, 
                -k_mounting_hole_spacing_width / 2.0
            ) , 
            BREEZY_DOUBLE_ARRAY(&scratch, 0.0), 
            &scratch
        ), 
        *bzy_scalar_vector_multn(
            k_hole_to_edge_distance, 
            x_unit_vector, 
            &scratch
        ), 
        &scratch
    );

    Points *top_hole_centers = base_points_to_top_points(
        base_hole_centers, wedge_angle, base_height, &scratch
    );

    const char *drills[base_hole_centers->num_vertices + 1];
    drills[base_hole_centers->num_vertices] = NULL;
    Vector3D *base_normal = bzy_polygon_normal(base_vertices, &scratch);
    Vector3D *top_normal = bzy_polygon_normal(top_vertices, &scratch);
    for (size_t i = 0; i < base_hole_centers->num_vertices; i++) {
        drills[i] = bzy_make_drill(
            "bolt_drill.s", 
            *base_normal, *base_hole_centers->vertices[i], 
            *top_normal, *top_hole_centers->vertices[i], 
            k_hole_diameter
        );
    }

    char *drill_set = bzy_make_union(
        "drills.c", drills, false
    );

    char *riser = BREEZY_MAKE_DIFFERENCE(basename, true, body, drill_set);

    free_arena(&scratch);

    return riser;
}

int main (int argc, char *argv[]) {
    assert(argc == 2);
    bzy_open_db(argv[1], "risers");

    make_riser(10, -0.5, "riser.r");

    bzy_close_db();
}
