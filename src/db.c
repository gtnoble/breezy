#include <stdalign.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

#include "vmath.h"
#include "bu/app.h"
#include "bn.h"
#include "raytrace.h"
#include "rt/geom.h"

#include "db.h"
#include "arena.h"
#include "bzy_math.h"

static SlickDb db;

void bzy_open_db(const char *filepath, const char *database_title) {
    assert(db.db_file == NULL);
    assert(db.arena.memory == NULL);
    assert(db.arena.offset == 0);
    assert(db.first_object == NULL);
    assert(db.latest_object == NULL);

    db.arena = make_arena(10000000);
    db.db_file = wdb_fopen(filepath);

    mk_id_units(db.db_file, database_title, "m");
}

void bzy_close_db(void) {
    db_close(db.db_file->dbip);
    free_arena(&db.arena);
    db.first_object = NULL;
    db.latest_object = NULL;
    db.db_file = NULL;
}

SlickObject *search_base_name(const char *basename, SlickObject *root_object) {
    if (root_object == NULL) {
        return NULL;
    }
    else if (strcmp(basename, root_object->basename) == 0) {
        return root_object;
    }
    else {
        return search_base_name(basename, root_object->next_object);
    }
}

SlickObject *push_new_object(const char *basename) {
    SlickObject *new_object = arena_allocate_aligned(
        sizeof(SlickObject), alignof(SlickObject),&db.arena
    );
    new_object->basename = arena_duplicate_string(basename, &db.arena);
    new_object->object_count = 1;

    if (db.latest_object == NULL) {
        assert(db.first_object == NULL);
        db.first_object = new_object;
        db.latest_object = new_object;
    }
    else {
        db.latest_object->next_object = new_object;
    }
    return new_object;
}

char *evaluate_name(SlickObject object, Arena *arena) {
    return arena_sprintf(arena, "%s%d", object.basename, object.object_count);
}

char *new_object(const char *basename) {
    SlickObject *matched_object = search_base_name(basename, db.first_object);

    SlickObject *object;
    if (matched_object == NULL) {
        object = push_new_object(basename);
    }
    else {
        matched_object->object_count++;
        object = matched_object;
    }

    return evaluate_name(*object, &db.arena);
}

char *bzy_make_sphere(const char *basename, const Vector3D vertex, double radius) {
    char *name = new_object(basename);
    mk_sph(db.db_file, name, vertex, radius);
    return name;
}

char *bzy_make_rpp(
    const char *basename, 
    const Vector3D min_point, 
    const Vector3D max_point
) {
    char *name = new_object(basename);
    mk_rpp(
        db.db_file, 
        name, 
        min_point, 
        max_point
    );
    return name;
}

char *bzy_make_arb8(
    const char *basename,
    const Points *quad1,
    const Points *quad2
) {
    Arena scratch = make_arena(10000000);

    Points *poly1 = bzy_polygonalize(quad1, &scratch);
    Points *poly2 = bzy_polygonalize(quad2, &scratch);

    char *name = new_object(basename);
    mk_arb8(
        db.db_file, 
        name,
        bzy_flatten_vectors(
            (Vector3D * const []) {
                poly1->vertices[0],
                poly1->vertices[1],
                poly1->vertices[2],
                poly1->vertices[3],
                poly2->vertices[0],
                poly2->vertices[1],
                poly2->vertices[2],
                poly2->vertices[3]
            },
            8, 
            &scratch
        )
    );

    return name;
}

char *bzy_make_rcc(
    const char *basename,
    const Vector3D base_vertex,
    const Vector3D height,
    double radius
) {
    char *name = new_object(basename);
    
    mk_rcc(db.db_file, name, base_vertex, height, radius);
    return name;
}

char *bzy_make_drill(
    const char *basename,
    const Vector3D start_normal,
    const Vector3D start_center,
    const Vector3D end_normal,
    const Vector3D end_center,
    double diameter
) {
    double radius = diameter / 2.0;

    Arena scratch = make_arena(10000000);

    Vector3D *direction = bzy_vector_subn(end_center, start_center, &scratch);

    Vector3D *start_clearance = bzy_hole_clearance(*direction, start_normal, radius, &scratch);
    Vector3D *end_clearance = bzy_hole_clearance(*direction, end_normal, radius, &scratch);

    Vector3D *drill_start = bzy_vector_subn(start_center, *start_clearance, &scratch);
    Vector3D *drill_end = bzy_vector_addn(end_center, *end_clearance, &scratch);
    Vector3D *drill_direction = bzy_vector_subn(*drill_end, *drill_start, &scratch);

    char *drill = bzy_make_rcc(
        basename, 
        *drill_start, 
        *drill_direction, 
        radius
    );

    free_arena(&scratch);

    return drill;
}


//@todo make names a specified-length array
char *bzy_make_combination(
    const char *basename, 
    const char *names[], 
    db_op_t operation, 
    bool is_region
) {
    char *name = new_object(basename);

    struct wmember wm_hd;
    BU_LIST_INIT(&wm_hd.l);

    for (size_t i = 0; names[i] != NULL; i++) {
        mk_addmember(
            names[i],
            &wm_hd.l, 
            NULL, 
            (i == 0 ? WMOP_UNION : operation)
        );
    }

    mk_lcomb(
        db.db_file, 
        name, &wm_hd, 
        is_region, 
        NULL, 
        NULL, 
        NULL, 
        0
    );

    return name;
}

char *bzy_make_union(const char *basename, const char *names[], bool is_region) {
    return bzy_make_combination(basename, names, WMOP_UNION, is_region);
}

char *bzy_make_difference(const char *basename, const char *names[], bool is_region) {
    return bzy_make_combination(basename, names, WMOP_SUBTRACT, is_region);
}


