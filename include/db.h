#include <stdbool.h>

#include <wdb.h>

#include "arena.h"
#include "bzy_math.h"

typedef struct ObjectName {
    char *basename;
    int object_count;
    struct ObjectName *next_object;
} SlickObject;

typedef struct SlickDb {
    struct rt_wdb *db_file;
    SlickObject *first_object;
    SlickObject *latest_object;
    Arena arena;
} SlickDb;

#define BREEZY_MAKE_UNION(basename, is_region, ...) \
    bzy_make_union(basename, (const char *[]) {__VA_ARGS__, NULL}, is_region)
#define BREEZY_MAKE_DIFFERENCE(basename, is_region, ...) \
    bzy_make_difference(basename, (const char *[]) {__VA_ARGS__, NULL}, is_region)

void bzy_open_db(const char *filepath, const char *database_title);
void bzy_close_db(void);

char *bzy_make_sphere(const char *basename, const Vector3D vertex, double radius);
char *bzy_make_rpp(
    const char *basename, 
    const Vector3D min_point, 
    const Vector3D max_point
);
char *bzy_make_arb8(
    const char *basename,
    const Points *quad1,
    const Points *quad2
);
char *bzy_make_combination(
    const char *basename, 
    const char *names[], 
    db_op_t operation, 
    bool is_region
);
char *bzy_make_union(const char *basename, const char *names[], bool is_region);
char *bzy_make_difference(const char *basename, const char *names[], bool is_region);
char *bzy_make_drill(
    const char *basename,
    const Vector3D start_normal,
    const Vector3D start_center,
    const Vector3D end_normal,
    const Vector3D end_center,
    double diameter
);
