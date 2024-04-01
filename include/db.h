#include <stdbool.h>

#include <stddef.h>
#include <wdb.h>

#include "arena.h"
#include "bzy_math.h"

typedef struct ObjectName {
    char *basename;
    int object_count;
    struct ObjectName *next_object;
} BzyObject;

typedef struct BzyDb {
    struct rt_wdb *db_file;
    BzyObject *first_object;
    BzyObject *latest_object;
    Arena arena;
} BzyDb;

#define BZY_MAKE_UNION(basename, is_region, ...) \
    bzy_make_union(basename, (const char *[]) {__VA_ARGS__, NULL}, is_region)
#define BZY_MAKE_DIFFERENCE(basename, is_region, ...) \
    bzy_make_difference(basename, (const char *[]) {__VA_ARGS__, NULL}, is_region)

void bzy_open_db(const char *filepath, const char *database_title);
void bzy_close_db(void);

typedef struct BzyStrings {
    size_t num_elements;
    char **strings;
} BzyStrings;

BzyStrings *bzy_new_strings(size_t num_strings, Arena *arena);
BzyStrings bzy_make_strings(char *strings[], size_t num_strings);
char **bzy_strings_get(size_t index, BzyStrings *strings);

#define BZY_MAKE_STRINGS(...) \
    bzy_make_strings( \
        (char *[]) {__VA_ARGS__}, \
        sizeof((char *[]) {__VA_ARGS__}) / sizeof(char *) \
    )

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
    BzyStrings names, 
    db_op_t operation, 
    bool is_region
);
char *bzy_make_union(const char *basename, BzyStrings names, bool is_region);
char *bzy_make_difference(const char *basename, BzyStrings names, bool is_region);
char *bzy_make_drill(
    const char *basename,
    const Vector3D start_normal,
    const Vector3D start_center,
    const Vector3D end_normal,
    const Vector3D end_center,
    double diameter
);
