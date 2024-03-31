#include <assert.h>

#include "db.h"
#include "bzy_math.h"

int main(int argc, char *argv[]) {

    assert(argc == 2);
    char *test_filename = argv[1];

    bzy_open_db(test_filename, "test_db");

    BREEZY_MAKE_UNION(
        "combo.r", 
        true, 
        bzy_make_rpp(
            "rpp.s", 
            (const Vector3D) {0.0, 0.0, 0.0}, 
            (const Vector3D) {2.0, 4.0, 2.5}
        ),
        bzy_make_sphere(
            "sph.s", 
            (const Vector3D) {1.0, 2.0, 3.0}, 
            0.75
        )
    );

    bzy_close_db();
}
