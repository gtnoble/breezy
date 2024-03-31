#/usr/bin/rc

gcc \
    src/arena.c src/db.c src/math.c test/test.c \
    -fsanitize=address \
    -fsanitize=leak \
    -fsanitize=undefined \
    -Wall \
    -Wextra \
    -Wpedantic \
    -Werror \
    -g \
    -Iinclude \
    -I/home/gtnoble/Software/BRL-CAD/brlcad/include/brlcad \
    -I/home/gtnoble/Software/BRL-CAD/brlcad/include \
    -L/home/gtnoble/Software/BRL-CAD/brlcad/lib \
    -o wdb_test \
    -lm -lwdb -lrt

gcc \
    src/arena.c src/db.c src/math.c test/skateboard_risers.c \
    -fsanitize=address \
    -fsanitize=leak \
    -fsanitize=undefined \
    -Wall \
    -Wextra \
    -Wpedantic \
    -Werror \
    -g \
    -Iinclude \
    -I/home/gtnoble/Software/BRL-CAD/brlcad/include/brlcad \
    -I/home/gtnoble/Software/BRL-CAD/brlcad/include \
    -L/home/gtnoble/Software/BRL-CAD/brlcad/lib \
    -o skateboard_risers \
    -lm -lwdb -lrt