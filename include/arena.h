
#ifndef BREEZY_ARENA
#define BREEZY_ARENA

#include <stdlib.h>

typedef struct {
    char *memory;
    size_t offset;
    size_t size;
} Arena;

Arena make_arena(size_t max_size);
void free_arena(Arena *arena);
void *arena_allocate(size_t size, Arena *arena);
void *arena_allocate_aligned(size_t size, size_t alignment, Arena *arena);
void arena_reset(Arena *arena);
size_t arena_checkpoint(const Arena *arena);
void arena_restore(size_t checkpoint, Arena *arena);
void *arena_allocate_zeros(size_t size, Arena *arena);
char *arena_duplicate_string(const char *string, Arena *arena);
char *arena_sprintf(Arena *arena, const char *format, ...);

#endif
