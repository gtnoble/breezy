#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>

#include "arena.h"

Arena make_arena(size_t max_size) {
    Arena arena;
    arena.memory = mmap(
        NULL, 
        max_size, 
        PROT_READ | PROT_WRITE, 
        MAP_SHARED | MAP_ANONYMOUS, 
        -1, 
        0
    );
    arena.offset = 0;
    arena.size = max_size; 
    return arena;
}

void free_arena(Arena *arena) {
    munmap(arena->memory, arena->size);
}

void *arena_allocate(size_t size, Arena *arena) {
    if (arena->offset + size > arena->size - 1) {
        fprintf(
            stderr, 
            "error: Out of memory, cannot allocate from arena: "
            "arena size: %zu: "
            "memory remaining: %zu: "
            "requested: %zu", 
            arena->size,
            arena->size - arena->offset,
            size
        );
        exit(1);
    }
    void *allocated_memory = &arena->memory[arena->offset];
    arena->offset += size;
    return allocated_memory;
}

void *arena_allocate_aligned(size_t size, size_t alignment, Arena *arena) {
    arena->offset += alignment - ((arena->offset % alignment) % alignment);
    return arena_allocate(size, arena);
}

void *arena_allocate_zeros(size_t size, Arena *arena) {
    char *allocated_memory = &arena->memory[arena->offset];
    for (size_t i = 0; i < size; i++) {
        allocated_memory[arena->offset + i] = 0;
    }
    arena->offset += size;
    return (void *) allocated_memory;
}

char *arena_duplicate_string(const char *string, Arena *arena) {
    size_t length = strlen(string);
    char *duplicated_string = arena_allocate(sizeof(char) * (length + 1), arena);
    strcpy(duplicated_string, string);
    return duplicated_string;
}

char *arena_sprintf(Arena *arena, const char *format, ...) {

    char *printed_string = arena->memory + arena->offset;

    va_list argptr;
    va_start(argptr, format);
    int num_chars_printed = vsprintf(printed_string, format, argptr);
    assert(num_chars_printed >= 0);

    arena->offset += num_chars_printed + 1;

    return printed_string;
}

void arena_reset(Arena *arena) {
    arena->offset = 0;
}

size_t arena_checkpoint(const Arena *arena) {
    return arena->offset;
}

void arena_restore(size_t checkpoint, Arena *arena) {
    arena->offset = checkpoint;
}


