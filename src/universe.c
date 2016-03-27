//
//  universe.c
//  SymUniverse - Universe utility functions
//
//  Created by J. Lowell Wofford on 3/25/16.
//  Copyright © 2016 J. Lowell Wofford. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <string.h>
#include "universe.h"
#include "SymUniverseConfig.h"

void vector_add(Vector *dst, Vector *a, Vector *b) {
    dst->x = a->x + b->x;
    dst->y = a->y + b->y;
    dst->z = a->z + b->z;
}

void vector_sub(Vector *dst, Vector *a, Vector *b) {
    dst->x = a->x - b->x;
    dst->y = a->y - b->y;
    dst->z = a->z - b->z;
}

double vector_dot(Vector *a, Vector *b) {
    return a->x * b->x + a->y * b->y + a->z * b->z;
}

int slice_free(Slice *s) {
    free(s->bodies);
    free(s);
    return 0;
}

Slice *slice_copy(Slice *s) {
    Slice *new = malloc(sizeof(Slice));
    memcpy(new, s, sizeof(Slice));
    new->bodies = malloc(sizeof(Particle) * new->nbody);
    memcpy(new->bodies, s->bodies, sizeof(Particle) * new->nbody);
    return new;
}

void slice_pack(Slice *s) {  // Repack particles (e.g. if some have been marked to delete)
    // TODO
}

Universe *universe_create(const char *path) {
    Universe *u = malloc(sizeof(Universe));
    u->path = path;
    u->fstream = fopen(path, "w+");
    u->is_open = 1;
    u->is_modified = 0;
    u->nslice = 0;
    u->slice_idx = calloc(1, sizeof(long));
    fwrite(&u->nslice, sizeof(uint64_t), 1, u->fstream);
    
    return u;
}

Universe *universe_open(const char *path) {
    Universe *u = malloc(sizeof(Universe));
    u->path = path;
    u->fstream = fopen(path, "r+");
    u->is_open = 1;
    u->is_modified = 0;
    fread(&u->nslice, sizeof(uint64_t), 1, u->fstream);
    
    u->slice_idx = malloc(sizeof(long)*u->nslice);
    fseek(u->fstream, -(sizeof(long)*u->nslice), SEEK_END);
    fread(u->slice_idx, sizeof(long), u->nslice, u->fstream);
    
    return u;
}

int universe_close(Universe *u) {
    fclose(u->fstream);
    u->is_open = 0;
    return 0;
}

int universe_free(Universe *u) {
    if(u->is_open == 1) {
        return -1;
    }
    free(u->slice_idx);
    free(u);
    return 0;
}

Slice *universe_get_slice(Universe *u, uint64_t slice) {
    Slice *s = malloc(sizeof(Slice));
    fseek(u->fstream, u->slice_idx[slice], SEEK_SET);
    fread(&s->time, sizeof(uint64_t) + sizeof(Vector), 2, u->fstream); // reads time, nbody and the two boundary vectors
    
    s->bodies = malloc(sizeof(Particle)*s->nbody);
    fread(s->bodies, sizeof(Particle), s->nbody, u->fstream);
    
    return s;
}

Slice *universe_get_first_slice(Universe *u) {
    return universe_get_slice(u, 0);
}

Slice *universe_get_last_slice(Universe *u) {
    return universe_get_slice(u, u->nslice - 1);
}

int universe_append_slice(Universe *u, Slice *s) {
    ++u->nslice;
    rewind(u->fstream);
    fwrite(&u->nslice, sizeof(uint64_t), 1, u->fstream);
    
    u->slice_idx = realloc(u->slice_idx, sizeof(long)*u->nslice);
    fseek(u->fstream, -(sizeof(long)*(u->nslice - 1)), SEEK_END);
    u->slice_idx[u->nslice - 1] = ftell(u->fstream);
    
    fwrite(&s->time, sizeof(uint64_t) + sizeof(Vector), 2, u->fstream); // writes time, nbody and the two boundary vectors
    fwrite(s->bodies, sizeof(Particle), s->nbody, u->fstream);
    fwrite(u->slice_idx, sizeof(long), u->nslice, u->fstream);
    
    return 0;
}
