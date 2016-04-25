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
#include <errno.h>
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

void vector_cross(Vector *dst, Vector *a, Vector *b) {
    dst->x = a->y * b->z - a->z * b->y;
    dst->y = a->z * b->x - a->x * b->z;
    dst->z = a->x * b->y - a->y * b->x;
}

int vector_equal(Vector *a, Vector *b) {
    if(a->x == b->x && a->y == b->y && a->z == b->z) { return 1; }
    return 0;
}

int slice_free(Slice *s) {
    free(s->bodies);
    free(s);
    return 0;
}

Slice *slice_copy(Slice *s) {
    Slice *new = malloc(sizeof(Slice));
    if(new == NULL) {
        return NULL;
    }
    memcpy(new, s, sizeof(Slice));
    new->bodies = malloc(sizeof(Particle) * new->nbody);
    if(new->bodies == NULL) {
        free(new);
        return NULL;
    }
    memcpy(new->bodies, s->bodies, sizeof(Particle) * new->nbody);
    return new;
}

int slice_pack(Slice *s) {  // Repack particles (e.g. if some have been marked to delete)
    Particle *new = malloc(sizeof(Particle) * s->nbody);
    if(new == NULL) {
        printf("Memory allocation error.\n");
        return 0;
    }
    uint64_t nnew = 0;
    
    for(int i = 0; i < s->nbody; i++) {
        s->bodies[i].flags &= ~PARTICLE_FLAG_CREATE;   // Remove CREATE flag.  It's actually faster to do this every time than check and remove.
        if(s->bodies[i].flags & PARTICLE_FLAG_DELETE) { continue; }
        memcpy(&new[nnew], &s->bodies[i], sizeof(Particle));
        nnew++;
    }
    realloc(new, sizeof(Particle) * nnew);
    if(new == NULL) {
        printf("Memory allocation error\n");
        return 0;
    }
    s->nbody = nnew;
    free(s->bodies);
    s->bodies = new;
    return 1;
}

void slice_clear_create(Slice *s) {
    for(int i = 0; i < s->nbody; i++) {
        s->bodies[i].flags &= ~PARTICLE_FLAG_CREATE;   // Remove CREATE flag.  It's actually faster to do this every time than check and remove.
    }
}

int slice_append_particle(Slice *s, Particle *p) {
    ++s->nbody;
    s->bodies = realloc(s->bodies, sizeof(Particle) * s->nbody);
    memcpy(&s->bodies[s->nbody-1], p, sizeof(Particle));
    return 0;
}

Universe *universe_create(const char *path) {
    UniverseHeader header;
    Universe *u = malloc(sizeof(Universe));
    if(u == NULL) {
        printf("Memory allocation error.\n");
        return NULL;
    }
    u->path = path;
    if(!access(path, W_OK)) {
        printf("Cannot create universe file: file already exists!\n");
        free(u);
        return NULL;
    }
    u->fstream = fopen(path, "w+");
    if(u->fstream == NULL) {
        printf("Could not open universe file, %s, because: %s\n", path, strerror(errno));
        free(u);
        return NULL;
    }
    u->is_open = 1;
    u->is_modified = 0;
    u->nslice = 0;
    u->slice_idx = calloc(1, sizeof(long));
    if(u->slice_idx == NULL) {
        printf("Memory allocation error.\n");
        free(u);
        return NULL;
    }
    strncpy(header.string, UNIVERSE_STRING, sizeof(header.string));
    header.version = UNIVERSE_VERSION;
    header.nslice = 0;
    fwrite(&header, sizeof(UniverseHeader), 1, u->fstream);
    
    return u;
}

Universe *universe_open(const char *path) {
    UniverseHeader header;
    Universe *u = malloc(sizeof(Universe));
    if(u == NULL) {
        printf("Memory allocation error.\n");
        return NULL;
    }
    u->path = path;
    u->fstream = fopen(path, "r+");
    if(u->fstream == NULL) {
        printf("Could not open universe file, %s, because: %s\n", path, strerror(errno));
        free(u);
        return NULL;
    }
    u->is_open = 1;
    u->is_modified = 0;
    fread(&header, sizeof(UniverseHeader), 1, u->fstream);
    if(strncmp(header.string, UNIVERSE_STRING, sizeof(header.string)) != 0) {
        printf("%s does not appear to be a valid Universe Data File!\n", path);
        free(u);
        return NULL;
    } else if(header.version != UNIVERSE_VERSION) {
        printf("Universe Data File version mismatch.  File is %d, we need %d.\n", header.version, UNIVERSE_VERSION);
        free(u);
        return NULL;
    }
    u->nslice = header.nslice;
    
    u->slice_idx = malloc(sizeof(long)*u->nslice);
    if(u->slice_idx == NULL) {
        printf("Memory allocation error.\n");
        fclose(u->fstream);
        free(u);
        return NULL;
    }
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
    if(s == NULL) {
        printf("Memory allocation error.\n");
        return NULL;
    }
    fseek(u->fstream, u->slice_idx[slice], SEEK_SET);
    fread(&s->time, sizeof(uint64_t) + sizeof(Vector), 2, u->fstream); // reads time, nbody and the two boundary vectors
    
    s->bodies = malloc(sizeof(Particle)*s->nbody);
    if(s->bodies == NULL) {
        printf("Memory allocation error.\n");
        free(s);
        return NULL;
    }
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
    
    UniverseHeader header;
    strncpy(header.string, UNIVERSE_STRING, sizeof(header.string));
    header.version = UNIVERSE_VERSION;
    header.nslice = u->nslice;
    
    rewind(u->fstream);
    fwrite(&header, sizeof(UniverseHeader), 1, u->fstream);
    
    u->slice_idx = realloc(u->slice_idx, sizeof(long)*u->nslice);
    if(u->slice_idx == NULL) {
        printf("Memory allocation error.\n");
        return 0;
    }
    fseek(u->fstream, -(sizeof(long)*(u->nslice - 1)), SEEK_END);
    u->slice_idx[u->nslice - 1] = ftell(u->fstream);
    
    fwrite(&s->time, sizeof(uint64_t) + sizeof(Vector), 2, u->fstream); // writes time, nbody and the two boundary vectors
    fwrite(s->bodies, sizeof(Particle), s->nbody, u->fstream);
    fwrite(u->slice_idx, sizeof(long), u->nslice, u->fstream);
    
    return 1;
}
