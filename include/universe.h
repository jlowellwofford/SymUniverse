//
//  universe.h
//  SymUniverse - structs/types of universes + utility functions
//
//  Created by J. Lowell Wofford on 3/25/16.
//  Copyright © 2016 J. Lowell Wofford. All rights reserved.
//

// Universe file format: (while a little complicated, this allows for universes with different size slices)
// Header:
//  - uint64_t nslice
// Slices:
//  - Slice
//  -- uint64_t time
//  -- uint64_t nbody
//  -- double bound_min[3]
//  -- double bound_max[3]
//  -- Particles:
//  --- uint64_t flags
//  --- double mass
//  --- double radius
//  --- double pos[3]
//  --- double vel[3]
//  --- double acc[3]
// Index:
//  long slice_pos[nslice]

#ifndef universe_h
#define universe_h

#include <stdint.h>

#define PARTICLE_FLAG_DELETE 1      // Indicates a particle is to be deleted
#define PARTICLE_FLAG_CREATE 2      // Keeps track of particles that weren't part of the original universe
#define PARTICLE_FLAG_NOCOLL 4      // This particle doesn't collide

#pragma pack(4) // Note: on most architectures pack(4) does nothing.  This is just to be certain.
typedef struct Vector {     // Simple 3-vector
    double x;
    double y;
    double z;
} Vector;

#pragma pack(4)
typedef struct Particle {   // Particle properties
    uint64_t    flags;      // _way_ more flag space than we need, but keeps packing nice.
    double      mass;
    double      charge;
    double      radius;
    Vector      pos;
    Vector      vel;
    Vector      acc;
} Particle;

#pragma pack(4)
typedef struct Slice {      // Time slice: time index + body count + body array
    uint64_t    time;
    uint64_t    nbody;
    Vector      bound_min;
    Vector      bound_max;
    Particle    *bodies;
} Slice;

typedef struct Universe {   // A universe: number of slices + slice array
    const char  *path;
    FILE        *fstream;
    char        is_open;
    char        is_modified;
    uint64_t    nslice;
    long        *slice_idx;
} Universe;

// Function definitions
int slice_free(Slice *s);
Slice *slice_copy(Slice *s);
void slice_pack(Slice *);
int slice_append_particle(Slice *s, Particle *p);

Universe *universe_create(const char *path);
Universe *universe_open(const char *path);
int universe_close(Universe *u);
int universe_free(Universe *u);
Slice *universe_get_slice(Universe *u, uint64_t slice);
Slice *universe_get_first_slice(Universe *u);
Slice *universe_get_last_slice(Universe *u);

int universe_append_slice(Universe *u, Slice *s );

#endif /* universe_h */
