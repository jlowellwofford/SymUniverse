//
//  boundaries.h
//  SymUniverse - Routines for resolving boundary conditions.
//
//  Created by J. Lowell Wofford on 3/27/16.
//
//

#ifndef boundaries_h
#define boundaries_h

#include "universe.h"

typedef enum {
    periodic,
    elastic,
    diffuse,
    none
} BoundaryType;

int boundary_periodic(Slice *s, Particle *p);
int boundary_elastic(Slice *s, Particle *p);
int boundary_diffuse(Slice *s, Particle *p);
int boundary_none(Slice *s, Particle *p);

#endif /* boundaries_h */
