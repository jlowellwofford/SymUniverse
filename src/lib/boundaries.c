//
//  boundaries.c
//  SymUniverse - Routines for resolving boundary conditions.
//
//  Created by J. Lowell Wofford on 3/25/16.
//  Copyright © 2016 J. Lowell Wofford. All rights reserved.
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "boundaries.h"
#include "universe.h"
#include "sym.h"
#include "SymUniverseConfig.h"

int boundary_periodic(Slice *s, Particle *p) {
    if(p->pos.x > s->bound_max.x) {
        double bounces;
        double off = modf((p->pos.x - s->bound_min.x) / (s->bound_max.x - s->bound_min.x), &bounces) * (s->bound_max.x - s->bound_min.x);
        p->pos.x = s->bound_min.x + off;
    } else if(p->pos.x < s->bound_min.x) {
        double bounces;
        double off = modf(-(p->pos.x - s->bound_min.x) / (s->bound_max.x - s->bound_min.x), &bounces) * (s->bound_max.x - s->bound_min.x);
        p->pos.x = s->bound_max.x - off;
    }
    
    if(p->pos.y > s->bound_max.y) {
        double bounces;
        double off = modf((p->pos.y - s->bound_min.y) / (s->bound_max.y - s->bound_min.y), &bounces) * (s->bound_max.y - s->bound_min.y);
        p->pos.y = s->bound_min.y + off;
    } else if(p->pos.y < s->bound_min.y) {
        double bounces;
        double off = modf(-(p->pos.y - s->bound_min.y) / (s->bound_max.y - s->bound_min.y), &bounces) * (s->bound_max.y - s->bound_min.y);
        p->pos.y = s->bound_max.y - off;
    }
    
    if(p->pos.z > s->bound_max.z) {
        double bounces;
        double off = modf((p->pos.z - s->bound_min.z) / (s->bound_max.z - s->bound_min.z), &bounces) * (s->bound_max.z - s->bound_min.z);
        p->pos.z = s->bound_min.z + off;
    } else if(p->pos.z < s->bound_min.z) {
        double bounces;
        double off = modf(-(p->pos.z - s->bound_min.z) / (s->bound_max.z - s->bound_min.z), &bounces) * (s->bound_max.z - s->bound_min.z);
        p->pos.z = s->bound_max.z - off;
    }
    return MOD_RET_OK;
}

int boundary_elastic(Slice *s, Particle *p) {  // Biggest complication here is we have to consider particles that travel more than 1 boundary length beyond
    if(p->pos.x > s->bound_max.x) {
        double bounces;
        double off = modf((p->pos.x - s->bound_min.x) / (s->bound_max.x - s->bound_min.x), &bounces) * (s->bound_max.x - s->bound_min.x);
        p->pos.x = ((int)bounces % 2) ? s->bound_min.x + off : s->bound_max.x - off;
    } else if(p->pos.x < s->bound_min.x) {
        double bounces;
        double off = modf(-(p->pos.x - s->bound_min.x) / (s->bound_max.x - s->bound_min.x), &bounces) * (s->bound_max.x - s->bound_min.x);
        p->pos.x = ((int)bounces % 2) ? s->bound_max.x - off : s->bound_min.x + off;
    }
    
    if(p->pos.y > s->bound_max.y) {
        double bounces;
        double off = modf((p->pos.y - s->bound_min.y) / (s->bound_max.y - s->bound_min.y), &bounces) * (s->bound_max.y - s->bound_min.y);
        p->pos.y = ((int)bounces % 2) ? s->bound_min.y + off : s->bound_max.y - off;
    } else if(p->pos.y < s->bound_min.y) {
        double bounces;
        double off = modf(-(p->pos.y - s->bound_min.y) / (s->bound_max.y - s->bound_min.y), &bounces) * (s->bound_max.y - s->bound_min.y);
        p->pos.y = ((int)bounces % 2) ? s->bound_max.y - off : s->bound_min.y + off;
    }
    
    if(p->pos.z > s->bound_max.z) {
        double bounces;
        double off = modf((p->pos.z - s->bound_min.z) / (s->bound_max.z - s->bound_min.z), &bounces) * (s->bound_max.z - s->bound_min.z);
        p->pos.z = ((int)bounces % 2) ? s->bound_min.z + off : s->bound_max.z - off;
    } else if(p->pos.z < s->bound_min.z) {
        double bounces;
        double off = modf(-(p->pos.z - s->bound_min.z) / (s->bound_max.z - s->bound_min.z), &bounces) * (s->bound_max.z - s->bound_min.z);
        p->pos.z = ((int)bounces % 2) ? s->bound_max.z - off : s->bound_min.z + off;
    }
    return MOD_RET_OK;
}

int boundary_diffuse(Slice *s, Particle *p) {
    if(p->pos.x > s->bound_max.x | p->pos.x < s->bound_min.x |
       p->pos.y > s->bound_max.y | p->pos.y < s->bound_min.y |
       p->pos.z > s->bound_max.z | p->pos.z < s->bound_min.z) {
        p->flags |= PARTICLE_FLAG_DELETE;
        return MOD_RET_PACK;
    }
    return MOD_RET_OK;
}

int boundary_none(Slice *s, Particle *p) {
    return MOD_RET_OK;
}