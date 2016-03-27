//
//  mod_dummy.c
//  SymUniverse - Dummy module that does nothing. Can be used as a template.
//
//  Created by J. Lowell Wofford on 3/25/16.
//  Copyright © 2016 J. Lowell Wofford. All rights reserved.
//

// -- MODULE RULES (these MUST be followed for consistency) --
// 1. Do not modify ps.  It's there for consistency only.
// 2. Do not define cfg globally.  Different pipeline instances keep track of different cfgs.
// 3. Never shrink p.  If you want to delete a particle, mark it with PARTICLE_FLAG_DELETE and return MOD_RET_PACK.  Packing will take care of it later.
// 4. If you want to add a particle, append it to the end and mark it with PARTICLE_FLAG_CREATE. (Or use helper int slice_append_particle(Slice *s, Particle *p))
// 5. If you open or allocate anything globally, clean it up in the finalizer.

// -- MODULE GUIDELINES (these should probably be followed) --
// 1. help() should give a full description of the module, module options, and where it should be placed in the pipeline.
// 2. The standard format for cfg_str should be a comma delimated list of either binary options or option=value pairs, e.g. option1,option2=value2,etc.
// 3. Avoid global allocations if possible.
// 4. Only perform more than one transformation if doing so has significant optimization, otherwise write a seperate module.
// 5. You should probably check for and ignore any particle marked PARTICLE_FLAG_DELETE.
// 6. You can safely assume that any particle not marked PARTICLE_FLAG_CREATE has a one-to-one correspondence in ps, and visa verse.
// 7. If you detect a fundamental error (e.g. bad physics like superluminal particles), print an informative message and return MOD_RET_ABRT.
// 8. _Anything_ you print should be prefixed with "[module_name] ".  You can use the macro MPRINTF defined in sym.h.
// 9. Force modules should have an option whether or not to zero acceleration before computing.  They should default to no.
// 10. Major optimizations should have separate modules, e.g. a standard serial version, an OpenMP parallel version, a AMR version, etc.
// 11. There's nothing wrong with modules having their own modules, if it makes sense.
// 12. Use existing helper functions when dealing with Universes and their inhabitants.  This keeps code/data structures consistent. (see universe.h)
// 13. Report asymptotic algorithm performance in help(), e.g. O(NlogN), O(N^2), O(N), etc.

// -- FINAL NOTE: A module can do just about anything.  Be clear to the user what to expect.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sym.h"
#include "universe.h"
#include "SymUniverseConfig.h"

#define EXPORT __attribute__((visibility("default")))

#define _INTEGRATION_METH_PRE       0   // Use ps velocities to calculate displacement
#define _INTEGRATION_METH_LEAPFROG  1   // Use s velocites to calculate displacement (Symplectic evolution)

#define DEFAULT_TIMESTEP            1.0 // Time interval per slice
#define DEFAULT_BOUNDARY_METH       _boundary_periodic
#define DEFAULT_INTEGRATION_METH    _integrate_leapfrog

#define _NOPT           3
#define _OPT_BOUNDARY   0
#define _OPT_METHOD     1
#define _OPT_TIMESTEP   2

static const char *_opt_str[_NOPT] = { "boundary", "method", "timestamp" };

static const char *name_str = "integrate";      // Name _must_ be unique

typedef struct {
    int (*boundary_method)(Slice *s, Particle *p);
    int (*integration_method)(Particle *p, double ts);
    double timestep;
} Config;

int _boundary_periodic(Slice *s, Particle *p) {
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

int _boundary_elastic(Slice *s, Particle *p) {  // Biggest complication here is we have to consider particles that travel more than 1 boundary length beyond
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

int _boundary_diffuse(Slice *s, Particle *p) {
    if(p->pos.x > s->bound_max.x | p->pos.x < s->bound_min.x |
       p->pos.y > s->bound_max.y | p->pos.y < s->bound_min.y |
       p->pos.z > s->bound_max.z | p->pos.z < s->bound_min.z) {
        p->flags |= PARTICLE_FLAG_DELETE;
        return MOD_RET_PACK;
    }
    return MOD_RET_OK;
}

int _integrate_pre(Particle *p, double ts) {
    p->pos.x += p->vel.x * ts;
    p->pos.y += p->vel.y * ts;
    p->pos.z += p->vel.z * ts;
    p->vel.x += p->acc.x * ts;
    p->vel.y += p->acc.y * ts;
    p->vel.z += p->acc.z * ts;
    return MOD_RET_OK;
}

int _integrate_leapfrog(Particle *p, double ts) {
    p->vel.x += p->acc.x * ts;
    p->vel.y += p->acc.y * ts;
    p->vel.z += p->acc.z * ts;
    p->pos.x += p->vel.x * ts;
    p->pos.y += p->vel.y * ts;
    p->pos.z += p->vel.z * ts;
    return MOD_RET_OK;
}

int _get_opt_idx(const char *opt_str) {
    for(int i = 0; i < _NOPT; i++) {
        if(strcmp(opt_str, _opt_str[i]) == 0) { return i; }
    }
    return -1;
}

__attribute__((constructor))
static void initializer(void) {             // Called when module is opened (dlopen)
    
}

__attribute__((destructor))
static void finalizer(void) {               // Called when module is closed (dlclose)
    
}

EXPORT
const char *name(void) {                    // Required for name reporting.
    return name_str;
}

EXPORT
void *init(char *cfg_str) {           // Called when added to the pipeline.  Note: the pipeline can have multiple instances of a module with different cfg.
    Config *cfg = malloc(sizeof(Config));
    cfg->boundary_method = DEFAULT_BOUNDARY_METH;
    cfg->integration_method = DEFAULT_INTEGRATION_METH;
    cfg->timestep = DEFAULT_TIMESTEP;
    
    while(cfg_str[0] != '\0') {
        char *val = strsep(&cfg_str, ",");
        char *opt = strsep(&val, "=");
        switch(_get_opt_idx(opt)) {
            case _OPT_BOUNDARY:
                if(strcmp(val, "periodic") == 0) {
                    cfg->boundary_method = _boundary_periodic;
                } else if(strcmp(val, "elastic") == 0) {
                    cfg->boundary_method = _boundary_elastic;
                } else if(strcmp(val, "diffuse") == 0) {
                    cfg->boundary_method = _boundary_diffuse;
                } else {
                    MPRINTF("boundary must take one of the options: periodic, elastic or diffuse.\n", NULL);
                    free(cfg);
                    return NULL;
                }
                break;
            case _OPT_METHOD:
                if(strcmp(val, "pre") == 0) {
                    cfg->integration_method = _integrate_pre;
                } else if(strcmp(val, "leapfrog") == 0 ) {
                    cfg->integration_method = _integrate_leapfrog;
                } else {
                    MPRINTF("method must take one of the options: pre or leapfrog.\n", NULL);
                    free(cfg);
                    return NULL;
                }
                break;
            case _OPT_TIMESTEP:
                cfg->timestep = strtod(val, NULL);
                if(cfg->timestep <= 0) {
                    MPRINTF("timestep must be greater than zero.\n", NULL);
                    free(cfg);
                    return NULL;
                }
                break;
            default:
                MPRINTF("Invalid argument, %s.  Valid options are: boundary=?,method=?,timestep=?\n", opt);
                free(cfg);
                return NULL;
        }
    }
    
    return (void *)cfg;                     // This void pointer refers to internal configuration.  Can be anything useful.
}

EXPORT
void deinit(Config *cfg) {                  // Called when pipeline is deconstructed.
    free(cfg);
}

EXPORT
void help(void) {
    MPRINTF("This module does timestep integration.\n", NULL);
    MPRINTF("This is a simple algorithm with O(N) asymptotic performance.\n", NULL);
    MPRINTF("Integration should typically happen after forces and before collision detection.\n", NULL);
    MPRINTF("Initialization parameters take the form: option1=value1,option2=value2,...\n", NULL);
    MPRINTF("Available options are:\n", NULL);
    MPRINTF("\t- boundary: boundary conditions (default: periodic)\n", NULL);
    MPRINTF("\t\t- periodic (particles pass from one side to the other, i.e. Asteroids (tm) style).\n", NULL);
    MPRINTF("\t\t- elastic (particles bounce off of walls elastically).\n", NULL);
    MPRINTF("\t\t- diffuse (particles escape the system and disappear).\n", NULL);
    MPRINTF("\t- method: integration method (default: leapfrog)\n", NULL);
    MPRINTF("\t\t- pre (particles move based on velocities in the previous slice, then velocities are adjusted.\n", NULL);
    MPRINTF("\t\t- leapfrog (particle velocities are adjusted, then positions are adjusted accordingly.  This preseverse symplectic evolution.\n", NULL);
    MPRINTF("\t- timestep: takes any double value greater than zero.  This is the time period between each slice. (default: 1.0)\n", NULL);
    MPRINTF("Example: -m integrate[boundary=periodic,method=leapfrog,timestep=0.001]\n", NULL);
}

EXPORT
int exec(Config *cfg, Slice *ps, Slice *s) {  // Main execution loop.  Maps (ps, s) -> s.  Should _not_ modify ps.  Uses cfg to specify pipeline params.
    int ret = MOD_RET_OK;
    for(int i = 0; i < s->nbody; i++) {
        if(s->bodies[i].flags & PARTICLE_FLAG_DELETE) { continue; }
        ret |= cfg->integration_method(&s->bodies[i], cfg->timestep);
        ret |= cfg->boundary_method(s, &s->bodies[i]);
    }
    return ret;                             // Return value can control flow of overall execution, see MOD_RET_*
}
