//
//  boundary.c
//  SymUniverse - This module enforces boundary conditions.  Use this if you disabled boundaries in integrate.
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
#include "boundaries.h"
#include "SymUniverseConfig.h"

#define EXPORT __attribute__((visibility("default")))

#define DEFAULT_BOUNDARY_METH       boundary_periodic

#define _NOPT           1
#define _OPT_BOUNDARY   0

static const char *_opt_str[_NOPT] = { "boundary" };

EXPORT
const char *name = "boundary";      // Name _must_ be unique

typedef struct {
    int (*boundary_method)(Slice *s, Particle *p);
} Config;

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
void *init(char *cfg_str) {           // Called when added to the pipeline.  Note: the pipeline can have multiple instances of a module with different cfg.
    Config *cfg = malloc(sizeof(Config));
    cfg->boundary_method = DEFAULT_BOUNDARY_METH;
    
    while(cfg_str != NULL && cfg_str[0] != '\0') {
        char *val = strsep(&cfg_str, ",");
        char *opt = strsep(&val, "=");
        switch(_get_opt_idx(opt)) {
            case _OPT_BOUNDARY:
                if(strcmp(val, "periodic") == 0) {
                    cfg->boundary_method = boundary_periodic;
                } else if(strcmp(val, "elastic") == 0) {
                    cfg->boundary_method = boundary_elastic;
                } else if(strcmp(val, "diffuse") == 0) {
                    cfg->boundary_method = boundary_diffuse;
                } else if(strcmp(val, "none") == 0) {
                    MPRINTF("Warning: you have chosen not to use boundary conditions. Make sure this is handled by another module!\n", NULL);
                    cfg->boundary_method = boundary_none;
                } else {
                    MPRINTF("boundary must take one of the options: periodic, elastic, diffuse or none.\n", NULL);
                    free(cfg);
                    return NULL;
                }
                break;
            default:
                MPRINTF("Invalid argument, %s.  Valid options are: boundary=?\n", opt);
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
    MPRINTF("This module enforces boundary conditions.\n", NULL);
    MPRINTF("This is a simple algorithm with O(N) asymptotic performance.\n", NULL);
    MPRINTF("This will often be near the end of your pipeline, and should happen after collision detection.\n", NULL);
    MPRINTF("Initialization parameters take the form: option1=value1,option2=value2,...\n", NULL);
    MPRINTF("Available options are:\n", NULL);
    MPRINTF("\t- boundary: boundary conditions (default: periodic)\n", NULL);
    MPRINTF("\t\t- periodic (particles pass from one side to the other, i.e. Asteroids (tm) style).\n", NULL);
    MPRINTF("\t\t- elastic (particles bounce off of walls elastically).\n", NULL);
    MPRINTF("\t\t- diffuse (particles escape the system and disappear).\n", NULL);
    MPRINTF("\t\t- none (no boundaries enforced.  Use this if you are going to use collision detection.", NULL);
    MPRINTF("\t\t\tWarning: you need to do boundary enforcement at some point.\n",NULL);
    MPRINTF("\t\t\tIf you've disabled it here, make sure another module does it!\n", NULL);
    MPRINTF("Example: -m boundary[boundary=periodic]\n", NULL);
}

EXPORT
int exec(Config *cfg, Slice *ps, Slice *s) {  // Main execution loop.  Maps (ps, s) -> s.  Should _not_ modify ps.  Uses cfg to specify pipeline params.
    int ret = MOD_RET_OK;
    for(int i = 0; i < s->nbody; i++) {
        if(s->bodies[i].flags & PARTICLE_FLAG_DELETE) { continue; }
        ret |= cfg->boundary_method(s, &s->bodies[i]);
    }
    return ret;                             // Return value can control flow of overall execution, see MOD_RET_*
}
