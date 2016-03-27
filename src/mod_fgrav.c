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
// 8. _Anything_ you print should be prefixed with "[module_name] ".
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

#define DEFAULT_CLEARA 0

EXPORT
const char *name = "fgrav";      // Name _must_ be unique

typedef struct {
    int cleara;
} Config;

__attribute__((constructor))
static void initializer(void) {             // Called when module is opened (dlopen)
    
}

__attribute__((destructor))
static void finalizer(void) {               // Called when module is closed (dlclose)
    
}

EXPORT
void *init(char *cfg_str) {                 // Called when added to the pipeline.  Note: the pipeline can have multiple instances of a module with different cfg.
    
    Config *cfg = malloc(sizeof(Config));
    cfg->cleara = DEFAULT_CLEARA;

    while(cfg_str != NULL && cfg_str[0] != '\0') {
        char *val = strsep(&cfg_str, ",");
        char *opt = strsep(&val, "=");
        if(strcmp(opt, "cleara") == 0) {
            cfg->cleara = atoi(val);
            if(cfg->cleara != 0 && cfg->cleara != 1) {
                MPRINTF("Option cleara accepts only 0 (disable) or 1 (enable)!\n", NULL);
                free(cfg);
                return NULL;
            }
        } else {
            MPRINTF("Option not recognized! See help (-h) for options.\n", NULL);
            free(cfg);
            return NULL;
        }
    }
    
    return (void *)cfg;                 // This void pointer refers to internal configuration.  Can be anything useful.
}

EXPORT
void deinit(Config *cfg) {                    // Called when pipeline is deconstructed.
    free(cfg);
}

EXPORT
void help(void) {
    MPRINTF("This module calculates gravitational acceleration.\n", NULL);
    MPRINTF("This is a simplistic algorithm with asymptotic performance of O(NlogN).\n", NULL);
    MPRINTF("Usually, force modules should come first in the pipeline, followed by integration and collision detection.\n", NULL);
    MPRINTF("There is one available option:\n", NULL);
    MPRINTF("\t- cleara: reset accelerations to zero before calculating?\n", NULL);
    MPRINTF("\t\tTakes two options: 0 to disable, 1 to enable.\n", NULL);
    MPRINTF("The first force module in the pipeline should set cleara=1.\n", NULL);
    MPRINTF("Example: -m fgrav[cleara=1]\n", NULL);
}

EXPORT
int exec(Config *cfg, Slice *ps, Slice *s) {  // Main execution loop.  Maps (ps, s) -> s.  Should _not_ modify ps.  Uses cfg to specify pipeline params.
    // Slightly less efficient to do this separately, but makes the code reusable later
    // (e.g. for a parallel version)
    if(cfg->cleara) {
        for(int i = 0; i < s->nbody; i++) {
            s->bodies[i].acc.x = 0;
            s->bodies[i].acc.y = 0;
            s->bodies[i].acc.z = 0;
        }
    }
    
    // The calculation loop
    for(int i = 0; i < s->nbody; i++) {
        for(int j = i + 1; j < s->nbody; j++) {
            Vector r;
            vector_sub(&r, &s->bodies[i].pos, &s->bodies[j].pos);
            double a = s->bodies[i].mass * s->bodies[j].mass / pow(vector_dot(&r, &r),1.5); // note: this pow() takes about 72% of total compute time
            s->bodies[i].acc.x += a * r.x;
            s->bodies[i].acc.y += a * r.y;
            s->bodies[i].acc.z += a * r.z;
            s->bodies[j].acc.x -= a * r.x;
            s->bodies[j].acc.y -= a * r.y;
            s->bodies[j].acc.z -= a * r.z;
        }
    }
    
    return MOD_RET_OK;                      // Return value can control flow of overall execution, see MOD_RET_*
}
