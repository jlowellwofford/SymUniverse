//
//  ptcollide.c
//  SymUniverse - Module for doing simple point-like collision detection and resolution. Does not guarantee conservation of angular momentum.
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

// FIXME: This module doesn't do anything yet!

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sym.h"
#include "universe.h"
#include "SymUniverseConfig.h"

#define EXPORT __attribute__((visibility("default")))

EXPORT
const char *name = "hscollide";      // Name _must_ be unique

typedef struct {                    // Not currently used, but there for later use
    int configed;
} Config;

__attribute__((constructor))
static void initializer(void) {             // Called when module is opened (dlopen)
    
}

__attribute__((destructor))
static void finalizer(void) {               // Called when module is closed (dlclose)
    
}

EXPORT
void *init(char *cfg_str) {                 // Called when added to the pipeline.  Note: the pipeline can have multiple instances of a module with different cfg.
    Config *cfg;
    if((cfg = malloc(sizeof(Config))) == NULL) {
        MPRINTF("Memory allocation failure.\n", NULL);
        return NULL;
    }
    cfg->configed = 1;
    return (void *)cfg;                     // This void pointer refers to internal configuration.  Can be anything useful.
}

EXPORT
void deinit(Config *cfg) {                    // Called when pipeline is deconstructed.
    free(cfg);
}

EXPORT
void help(void) {
    MPRINTF("This module resolves point-like collisions.\n", NULL);
    MPRINTF("Warning: this modules does NOT guarantee conservation of angular momentum.\n", NULL);
    MPRINTF("This simple algorithm is O(N^2), and takes no options.\n", NULL);
    MPRINTF("Typically this should be placed after forces and integration.\n", NULL);
    MPRINTF("For best results, disable boundary detection in integrate and use the boundary module after this.\n", NULL);
}

// !!!: It would be nice to have a recursive, time-ordered collision detection option.
// FIXME: This algorithm is not yet implemented!
EXPORT
int exec(Config *cfg, Slice *ps, Slice *s) {  // Main execution loop.  Maps (ps, s) -> s.  Should _not_ modify ps.  Uses cfg to specify pipeline params.
    for(int i = 0; i < ps->nbody; i++) {
        if(s->bodies[i].flags & (PARTICLE_FLAG_DELETE | PARTICLE_FLAG_CREATE)) { continue; }
        for(int j = i+1; j < ps->nbody; j++) {
            if(s->bodies[j].flags & (PARTICLE_FLAG_DELETE | PARTICLE_FLAG_CREATE)) { continue; }

        }
    }
    
    return MOD_RET_OK;                      // Return value can control flow of overall execution, see MOD_RET_*
}
