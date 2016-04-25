//
//  fgrav.c
//  SymUniverse - This module computes gravitational accelerations.
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

#define DEFAULT_G 1
#define DEFAULT_CLEARA 0
#define DEFAULT_PLUMMER2 0      // Used to impose a "Plummer sphere".  Shouldn't be needed if we're doing hscollide.

EXPORT
const char *name = "fgrav";      // Name _must_ be unique

typedef struct {
    int cleara;
    double plummer2;            // Plummer distance squared (we never use the un-squared version)
    double G;
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
        printf("Memory allocation failure.\n");
        return NULL;
    }
    cfg->cleara = DEFAULT_CLEARA;
    cfg->plummer2 = DEFAULT_PLUMMER2;
    cfg->G = DEFAULT_G;

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
        } else if(strcmp(opt, "plummer") == 0) {
            cfg->plummer2 = pow(strtod(val, NULL), 2);
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
    MPRINTF("This is a simplistic algorithm with asymptotic performance of O(N^2).\n", NULL);
    MPRINTF("Usually, force modules should come first in the pipeline, followed by integration and collision detection.\n", NULL);
    MPRINTF("There two available options:\n", NULL);
    MPRINTF("\t- cleara: reset accelerations to zero before calculating?\n", NULL);
    MPRINTF("\t\tTakes two options: 0 to disable, 1 to enable.\n", NULL);
    MPRINTF("\t\tThe first force module in the pipeline should set cleara=1.\n", NULL);
    MPRINTF("\t- plummer: Set a plummer distance for potential softening.\n", NULL);
    MPRINTF("\t\tTakes a double value.  Should be used if we're dealing with point particles.\n", NULL);
    MPRINTF("\t\tThis shouldn't be necessary if we're using particles with physical size, e.g. hscollide.\n", NULL);
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
        if(s->bodies[i].flags & PARTICLE_FLAG_DELETE) { continue; }
        for(int j = i + 1; j < s->nbody; j++) {
            if(s->bodies[j].flags & PARTICLE_FLAG_DELETE) { continue; }
            Vector r;
            vector_sub(&r, &s->bodies[i].pos, &s->bodies[j].pos);
            double f = pow(vector_dot(&r, &r) + cfg->plummer2,-1.5);    // note: this pow() takes about 75% of total compute time
                                                                        // consider pre-computing a table?
            s->bodies[i].acc.x -= cfg->G * s->bodies[j].mass * f * r.x;
            s->bodies[i].acc.y -= cfg->G * s->bodies[j].mass * f * r.y;
            s->bodies[i].acc.z -= cfg->G * s->bodies[j].mass * f * r.z;
            s->bodies[j].acc.x += cfg->G * s->bodies[i].mass * f * r.x;
            s->bodies[j].acc.y += cfg->G * s->bodies[i].mass * f * r.y;
            s->bodies[j].acc.z += cfg->G * s->bodies[i].mass * f * r.z;
        }
    }
    
    return MOD_RET_OK;                      // Return value can control flow of overall execution, see MOD_RET_*
}
