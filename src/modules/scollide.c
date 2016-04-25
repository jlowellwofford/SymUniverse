//
//  scollide.c
//  SymUniverse - Module for doing simple sphere collision detection and resolution.  Not very physically accurate.
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
const char *name = "scollide";      // Name _must_ be unique

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
        printf("Memory allocation failure.\n");
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
    MPRINTF("This module resolves sphere collisions.\n", NULL);
    MPRINTF("This simple algorithm is O(N^2), and takes no options.\n", NULL);
    MPRINTF("Note: this module isn't very good about conserving physical quantities, but the differences should average out over time.\n", NULL);
    MPRINTF("Typically this should be placed after forces and integration.\n", NULL);
    MPRINTF("For best results, disable boundary detection in integrate and use the boundary module after this.\n", NULL);
}

// !!!: It would be nice to have a recursive, time-ordered collision detection option.
// FIXME: This algorithm is not yet implemented!
EXPORT
int exec(Config *cfg, Slice *ps, Slice *s) {  // Main execution loop.  Maps (ps, s) -> s.  Should _not_ modify ps.  Uses cfg to specify pipeline params.
    Vector p, pp, v;
    double ts = 0;
    if(s->nbody < 1) { return MOD_RET_OK; }
    for(int i = 0; i < ps->nbody; i++) {    // A little odd.  Basically, we want to account for the fact that some particles may have zero v.x;
        if(s->bodies[i].vel.x == 0) { continue; }
        ts = (s->bodies[i].pos.x - ps->bodies[i].pos.x) / s->bodies[i].vel.x; // Reverse engineer the timestep.  Only need to do this once.
        break;
    }
    if(!ts) { return MOD_RET_OK; }  // Everything in our sym seems to be sitting still! (at least in the x direction)
    
    int ccount = 0;
    for(int i = 0; i < ps->nbody; i++) {
        if(s->bodies[i].flags & (PARTICLE_FLAG_DELETE | PARTICLE_FLAG_CREATE)) { continue; }
        for(int j = i+1; j < ps->nbody; j++) {
            if(s->bodies[j].flags & (PARTICLE_FLAG_DELETE | PARTICLE_FLAG_CREATE)) { continue; }
            
            // 1. move into the rest frame of body i
            vector_sub(&p, &s->bodies[j].pos, &ps->bodies[i].pos);
            vector_sub(&pp, &ps->bodies[j].pos, &ps->bodies[i].pos);
            vector_sub(&v, &s->bodies[j].vel, &s->bodies[i].vel);
            
            // TODO 1b. early collision elimination
            //if(0) { continue; }
            
            // 2. move to scattering frame
            double v2 = vector_dot(&v, &v);
            if(v2 == 0) { continue; }
            double t0 = - vector_dot(&v, &pp) / (v2 * ts);  // time to perpendicular bisector to infinite line defined by pp -> p

            Vector chi;                                     // vector defining bisector
            chi.x = v.x * t0 + pp.x;
            chi.y = v.y * t0 + pp.y;
            chi.z = v.z * t0 + pp.z;
            double b = sqrt(vector_dot(&chi, &chi));        // impact parameter
            double R = s->bodies[i].radius + s->bodies[j].radius;
            if(b > R) { continue; }  // definitely don't collide
            
            double vel = sqrt(v2);                          // note: this is always >= 0
            double xi = - vel * t0;
            if(xi >= 0) { continue; }                        // particle is moving away
            double xf = vel * (ts - t0);
            if(xf < - R) { continue; }
            
            // 3. past this point, we know there's a collision
            ++ccount;
            // The "theory" here is that v1f should be along the line connecting the centers at the point of collision.
            // The magnitude of v1f is then given by v1f = 2*vel*j.mass/(i.mass+j.mass)*cos(theta)
            // We then just adjust v1f to conserve linear momentum.
            double tc = - (xi + R) / vel;                   // time of collision
            
            double vifx = 2*vel*s->bodies[j].mass/(s->bodies[j].mass + s->bodies[i].mass)*(R*R-b*b)/(R*R);
            double vify = - 2*vel*s->bodies[j].mass/(s->bodies[j].mass + s->bodies[i].mass)*sqrt(R*R-b*b)/(R*R)*b;
            double vjfx = - s->bodies[i].mass/s->bodies[j].mass*vifx + vel;
            double vjfy = - s->bodies[i].mass/s->bodies[j].mass*vify;
            double xifx = vifx * (ts - tc);
            double xify = vify * (ts - tc);
            double xjfx = vjfx * (ts - tc) - sqrt(R*R-b*b);
            double xjfy = vjfx * (ts - tc) + b;
            
            // 4. return to sym frame
            if(b != 0) {
                Vector ux, uy; //,uz;      // together as a row matrix, these give the inverse rotation matrix. (uz is unused)
                ux.x = v.x / vel;
                ux.y = v.y / vel;
                ux.z = v.z / vel;
                uy.x = chi.x / b;
                uy.y = chi.y / b;
                uy.z = chi.z / b;
                // vector_cross(&uz, &ux, &uy);
                // rotate + translate in one step
                s->bodies[j].pos.x = xjfx * ux.x + xjfy * uy.x + s->bodies[i].pos.x;
                s->bodies[j].pos.y = xjfx * ux.y + xjfy * uy.y + s->bodies[i].pos.y;
                s->bodies[j].pos.z = xjfx * ux.z + xjfy * uy.z + s->bodies[i].pos.z;
                s->bodies[i].pos.x += xifx * ux.x + xify * uy.x;
                s->bodies[i].pos.y += xifx * ux.y + xify * uy.y;
                s->bodies[i].pos.z += xifx * ux.z + xify * uy.z;
                
                s->bodies[j].vel.x = vjfx * ux.x + vjfy * uy.x + s->bodies[i].vel.x;
                s->bodies[j].vel.y = vjfx * ux.y + vjfy * uy.y + s->bodies[i].vel.y;
                s->bodies[j].vel.z = vjfx * ux.z + vjfy * uy.z + s->bodies[i].vel.z;
                s->bodies[i].vel.x += vifx * ux.x + vify * uy.x;
                s->bodies[i].vel.y += vifx * ux.y + vify * uy.y;
                s->bodies[i].vel.z += vifx * ux.z + vify * uy.z;
            } else {                    // in the off chance b = 0, we have a divide by zero, we have to do things differently
                MPRINTF("Oops, we got a collision we couldn't handle.\n", NULL);
//                Vector uv;
//                uv.x = v.x / vel;       // Both position & velocity will lie relative to the velocity unit vector
//                uv.y = v.y / vel;
//                uv.z = v.z / vel;
//                
//                s->bodies[j].pos.x = x2f * uv.x + s->bodies[i].pos.x;
//                s->bodies[j].pos.y = x2f * uv.y + s->bodies[i].pos.y;
//                s->bodies[j].pos.z = x2f * uv.z + s->bodies[i].pos.z;
//                s->bodies[i].pos.x += x1f * uv.x;
//                s->bodies[i].pos.y += x1f * uv.y;
//                s->bodies[i].pos.z += x1f * uv.z;
//                s->bodies[j].vel.x = v2f * uv.x + s->bodies[i].vel.x;
//                s->bodies[j].vel.y = v2f * uv.y + s->bodies[i].vel.y;
//                s->bodies[j].vel.z = v2f * uv.z + s->bodies[i].vel.z;
//                s->bodies[i].vel.x += v1f * uv.x;
//                s->bodies[i].vel.y += v1f * uv.y;
//                s->bodies[i].vel.z += v1f * uv.z;
            }
            
        }
    }
    
    MPRINTF("Processed %d collisions.\n", ccount);
    
    return MOD_RET_OK;                      // Return value can control flow of overall execution, see MOD_RET_*
}
