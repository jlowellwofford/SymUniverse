//
//  utocsv.c - Universe to CSV converter
//  SymUniverse
//
//  Created by J. Lowell Wofford on 4/14/16.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "universe.h"

int main(int argc, char *argv[]) {
    
    if(argc < 3 || argc > 4) {
        printf("Usage: %s [interval] <universe_file> <out_file>\n", argv[0]);
        return 1;
    }
    
    int argn = 1;
    int interval = 1;
    if(argc == 4) {
        interval = atoi(argv[argn++]);
    }
    
    Universe *u = universe_open(argv[argn++]);
    if(u == NULL) {
        printf("Unable to open universe file: %s\n", argv[1]);
        return 1;
    }
    FILE *o = fopen(argv[argn++], "w+");
    if(o == NULL) {
        printf("Unable to open output file: %s\n", argv[2]);
        universe_close(u);
        return 1;
    }
    
    fprintf(o, "time,min.x,min.y,min.z,max.x,max.y,max.z,flags,uflags,mass,charge,radius,"
            "pos.x,pos.y,pos.z,vel.x,vel.y,vel.z,acc.x,acc.y,acc.z\n");
            
    for(int k = 0; k < floor((float)u->nslice/interval); k++) {
        int i = k * interval;
        Slice *s = universe_get_slice(u, i);
        if(s == NULL) {
            printf("Unable to read slice %d.\n", i);
            fclose(o);
            universe_close(u);
            return 1;
        }
        for(int j = 0; j < s->nbody; j++) {
            fprintf(o,
                    "%llu,"
                    "%g,%g,%g,"
                    "%g,%g,%g,"
                    "%x,%x,"
                    "%g,%g,%g,"
                    "%g,%g,%g,"
                    "%g,%g,%g,"
                    "%g,%g,%g\n",
                    s->time,
                    s->bound_min.x, s->bound_min.y, s->bound_min.z,
                    s->bound_max.x, s->bound_max.y, s->bound_max.z,
                    s->bodies[j].flags, s->bodies[j].uflags,
                    s->bodies[j].mass, s->bodies[j].charge, s->bodies[j].radius,
                    s->bodies[j].pos.x, s->bodies[j].pos.y, s->bodies[j].pos.z,
                    s->bodies[j].vel.x, s->bodies[j].vel.y, s->bodies[j].vel.z,
                    s->bodies[j].acc.x, s->bodies[j].acc.y, s->bodies[j].acc.z
                    );
        }
        slice_free(s);
    }
    
    fclose(o);
    universe_close(u);
    return 0;
}
