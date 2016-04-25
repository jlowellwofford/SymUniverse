//
//  ufromcsv.c - Convert CSV to Universe
//  SymUniverse
//
//  Created by J. Lowell Wofford on 4/14/16.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "universe.h"

#define NUM_FIELDS 21

int main(int argc, char *argv[]) {
    if(argc != 3) {
        printf("Usage: %s <in_file> <universe_file>\n", argv[0]);
        return 1;
    }
    FILE *i = fopen(argv[1], "r");
    if(i == NULL) {
        printf("Unable to open input file: %s\n", argv[1]);
        return 1;
    }
    Universe *u = universe_create(argv[2]);
    if(u == NULL) {
        printf("Unable to create universe file: %s\n", argv[2]);
        fclose(i);
        return 1;
    }
    
    Slice s;
    s.bodies = malloc(sizeof(Particle));
    s.time = 0;
    s.nbody = 0;
    s.bound_min.x = 0;
    char buff[4096];
    int line = 0;
    while(fgets(buff, 4096, i) != NULL) {
        ++line;
        Slice ts;
        Particle p;
        
        int ret = sscanf(buff,
                        "%llu,"
                        "%lg,%lg,%lg,"
                        "%lg,%lg,%lg,"
                        "%x,%x,"
                        "%lg,%lg,%lg,"
                        "%lg,%lg,%lg,"
                        "%lg,%lg,%lg,"
                        "%lg,%lg,%lg\n",
                        &ts.time,
                        &ts.bound_min.x, &ts.bound_min.y, &ts.bound_min.z,
                        &ts.bound_max.x, &ts.bound_max.y, &ts.bound_max.z,
                        &p.flags, &p.uflags,
                        &p.mass, &p.charge, &p.radius,
                        &p.pos.x, &p.pos.y, &p.pos.z,
                        &p.vel.x, &p.vel.y, &p.vel.z,
                        &p.acc.x, &p.acc.y, &p.acc.z
                        );
        if(ret != NUM_FIELDS) {
            printf("Error parsing line %d, expected %d, got %d.\n", line, NUM_FIELDS, ret);
            continue;
        }
        if(ts.time != s.time) { // new slice
            universe_append_slice(u, &s);
            s.time = ts.time;
            s.nbody = 0;
            s.bound_min.x = 0;
        }
        s.bound_min.x = ts.bound_min.x; s.bound_min.y = ts.bound_min.y; s.bound_min.z = ts.bound_min.z;
        s.bound_max.x = ts.bound_max.x; s.bound_max.y = ts.bound_max.y; s.bound_max.z = ts.bound_max.z;
        slice_append_particle(&s, &p);
    }
    universe_append_slice(u, &s);
    
    free(s.bodies);

    fclose(i);
    universe_close(u);
    return 0;
}