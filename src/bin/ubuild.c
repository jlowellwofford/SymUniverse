//
//  ubuild.c
//  SymUniverse - Universe Builder
//
//  Created by J. Lowell Wofford on 3/25/16.
//  Copyright © 2016 J. Lowell Wofford. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <string.h>
#include <math.h>
#include "universe.h"
#include "SymUniverseConfig.h"

#define DEFAULT_OUT "out.univ"
#define DEFAULT_TEMP 300
#define DEFAULT_NBODY 1000
#define DEFAULT_MASS 938            // Mass of the proton in MeV
#define DEFAULT_RADIUS 5E-11        // Bohr radius (m)
#define DEFAULT_CHARGE 0
#define DEFAULT_BOUND_MIN 0
#define DEFAULT_BOUND_MAX 1
#define DEFAULT_BOX_MIN 0
#define DEFAULT_BOX_MAX 1
#define DEFAULT_FLAGS 0

struct {
    const char  *out_file;
    double      temp;
    uint64_t    nbody;
    double      mass;
    double      radius;
    double      charge;
    Vector      bound_min;
    Vector      bound_max;
    Vector      box_min;
    Vector      box_max;
    uint64_t    flags;
} cfg;

void help(const char *cmd) {
    // Print help
    printf("\n"
           "This utility creates an initial universe file for SymUniverse.\n"
           "It assumes the following:\n"
           "\t- Monotomic gas (equal mass and radius for each particle).\n"
           "\t- Zero initial acceleration.\n"
           "\t- Uniform, random spacial configuration.\n"
           "\t- Maxwell-Boltzmann distribution of velocities.\n"
           "\n"
           "Usage: %s [options]\n"
           "Options:\n"
           "\t-h : Print this help.\n"
           "\t-o <file_name> : Output file (default: %s)\n"
           "\t-T <temp> : Distribution temperature (default: %d)\n"
           "\t-n <num> : Number of particles (default: %u)\n"
           "\t-m <mass> : Mass of particles (default: %d)\n"
           "\t-c <charge> : Charge of the particles (default: %d)\n"
           "\t-r <radius> : Radius of particles (default: %f)\n"
           "\t-b <bound_spec> : Boundary corners, format xmin,ymin,zmin,xmax,ymax,zmax of min,max for a cube.\n"
           "\t-B <box_spec> : Box corners, format xmin,ymin,zmin,xmax,ymax,zmax of min,max for a cube.\n"
           "\t\tThe box allows particles to be contained in a box, even if the boundary is larger.\n"
           "\t-f <flags> : Specify particle flags (see header file) (default: %ull)\n"
           "\n",
           cmd, DEFAULT_OUT, DEFAULT_TEMP, DEFAULT_NBODY, DEFAULT_MASS, DEFAULT_CHARGE,DEFAULT_RADIUS,DEFAULT_FLAGS
    );
}

void maxwell_boltzmann(Particle *p) {
    double r = ((double)rand())/RAND_MAX;
    double sigma = sqrt(cfg.temp / p->mass);
    // the following comes from the inverse series of the m-b distribution up to 12th order
    double vt = sigma * (1.5499 * pow(r, 1/3) + 0.375994 * r + 0.201312 * pow(r, 5/3) + 0.134632 * pow(r, 7/3) + 0.100084 * pow(r,3));
    // pick a random orientation
    double theta = ((double)rand())/RAND_MAX * M_PI;
    double phi = ((double)rand())/RAND_MAX * 2 * M_PI;
    p->vel.x = vt * sin(theta) * cos(phi);
    p->vel.y = vt * sin(theta) * sin(phi);
    p->vel.z = vt * cos(theta);
}

void min_max_to_vector(char *minmax, Vector *min, Vector *max) {
    int argc = 0;
    char *argv[6];
    do {
        argv[argc] = strsep(&minmax, ",");
        ++argc;
    } while(argc < 6 && argv[argc - 1][0] != '\0');
    switch (argc) {
        case 2:
            min->x = min->y = min->z = strtod(argv[0], NULL);
            max->x = max->y = max->z = strtod(argv[1], NULL);
            break;
        case 6:
            min->x = strtod(argv[0], NULL);
            min->y = strtod(argv[1], NULL);
            min->y = strtod(argv[2], NULL);
            max->x = strtod(argv[3], NULL);
            max->y = strtod(argv[4], NULL);
            max->z = strtod(argv[5], NULL);
            break;
        default:
            printf("Min-Max specification must be of form: minx,miny,minz,maxx,maxy,maxz or min,max\n");
            exit(-1);
    }
}

int main(int argc, const char *argv[]) {   // Entry point

    cfg.out_file = DEFAULT_OUT;
    cfg.temp = DEFAULT_TEMP;
    cfg.nbody = DEFAULT_NBODY;
    cfg.mass = DEFAULT_MASS;
    cfg.radius = DEFAULT_RADIUS;
    cfg.charge = DEFAULT_CHARGE;
    cfg.bound_min.x = DEFAULT_BOUND_MIN;
    cfg.bound_min.y = DEFAULT_BOUND_MIN;
    cfg.bound_min.z = DEFAULT_BOUND_MIN;
    cfg.bound_max.x = DEFAULT_BOUND_MAX;
    cfg.bound_max.y = DEFAULT_BOUND_MAX;
    cfg.bound_max.z = DEFAULT_BOUND_MAX;
    cfg.box_min.x = DEFAULT_BOX_MIN;
    cfg.box_min.y = DEFAULT_BOX_MIN;
    cfg.box_min.z = DEFAULT_BOX_MIN;
    cfg.box_max.x = DEFAULT_BOX_MAX;
    cfg.box_max.y = DEFAULT_BOX_MAX;
    cfg.box_max.z = DEFAULT_BOX_MAX;
    cfg.flags = DEFAULT_FLAGS;
    
    printf("\n%s Version %d.%d by %s\n",
           SymUniverse_PROJECT_NAME, SymUniverse_VERSION_MAJOR,
           SymUniverse_VERSION_MINOR, SymUniverse_PROJECT_AUTHOR);
    
    // Parse args
    int ch;
    while((ch = getopt(argc, (char * const *)argv, "ho:T:n:m:r:c:b:B:f:")) != -1) {
        switch (ch) {   // TODO: sanity checks
            case 'o':
                cfg.out_file = optarg;
                break;
            case 'T':
                cfg.temp = strtod(optarg, NULL);
                break;
            case 'n':
                cfg.nbody = strtouq(optarg, NULL, 0);
                break;
            case 'm':
                cfg.mass = strtod(optarg, NULL);
                break;
            case 'r':
                cfg.radius = strtod(optarg, NULL);
                break;
            case 'c':
                cfg.charge = strtod(optarg, NULL);
                break;
            case 'b':
                min_max_to_vector(optarg, &cfg.bound_min, &cfg.bound_max);
                break;
            case 'B':
                min_max_to_vector(optarg, &cfg.box_min, &cfg.box_max);
                break;
            case 'f':
                cfg.flags = strtouq(optarg, NULL, 0);
                break;
            case '?':
            case 'h':
            default:
                help(argv[0]);
                exit(-1);
        }
    }

    // Generate universe
    Slice s;
    s.nbody = cfg.nbody;
    s.time = 0;
    s.bound_max.x = cfg.bound_max.x;
    s.bound_max.y = cfg.bound_max.y;
    s.bound_max.z = cfg.bound_max.z;
    s.bound_min.x = cfg.bound_min.x;
    s.bound_min.y = cfg.bound_min.y;
    s.bound_min.z = cfg.bound_min.z;
    s.bodies = (Particle *)calloc(cfg.nbody,sizeof(Particle));
    
    for(uint64_t b = 0; b < cfg.nbody; b++) {
        s.bodies[b].mass = cfg.mass;
        s.bodies[b].radius = cfg.radius;
        s.bodies[b].charge = cfg.charge;
        s.bodies[b].flags = cfg.flags;
        // Accel is automatically 0, since we used calloc
        sranddev();
        s.bodies[b].pos.x = ((double)rand())/RAND_MAX * (cfg.box_max.x - cfg.box_min.x) + cfg.box_min.x;
        s.bodies[b].pos.y = ((double)rand())/RAND_MAX * (cfg.box_max.y - cfg.box_min.y) + cfg.box_min.y;
        s.bodies[b].pos.z = ((double)rand())/RAND_MAX * (cfg.box_max.z - cfg.box_min.z) + cfg.box_min.z;
        maxwell_boltzmann(&s.bodies[b]);
    }
    // Write universe file
    Universe *u = universe_create(cfg.out_file);
    universe_append_slice(u, &s);
    
    free(s.bodies);
    universe_close(u);
    universe_free(u);
    
    return 0;
}
