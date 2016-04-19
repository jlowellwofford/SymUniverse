//
//  sym.c
//  SymUniverse - Universe simulator
//
//  The function of gravgas is very simple.  It just iterates through timesteps, doing the following:
//      1. Copy current slice to a new slice
//      2. Modify the new slide through a series of modular transformations
//      3. Append the slice to the universe file
//  All of the physics happens in the transformation modules.
//
//  Created by J. Lowell Wofford on 3/25/16.
//  Copyright © 2016 J. Lowell Wofford. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <string.h>
#include <dlfcn.h>
#include <signal.h>
#include "universe.h"
#include "sym.h"
#include "SymUniverseConfig.h"

#define DEFAULT_IN_FILE "in.univ"
#define DEFAULT_OUT_FILE "out.univ"
#define DEFAULT_TIMESTEPS -1

struct {
    const char      *in_file;
    const char      *out_file;
    const char      *module_path;
    int             nmodules;
    Module          *modules;
    int             npipeline;
    Module          *pipeline;
    Universe        *universe;
    int             timesteps;
} cfg;
int exit_loop;
int sigint_caught;

void help(const char *cmd) {
    printf("\n"
           "This utility is used to simulate a gas with various parameters.\n"
           "In reality, the program just does the following:\n"
           "\t1. Load a Universe file and pull out the last Slice.\n"
           "\t2. Copy current Slice to new Slice.\n"
           "\t3. Perform a set of transforms defined by modules on the new Slice.\n"
           "\t4. Append the new Slice to the Universe and set it as current.\n"
           "\t5. Goto #2.\n"
           "All of the physics is handled in various modules.\n"
           "\n"
           "Usage: %s [options]\n"
           "Options:\n"
           "\t-h : Print this help.\n"
           "\t-i <file> : In universe file (default: %s).\n"
           "\t-o <file> : Out universe file (default: %s).\n"
           "\t-M <dir> : Directory to modules (default: %s).\n"
           "\t-t <steps> : Number of timesteps. -1 = infinite (default: %d)\n"
           "\t\t Simulations can always be safely stopped with Ctrl^c (a second Ctrl^c forces immedate abort).\n"
           "\n"
           "-- Module Help --\n"
           "To add modules to the pipeline, use the syntax:\n"
           "\t-m <mod_name>[module,option,string]\n"
           "\te.g. -m integrate[method=leapfrog,boundary=periodic]\n"
           "Modules are added to the pipeline in command line order.\n"
           "As a general rule add modules in this order: forces, integrate, collisions.\n"
           "Details on specific modules and options are below:\n\n",
           cmd, DEFAULT_IN_FILE, DEFAULT_OUT_FILE, DEFAULT_MODULE_PATH, DEFAULT_TIMESTEPS
    );
    for(int i = 0; i < cfg.nmodules; i++) {
        printf("Module name: %s\n", cfg.modules[i].name);
        cfg.modules[i].help();
        printf("\n");
    }
}

int load_module(char *path) {
    ++cfg.nmodules;
    if((cfg.modules = realloc(cfg.modules, sizeof(Module) * cfg.nmodules)) == NULL) {
        printf("Memory allocation failure.\n");
        return 0;
    }
    Module *m = &cfg.modules[cfg.nmodules - 1];

    m->handle = dlopen(path, RTLD_LAZY);
    if(m->handle == NULL) {
        printf("Could not open module at %s.\n%s\n", path, dlerror());
        exit(-1);
    }
    m->cfg = NULL;
    m->name = *(char **)dlsym(m->handle, "name");
    m->init = dlsym(m->handle, "init");
    m->deinit = dlsym(m->handle, "deinit");
    m->help = dlsym(m->handle, "help");
    m->exec = dlsym(m->handle, "exec");
    return 1;
}

void load_modules() {
    DIR *d;
    struct dirent *dir;
    d = opendir(cfg.module_path);
    if(d) {
        while((dir = readdir(d)) != NULL) {
            if(dir->d_type != DT_REG) { continue; }
            if(dir->d_namlen < 5) { continue; }
            if(strncmp(&dir->d_name[dir->d_namlen - 4], ".mod", 4) != 0) { continue; }
            char *full_path = calloc(strlen(cfg.module_path) + dir->d_namlen + 2, sizeof(char));
            if(full_path == NULL) {
                printf("Memory allocation failure.\n");
                exit(-1);
            }
            sprintf(full_path, "%s/%s", cfg.module_path, dir->d_name);
            if(!load_module(full_path)) {
                exit(-1);
            }
            free(full_path);
        }
    }
}

void unload_modules() { // Should be called after pipeline is destroyed!
    for(int i = 0; i < cfg.nmodules; i++) {
        dlclose(cfg.modules[i].handle);
    }
    cfg.nmodules = 0;
    free(cfg.modules);
}

void free_pipeline() {
    for(int i = 0; i < cfg.npipeline; i++) {
        cfg.pipeline[i].deinit(cfg.pipeline[i].cfg);
    }
    cfg.npipeline = 0;
    free(cfg.pipeline);
}

Module *find_module_by_name(const char *name) {
    for(int i = 0; i < cfg.nmodules; i++) {
        if(strcmp(name, cfg.modules[i].name) == 0) {
            return &cfg.modules[i];
        }
    }
    return NULL;
}

void init_pipeline(int argc, char *argv[]) {
    cfg.pipeline = realloc(cfg.pipeline, sizeof(Module) * argc);
    if(cfg.pipeline == NULL) {
        printf("Memory allocation error.\n");
        exit(-1);
    }
    printf("Pipeline: ");
    for(int i = 0; i < argc; i++) {
        char *margv = argv[i];
        char *mname = strsep(&margv, "[");
        margv = strsep(&margv, "]");
        Module *m = find_module_by_name(mname);
        memcpy(&cfg.pipeline[i], m, sizeof(Module));
        if((cfg.pipeline[i].cfg = cfg.pipeline[i].init(margv)) == NULL) {
            printf("Initialization of pipeline module, %s, failed!\n", mname);
            exit(-1);
        }
        ++cfg.npipeline;
        printf(" %s %s", mname, (argc == i + 1) ? "" : "->");
    }
    printf("\n");
}

void catch_SIGINT(int sig) {
    ++sigint_caught;
    if(sigint_caught == 1) {
        printf("\nCaught SIGINT, will exit after this loop iteration. SIGINT again to exit now.\n");
        exit_loop = 1;
    } else {
        printf("\nCaught second SIGINT, exiting now!\n");
        exit(0);
    }
}

void _universe_close() {    // Wrapper for atexit
    universe_close(cfg.universe);
}

int main(int argc, const char * argv[]) {   // Entry point
    // Parse args
    cfg.in_file = DEFAULT_IN_FILE;
    cfg.out_file = DEFAULT_OUT_FILE;
    cfg.module_path = DEFAULT_MODULE_PATH;
    cfg.timesteps = DEFAULT_TIMESTEPS;
    cfg.nmodules = 0;
    cfg.npipeline = 0;
    cfg.modules = malloc(sizeof(Module));
    if(cfg.modules == NULL) {
        printf("Memory allocation error.\n");
        exit(-1);
    }
    atexit(unload_modules);
    cfg.pipeline = malloc(sizeof(Module));
    if(cfg.pipeline == NULL) {
        printf("Memory allocation error.\n");
        exit(-1);
    }
    atexit(free_pipeline);

    printf("\n%s Version %d.%d by %s\n",
            SymUniverse_PROJECT_NAME, SymUniverse_VERSION_MAJOR, 
            SymUniverse_VERSION_MINOR, SymUniverse_PROJECT_AUTHOR);
    
    int ch, h = 0, pipe_argc = 0;
    char **pipe_argv = malloc(sizeof(char *));
    if(pipe_argv == NULL) {
        printf("Memory allocation error.\n");
        exit(-1);
    }
    while((ch = getopt(argc, (char * const *)argv, "hi:o:M:p:m:t:")) != -1) {
        switch (ch) {
            case 'i':
                cfg.in_file = optarg;
                break;
            case 'o':
                cfg.out_file = optarg;
                break;
            case 'M':
                cfg.module_path = optarg;
                break;
            case 'm':
                ++pipe_argc;
                pipe_argv = realloc(pipe_argv, sizeof(char *) * pipe_argc);
                if(pipe_argv == NULL) {
                    printf("Memory allocation error.\n");
                    exit(-1);
                }
                pipe_argv[pipe_argc - 1] = optarg;
                break;
            case 't':
                cfg.timesteps = atoi(optarg);
                break;
            case '?':
            case 'h':
            default:
                h = 1;
        }
    }
    
    load_modules();     // Load modules
    if(h != 0) {     // We want to load modules before displaying help so we can give mod help
        help(argv[0]);
        exit(-1);
    }
    // Initialize pipeline
    init_pipeline(pipe_argc, pipe_argv);
    free(pipe_argv);
    
    // Load Univere(s)
    if(access(cfg.out_file, W_OK) == 0) {
        printf("Opening existing file for output: %s\n", cfg.out_file);
        cfg.universe = universe_open(cfg.out_file);
    } else {
        printf("Creating new file for output: %s\n", cfg.out_file);
        cfg.universe = universe_create(cfg.out_file);

    }
    if(cfg.universe == NULL) {
        exit(-1);
    }
    atexit(_universe_close);
    
    Slice *pslice, *slice;
    
    if(strcmp(cfg.out_file, cfg.in_file) == 0) {
        printf("Output and Input files are the same, resuming from last slice.\n");
        pslice = universe_get_last_slice(cfg.universe);
    } else {
        printf("Reading initial configuration from: %s\n", cfg.in_file);
        Universe *iu = universe_open(cfg.in_file);
        pslice = universe_get_last_slice(iu);
        universe_close(iu);
    }
    if(pslice == NULL) {
        exit(-1);
    }
    slice = slice_copy(pslice);
    
    // Main loop
    int loop_idx = 0;
    exit_loop = (cfg.timesteps == 0) ? 1 : 0;
    sigint_caught = 0;
    signal(SIGINT, catch_SIGINT);
    setbuf(stdout, NULL);
    while(exit_loop == 0) {
        printf("\033[2K\rTimestep: %d/%d", loop_idx + 1, cfg.timesteps);
        int ret = 0;
        for(int i = 0; i < cfg.npipeline; i++) {
            ret |= cfg.pipeline[i].exec(cfg.pipeline[i].cfg, pslice, slice);
        }
        
        if(ret & MOD_RET_ABRT) { exit(-1); }            // Module requested abort without append
        if(ret & MOD_RET_EXIT) { exit_loop = 1; }       // Module requested exit
        if(ret & MOD_RET_PACK) {                        // Module thinks we need to repack
            if(!slice_pack(slice)) {
                exit(-1);                               // slice_pack failed, probably due to memory allocation.
            }
        }
        else { slice_clear_create(slice); }             // This is redundant if we run slice_pack
                                                        // ???: Could make clear_create based on ret value.
        
        if(!universe_append_slice(cfg.universe, slice)) {
            exit(-1);
        }
        slice_free(pslice);
        pslice = slice;
        slice = slice_copy(pslice);
        if(slice == NULL) {
            printf("Memory allocation error.\n");
            exit(-1);
        }
        ++slice->time;
        ++loop_idx;
        if(loop_idx >= cfg.timesteps) {
            exit_loop = 1;
        }
    }
    printf("\n");
    
    // TODO: atexit doesn't take care of these.
    slice_free(pslice);
    slice_free(slice);
    
    return 0;
}
