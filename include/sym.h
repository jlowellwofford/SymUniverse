//
//  sym.h
//  SymUniverse - misc structure/type declarations
//
//  Created by J. Lowell Wofford on 3/25/16.
//  Copyright © 2016 J. Lowell Wofford. All rights reserved.
//

#ifndef sym_h 
#define sym_h 

#include "universe.h"

#define DEFAULT_MODULE_PATH "modules/"

#define MPRINTF(f_, ...) printf("[%s] " f_, name, __VA_ARGS__)

// Module return flags
#define MOD_RET_OK   0      // Execution was OK, do nothing special
#define MOD_RET_ABRT 1      // Abort without appending current slice (i.e. something is wrong)
#define MOD_RET_EXIT 2      // Exit after appending slice (e.g. we met a finalization condition)
#define MOD_RET_PACK 4      // Usually means module marked some particles for deletion

// Modules must have the following symbols:
// 1. name (function)   returns the name of the module
// 2. init (function)   initialize module parameters for use in pipeline
// 3. help (function)   prints help info for the module
// 4. exec (function)   execute the module transorm

typedef struct {
    void        *handle;
    void        *cfg;
    const char  *name;
    void        (*help)(void);
    void        *(*init)(char *cfg_str);
    void        (*deinit)(void *cfg);
    int         (*exec)(void *cfg, Slice *ps, Slice *s);
} Module;

#endif /* sym_h */
