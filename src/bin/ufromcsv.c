//
//  ufromcsv.c - Convert CSV to Universe
//  SymUniverse
//
//  Created by J. Lowell Wofford on 4/14/16.
//
//

#include <stdio.h>
#include "universe.h"

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
    
    // we always assume the first line is the header
    
    
    fclose(i);
    universe_close(u);
    return 0;
}