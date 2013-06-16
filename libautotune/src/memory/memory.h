#ifndef MEMORY_H
#define MEMORY_H

#include <stdlib.h>

int my_mem();

void *my_malloc(size_t size);
void *my_realloc(void *ptr, size_t size);

void *my_calloc(int n, int size);        
void *my_2d_calloc(int n1, int n2, int size);        
void *my_3d_calloc(int n1, int n2, int n3, int size);

void my_2d_free(int n1, void **pptr);
void my_3d_free(int n1, int n2, void ***ppptr);

#endif
