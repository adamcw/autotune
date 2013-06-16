#ifndef BIN_BHEAP_H
#define BIN_BHEAP_H

#include "../dllist/dllist.h"

/** \define group bheap Binary Heap 
 * @{
 */

typedef struct bheap BHEAP;

struct bheap {
	int max_elem;
	int num_elem;
	void **k;
};

BHEAP *bh_create(int init_max_elem);
void bh_reset(BHEAP *h, void (*free_key)(void *k));
BHEAP *bh_copy(BHEAP *h, void *(*copy_key)(void *key), void (*set_key_parent)(void *k, BHEAP *h));

void bh_insert_key(BHEAP *h, void *k, int (*key_lt)(void *k1, void *k2), void (*swap_keys)(BHEAP *h, int i1, int i2), void (*set_key_index)(void *k, int i), void (*set_key_parent)(void *k, BHEAP *h)); 

void *bh_get_top(BHEAP *h);

void bh_bubble_down_key(BHEAP *h, int i, int (*key_lt)(void *k1, void *k2), void (*swap_keys)(BHEAP *h, int i1, int i2));
void bh_bubble_up_key(BHEAP *h, int i, int (*key_lt)(void *k1, void *k2), void (*swap_keys)(BHEAP *h, int i1, int i2));

void *bh_remove_key(BHEAP *h, int i, int (*key_lt)(void *k1, void *k2), void (*swap_keys)(BHEAP *h, int i1, int i2), void (*set_key_index)(void *k, int i), void (*set_key_parent)(void *k, BHEAP *h));

void bh_validate(BHEAP *h, int (*key_lt)(void *k1, void *k2), void (*print_key)(void *k));

void bh_print_top(BHEAP *h, void (*print_key)(void *k)); 
void bh_print(BHEAP *h, void (*print_key)(void *k));

void bh_free(BHEAP *h, void (*free_key)(void *k));

/** @} */

#endif
