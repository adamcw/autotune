#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "bheap.h"
#include "../memory/memory.h"
#include "../dllist/dllist.h"

#define TRUE 1
#define FALSE 0

#define PARENT(i) ((i)>>1)
#define LEFT(i) ((i)<<1)
#define RIGHT(i) (((i)<<1)+1)

#define ASSERT(cond) assert((cond))
//#define ASSERT(cond)

/**
 * \brief Creates a new minimum \ref heap
 * 
 * \param[in] init_max_elem The initial number of elements to allocate memory
 * for.
 *
 * \return The newly created \ref heap
 */
BHEAP *bh_create(int init_max_elem) {
	BHEAP *h;

	h = (BHEAP *)my_malloc(sizeof(BHEAP));
	
	h->max_elem = init_max_elem;
	h->num_elem = 0;
	h->k = (void **)my_malloc(init_max_elem*sizeof(void *));

	return h;
}

/**
 * \brief Resets a binary heap to its initial state
 * 
 * \param[in] h The \ref bheap to be reset
 * \param[in] free_key A function to be called on each key in the heap to free
 * and memory associated with it. Can be NULL.
 */
void bh_reset(BHEAP *h, void (*free_key)(void *k)) {
	int i;

	if (free_key != NULL) {
		for (i = 1; i <= h->num_elem; i++) {
			free_key(h->k[i]);
		}
	}

	h->num_elem = 0;
}

/**
 * \brief Copies a \ref heap
 * 
 * \param[in] h The \ref heap to be copied
 * \param[in] copy_key A function that accepts a key and returns a copy of it
 * \param[in] set_key_parent A function that takes a key and a \ref heap and set a
 * line from the key to the \ref heap. 
 *
 * \return The newly copied \ref heap
 */
BHEAP *bh_copy(BHEAP *h, void *(*copy_key)(void *key), void (*set_key_parent)(void *k, BHEAP *h)) {
	BHEAP *h2;
	int i;

	h2  = bh_create(h->max_elem);
	for (i = 1; i <= h->num_elem; i++) {
		h2->k[i] = copy_key(h->k[i]);
		if (set_key_parent != NULL) {
			set_key_parent(h2->k[i], h2);
		}
	}
	h2->num_elem = h->num_elem;
	return h2;
}

/**
 * \brief Inserts a key into a \ref heap
 * 
 * \param[in] h The \ref heap to insert the key into
 * \param[in] k The key to insert
 * \param[in] key_lt A function that takes two keys and returns 1 if the first
 * is less than the second, 0 otherwise.
 * \param[in] swap_keys A function that takes a heap and two indices and swaps
 * the positions of the two indices.
 * \param[in] set_key_index A function that takes a key and an index and sets a
 * link from the key to the index.
 * \param[in] set_key_parent A function that takes a key and a \ref heap and set a
 * line from the key to the \ref heap. 
 */
void bh_insert_key(BHEAP *h, void *k, int (*key_lt)(void *k1, void *k2), void (*swap_keys)(BHEAP *h, int i1, int i2), void (*set_key_index)(void *k, int i), void (*set_key_parent)(void *k, BHEAP *h)) {
	h->num_elem++;

    // Resize the keys array if required
	if (h->num_elem >= h->max_elem) {
		h->max_elem = h->num_elem;
		h->max_elem <<= 1;
		h->k = (void **)my_realloc(h->k, h->max_elem*sizeof(void *));
	}
	
    // Insert the key into the keys array
	h->k[h->num_elem] = k;
	if (set_key_index != NULL) {
		set_key_index(k, h->num_elem);
	}
	if (set_key_parent != NULL) {
		set_key_parent(k, h);
	}

	bh_bubble_up_key(h, h->num_elem, key_lt, swap_keys);
}

/**
 * \brief Gets the first key in a \ref heap
 * 
 * \param[in] h The \ref heap to get the first key of
 *
 * \return The key if the heap is non-empty, NULL otherwise.
 */
void *bh_get_top(BHEAP *h) {
	if (h->num_elem == 0) {
		return NULL;
	}
	return h->k[1];
}

/**
 * \brief Remove a key from a \ref heap
 * 
 * \param[in] h The \ref heap to remove the key from
 * \param[in] i The index of the key to remove in the \ref heap
 * \param[in] key_lt A function that takes two keys and returns 1 if the first
 * is less than the second, 0 otherwise.
 * \param[in] swap_keys A function that takes a heap and two indices and swaps
 * the positions of the two indices.
 * \param[in] set_key_index A function that takes a key and an index and sets a
 * link from the key to the index.
 * \param[in] set_key_parent A function that takes a key and a \ref heap and set a
 * line from the key to the \ref heap. 
 */
void *bh_remove_key(BHEAP *h, int i, int (*key_lt)(void *k1, void *k2), void (*swap_keys)(BHEAP *h, int i1, int i2), void (*set_key_index)(void *k, int i), void (*set_key_parent)(void *k, BHEAP *h)) {
	void *key;
	int end;

	ASSERT(i <= h->num_elem);
	ASSERT(i > 0);

	end = h->num_elem;

    key = h->k[i];
    h->k[i] = h->k[end];
    if (set_key_index != NULL) {
        set_key_index(h->k[i], i);
    }
    h->num_elem--;
    
    if (i != end) {
        if ((i > 1) && (key_lt(h->k[i], h->k[(i>>1)]) == TRUE)) {
            bh_bubble_up_key(h, i, key_lt, swap_keys);
        }
        else {
            bh_bubble_down_key(h, i, key_lt, swap_keys);
        }
    }
	
	if (set_key_parent != NULL) {
		set_key_parent(key, NULL);
	}
	if (set_key_index != NULL) {
		set_key_index(key, 0);
	}

	return key;
}

/**
 * \brief Bubbles a key down a \ref heap 
 * 
 * \param[in] h The \ref heap to bubble the key in
 * \param[in] i The index of the key to bubble in the \ref heap
 * \param[in] key_lt A function that takes two keys and returns 1 if the first
 * is less than the second, 0 otherwise.
 * \param[in] swap_keys A function that takes a heap and two indices and swaps
 * the positions of the two indices.
 */
void bh_bubble_down_key(BHEAP *h, int i, int (*key_lt)(void *k1, void *k2), void (*swap_keys)(BHEAP *h, int i1, int i2)) {
	int i_L, i_R, i_min, end;

	end = h->num_elem;
	if (end == 0) {
		return;
	}

	i_L = LEFT(i);
	i_R = RIGHT(i);

	if ((i_L <= end) && key_lt(h->k[i_L], h->k[i])) {
		i_min = i_L;
	}
	else {
		i_min = i;
	}

	if ((i_R <= end) && key_lt(h->k[i_R], h->k[i_min])) {
		i_min = i_R;
	}

	if (i_min != i) {
		swap_keys(h, i, i_min);
		bh_bubble_down_key(h, i_min, key_lt, swap_keys);
	}
}

/**
 * \brief Bubbles a key up a \ref heap 
 * 
 * \param[in] h The \ref heap to bubble the key in
 * \param[in] i The index of the key to bubble in the \ref heap
 * \param[in] key_lt A function that takes two keys and returns 1 if the first
 * is less than the second, 0 otherwise.
 * \param[in] swap_keys A function that takes a heap and two indices and swaps
 * the positions of the two indices.
 */
void bh_bubble_up_key(BHEAP *h, int i, int (*key_lt)(void *k1, void *k2), void (*swap_keys)(BHEAP *h, int i1, int i2)) {
	int i_P = PARENT(i);

	while (i > 1 && key_lt(h->k[i], h->k[i_P])) {
		swap_keys(h, i, i_P);
		i = PARENT(i);
		i_P = PARENT(i_P);
	}
}

/**
 * \brief Validates the properties of a \ref heap
 * 
 * \param[in] h The \ref heap to validate
 * \param[in] key_lt A function that takes two keys and returns 1 if the first
 * is less than the second, 0 otherwise.
 * \param[in] print_key A function that takes a key and prints its 
 */
void bh_validate(BHEAP *h, int (*key_lt)(void *k1, void *k2), void (*print_key)(void *k)) {
	int i;
	int elems;

	elems = h->num_elem;

	for (i = 1; i <= (elems / 2); i++) {
		if (key_lt(h->k[i<<1], h->k[i]) == TRUE) {
			printf("Assert: The left child of a heap key is less than the parent\n");
			printf("The current key is:\n");
			print_key(h->k[i]);
			printf("The left child key is:\n");
			print_key(h->k[i<<1]);
			printf("Dumping the heap\n");
			bh_print(h, print_key);
			assert(0);
		}
		if ((i<<1)+1 <= elems && key_lt(h->k[(i<<1)+1], h->k[i]) == TRUE) {
			printf("Assert: The right child of a heap key is less than the parent\n");
			printf("The current key is:\n");
			print_key(h->k[i]);
			printf("The right child key is:\n");
			print_key(h->k[(i<<1)+1]);
			printf("Dumping the heap\n");
			bh_print(h, print_key);
			assert(0);
		}
	}
}

/**
 * \brief Prints the key at the top of a \heap
 * 
 * An alias for ::bh_get_top then calling print_key on the result.
 *
 * \param[in] h The \ref heap to print the top of
 * \param[in] print_key A function that takes a key and prints it 
 */
void bh_print_top(BHEAP *h, void (*print_key)(void *k)) {
	print_key(bh_get_top(h));
}

/**
 * \brief Prints a \ref heap
 * 
 * \param[in] h The \ref heap to print 
 * \param[in] print_key A function that takes a key and prints it 
 */
void bh_print(BHEAP *h, void (*print_key)(void *k)) {
	int i;
	if (h->k != NULL) {
		printf("Heap of %d Keys\n", h->num_elem);
		printf("k[]\n");
		for (i = 1; i <= h->num_elem; i++) {
			print_key(h->k[i]);
		}
	}
}

/**
 * \brief Frees a heap
 * 
 * \param[in] h The \ref heap to free 
 * \param[in] free_key A function that takes a key and frees the memory associated with it 
 */
void bh_free(BHEAP *h, void (*free_key)(void *k)) {
	int i;

	if (h == NULL) {
		return;
	}

	if (h->k != NULL) {
		if (free_key != NULL) {
			for (i = 1; i <= h->num_elem; i++) {
				free_key(h->k[i]);
			}
		}
		free(h->k);
	}

    free(h);
}
