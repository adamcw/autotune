#ifndef HASH_TABLE_H
#define HASH_TABLE_H

#include "../dllist/dllist.h"

/** \defgroup ht Hash Table 
 * @{ */

/**
 * \brief A hash table 
 */
typedef struct hash_table HT;
struct hash_table {
    /** An array of double-linked lists to store the keys in */
	DLL_NODE **table;

    /** The length (size) of the hash table */
	int length;

    /** The number of elements currently in the hash table */
	int num_elem;
};

HT *ht_create(int len);
DLL_NODE *ht_insert_key(HT *ht, int id, void *key);
HT *ht_copy(HT *ht, void *(*copy_key)(void *));
DLL_NODE *ht_hash_lookup(HT *ht, int id); 
DLL_NODE *ht_find_node(HT *ht, int id, void *key);
DLL_NODE *ht_find_node_fn(HT *ht, int id, void *key, int (*comparator)(void *key1, void *key2));
void ht_delete_key(HT *ht, int id, void *key, void (*free_key)(void *key));
void ht_print(HT *ht, void (*print_key)(void *key));
void ht_free(HT *ht, void (*free_key)(void *key));

/** @} */

#endif
