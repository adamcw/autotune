#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "../memory/memory.h"
#include "hash_table.h"

#define WARNING 0

/**
 * \brief Creates a new hash table
 *
 * \param[in] len The length of the hash table 
 *
 * \return The newly created hash table
 */
HT *ht_create(int len) {
	HT *ht;

	ht = (HT *)my_malloc(sizeof(HT));

	ht->table = (DLL_NODE **)my_calloc(len, sizeof(DLL_NODE *));
	ht->length = len;
	ht->num_elem = 0;

	return ht;
}

/**
 * \brief Adds an element to a hash table
 *
 * \param[in] ht The hash table to be inserted into
 * \param[in] id The index (hash) of where the object should be inserted 
 * \param[in] key The object which will be inserted into the hash table
 *
 * \return A DLL_NODE of the newly inserted element.
 */
DLL_NODE *ht_insert_key(HT *ht, int id, void *key) {
	int i;

	i = id % ht->length;
	ht->table[i] = dll_insert(ht->table[i], key, NULL);
	ht->num_elem++;
	if (WARNING && ht->num_elem > ht->length) {
		printf("warning: ht->num_elem = %d > ht->length = %d\n", ht->num_elem, ht->length);
	}

	return ht->table[i];
}

/**
 * \brief Copies a hash table
 *
 * \param[in] ht The hash table to be inserted into
 * \param[in] copy_key A function that takes a key and copies it.
 *
 * \return The newly copied hash table
 */
HT *ht_copy(HT *ht, void *(*copy_key)(void *)) {
	HT *ht2;
	int i;

	ht2 = ht_create(ht->length);

	for (i=0; i<ht->length; i++) {
		ht2->table[i] = dll_copy(ht->table[i], copy_key, NULL);
	}

	ht2->num_elem = ht->num_elem;

	return ht2;
}

/**
 * \brief Returns the linked list associated with a hash
 * 
 * \param[in] ht The hash table to lookup
 * \param[in] id The id (hash) to search for
 *
 * \return The DLL_NODE at the head of the dll containing all the elements
 * matching the search for hash.
 */
DLL_NODE *ht_hash_lookup(HT *ht, int id) {
	return ht->table[id % ht->length];
}

/**
 * \brief Finds a node in a hash table
 *
 * Finds a specific node in a hash table. Requires the knowledge of both the
 * index (hash) for the key, as well as the actual key. 
 *
 * \param[in] ht The hash table to be searched
 * \param[in] id The index (hash) of where the object should be
 * \param[in] key The key of the node that is being searched for
 *
 * \sa ht_find_node_fn
 *
 * \return The found node if found, otherwise NULL.
 */
DLL_NODE *ht_find_node(HT *ht, int id, void *key) {
	return dll_find_node(ht->table[id % ht->length], key);
}

/**
 * \brief Finds a node in a hash table using a comparator function
 *
 * Finds a specific node in a hash table. Requires the knowledge of both the
 * index (hash) for the key, as well as the actual key. 
 *
 * \param[in] ht The hash table to be searched
 * \param[in] id The index (hash) of where the object should be
 * \param[in] key The key of the node that is being searched for
 * \param[in] comparator A function that takes two keys and returns 1 if they
 * are equal and 0 otherwise.
 *
 * \sa ht_find_node
 *
 * \return The found node if found, otherwise NULL.
 */
DLL_NODE *ht_find_node_fn(HT *ht, int id, void *key, int (*comparator)(void *key1, void *key2)) {
	return dll_find_node_fn(ht->table[id % ht->length], key, comparator);
}

/**
 * \brief Deletes a key from a hash table
 *
 * Deletes a key from a hash table. Also removes the node the key exists in.
 * Will do nothing if the key does not exist in the hash table.
 *
 * \param[in] ht The hash table the key exists in
 * \param[in] id The index (hash) of where the object should be
 * \param[in] key The key of the node that is to be deleted
 * \param[in] free_key A function that takes a key and frees any memory
 * associated with it. Can be NULL to avoid freeing the key.
 */
void ht_delete_key(HT *ht, int id, void *key, void (*free_key)(void *key)) {
	int i;
	DLL_NODE *n;

	i = id % ht->length;
	n = dll_find_node(ht->table[i], key);

	if (n != NULL) {
		ht->table[i] = dll_delete_node(ht->table[i], n, free_key);
		ht->num_elem--;
	}
}

/**
 * \brief Prints a hash table
 *
 * \param[in] ht The hash table to be printed
 * \param[in] print_key A function that accepts the key from a hash table node that 
 * specifies how to print each individual node.
 */
void ht_print(HT *ht, void (*print_key)(void *key)) {
	int i;

	for (i = 0; i < ht->length; i++) {
		dll_print(ht->table[i], print_key);
	}
}

/**
 * \brief Frees a hash table
 *
 * \param[in] ht The hash table to be deleted
 * \param[in] free_key A function that takes a key and frees any memory
 * associated with it. Can be NULL to avoid freeing the key.
 */
void ht_free(HT *ht, void (*free_key)(void *key)) {
	int i;

	if (ht == NULL) return;

	for (i = 0; i < ht->length; i++) {
		dll_free(ht->table[i], free_key);
	}

	free(ht->table);
	free(ht);
}
