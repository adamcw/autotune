#ifndef DLLIST_H
#define DLLIST_H

/** \defgroup dll Double-Linked List 
 * @{ */

/**
 * \brief A node for a double-linked list
 */
typedef struct dll_node DLL_NODE;
struct dll_node {
	/** A pointer to the next node in the linked list. NULL if there is no next node */
	DLL_NODE *next;
	
    /** A pointer to the previous node in the linked list. NULL if there is no previous node */
	DLL_NODE *prev;
	
    /** A void pointer to an arbitary key to be stored in the node */
    void *key;
};

DLL_NODE *dll_create();
DLL_NODE *dll_insert(DLL_NODE *head, void *key, void (*set_key_node)(void *key, DLL_NODE *n));
DLL_NODE *dll_delete(DLL_NODE *head, void *key, DLL_NODE *(*get_key_node)(void *key), void (*free_key)(void *));
DLL_NODE *dll_create_node(void *key);
DLL_NODE *dll_insert_node(DLL_NODE *head, DLL_NODE *n);
DLL_NODE *dll_copy(DLL_NODE *head, void *(*copy_key)(void *key), void (*set_key_node)(void *key, DLL_NODE *n));
DLL_NODE *dll_delete_head(DLL_NODE *head, void (*free_key)(void *));
DLL_NODE *dll_delete_node(DLL_NODE *head, DLL_NODE *n, void (*free_key)(void *));
DLL_NODE *dll_remove_node(DLL_NODE *head, DLL_NODE *n);
DLL_NODE *dll_find_node(DLL_NODE *head, void *key);
DLL_NODE *dll_find_node_fn(DLL_NODE *head, void *key, int (*comparator)(void *key1, void *key2));
void dll_print(DLL_NODE *head, void (*print_key)(void *));
void dll_free(DLL_NODE *head, void (*free_key)(void *));

/** @} */

#endif
