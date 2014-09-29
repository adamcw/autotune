#ifndef LLIST_H
#define LLIST_H

/** \defgroup ll Single-Linked List 
 * @{ */

/**
 * \brief A node for a single-linked list 
 */
typedef struct ll_node LL_NODE;
struct ll_node {
    /** A pointer to the next node in the linked list. NULL if there is no next node */
	struct ll_node *next; 

	/** A void pointer to an arbitary key to be stored in the node */
    void *key;
};

LL_NODE *ll_create();
LL_NODE *ll_insert(LL_NODE *head, void *key);
LL_NODE *ll_insert_sorted(LL_NODE *head, void *key, int (*key_lt)(void *k1, void *k2));
int ll_eq(LL_NODE *head1, LL_NODE *head2, int (*key_eq)(void *k1, void *k2));
LL_NODE *ll_copy(LL_NODE *head, void *(*copy_key)(void *));
LL_NODE *ll_delete(LL_NODE *head);
LL_NODE *ll_delete_node(LL_NODE *head, void *key);
void ll_print(LL_NODE *head, void (*print_key)(void *));
void ll_free(LL_NODE *head, void (*free_key)(void *));

/** @} */

#endif
