#ifndef CDLLIST_H
#define CDLLIST_H

/** \defgroup cdll Circular Double-Linked List 
 * @{ */

/**
 * \brief A node for a circular double-linked list
 */
typedef struct cdll_node CDLL_NODE;
struct cdll_node {
	/** A pointer to the next node in the linked list. NULL if there is no next node */
	CDLL_NODE *next;
	
    /** A pointer to the previous node in the linked list. NULL if there is no previous node */
	CDLL_NODE *prev;
	
    /** A void pointer to an arbitary key to be stored in the node */
    void *key;
};

CDLL_NODE *cdll_create();
CDLL_NODE *cdll_create_node(void *key);
void cdll_insert_head(CDLL_NODE *sent, void *key, void (*set_key_node)(void *key, CDLL_NODE *n));
void cdll_insert_tail(CDLL_NODE *sent, void *key, void (*set_key_node)(void *key, CDLL_NODE *n));
void cdll_insert_node_head(CDLL_NODE *sent, CDLL_NODE *n);
void cdll_insert_node_tail(CDLL_NODE *sent, CDLL_NODE *n);
CDLL_NODE *cdll_copy(CDLL_NODE *sent, void *(*copy_key)(void *), void (*set_key_node)(void *key, CDLL_NODE *n));
void cdll_delete(void *key, CDLL_NODE *(*get_key_node)(void *key), void (*free_key)(void *));
void cdll_delete_head(CDLL_NODE *sent, void (*free_key)(void *));
void cdll_delete_tail(CDLL_NODE *sent, void (*free_key)(void *));
void cdll_delete_node(CDLL_NODE *n, void (*free_key)(void *));
void cdll_join(CDLL_NODE *n1, CDLL_NODE *n2);
void cdll_extract_node(CDLL_NODE *n);
void cdll_print(CDLL_NODE *sent, void (*print_key)(void *));
void cdll_reverse_print(CDLL_NODE *sent, void (*print_key)(void *));
void cdll_free(CDLL_NODE *sent, void (*free_key)(void *));
int cdll_length(CDLL_NODE *sent);

/** @} */

#endif
