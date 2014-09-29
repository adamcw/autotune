#include <stdio.h>
#include <stdlib.h>
#include "cdllist.h"
#include "../memory/memory.h"

/**
 * \brief Creates a new circular double-linked list
 *
 * Allocates a new head node known as the "sentinel". Next and previous are
 * initialised to point to this sentinel. 
 *
 * \returns The sentinel for the new linked list
 */
CDLL_NODE *cdll_create() {
	CDLL_NODE *sent;

	sent = (CDLL_NODE *) my_malloc(sizeof(CDLL_NODE));
	sent->next = sent;
	sent->prev = sent;
	sent->key = NULL;

	return sent;
}

/**
 * \brief Creates a new node for a linked list 
 *
 * Creates a new node for a linked list, but does not insert it into any list.
 *
 * \param[in] key The key to be stored in the linked list. 
 *
 * \return A new CDLL_NODE with next and previous initialised to NULL.
 */
CDLL_NODE *cdll_create_node(void *key) {
	CDLL_NODE *n;

	n = (CDLL_NODE *) my_malloc(sizeof(CDLL_NODE));
	n->next = NULL;
	n->prev = NULL;
	n->key = key;

	return n;
}

/**
 * \brief Creates and inserts a node at the head of the linked list
 *
 * \param[in,out] sent The sentinel of a linked list
 * \param[in] key A key to be stored in the new node
 * \param[in] set_key_node A function that takes the key given and the newly
 * created node and links the node to the key. Can be NULL if the key doesn't
 * need to know which node it is in.
 *
 * \sa cdll_insert_tail, cdll_insert_node_head, cdll_insert_node_tail
 */
void cdll_insert_head(CDLL_NODE *sent, void *key, void(*set_key_node)(void *key, CDLL_NODE *n)) {
	CDLL_NODE *n;

	n = cdll_create_node(key);
	if (set_key_node != NULL) set_key_node(key, n);
	cdll_insert_node_head(sent, n);
}

/**
 * \brief Creates and inserts a node at the tail of the linked list
 *
 * \param[in,out] sent The sentinel of a linked list
 * \param[in] key A key to be stored in the new node
 * \param[in] set_key_node A function that takes the key given and the newly
 * created node and links the node to the key. Can be NULL if the key doesn't
 * need to know which node it is in.
 *
 * \sa cdll_insert_head, cdll_insert_node_head, cdll_insert_node_tail
 */
void cdll_insert_tail(CDLL_NODE *sent, void *key, void(*set_key_node)(void *key, CDLL_NODE *n)) {
	CDLL_NODE *n;

	n = cdll_create_node(key);
	if (set_key_node != NULL) set_key_node(key, n);
	cdll_insert_node_tail(sent, n);
}

/**
 * \brief Inserts a node at the head of the linked list
 *
 * \param[in,out] sent The sentinel of a linked list
 * \param[in] n The node to be inserted
 *
 * \sa cdll_insert_head, cdll_insert_tail, cdll_insert_node_tail
 */
void cdll_insert_node_head(CDLL_NODE *sent, CDLL_NODE *n) {
	n->next = sent->next;
	n->prev = sent;
	sent->next = n;
	n->next->prev = n;
}

/**
 * \brief Inserts a node at the tail of the linked list
 *
 * \param[in,out] sent The sentinel of a linked list
 * \param[in] n The node to be inserted
 *
 * \sa cdll_insert_head, cdll_insert_tail, cdll_insert_node_head
 */
void cdll_insert_node_tail(CDLL_NODE *sent, CDLL_NODE *n) {
	n->next = sent;
	n->prev = sent->prev;
	sent->prev = n;
	n->prev->next = n;
}

/**
 * \brief Deletes a node from a linked list based on key
 *
 * \param[in] key The key to be deleted from the linked list
 * \param[in] get_key_node A function that takes a key and returns the node 
 * the given key belongs to.
 * \param[in] free_key A function that accepts the key from a linked list node
 * that specifies how to free any memory associated with the key. Can be NULL
 * if the key does not need to be freed.
 * *
 * \sa cdll_delete_head, cdll_delete_tail, cdll_delete_node
 */
void cdll_delete(void *key, CDLL_NODE *(*get_key_node)(void *key), void (*free_key)(void *)) {
	CDLL_NODE *n;

	n = get_key_node(key);
	cdll_delete_node(n, free_key);
}

/**
 * \brief Deletes the node at the head of a linked list
 *
 * \param[in,out] sent The sentinel of a linked list
 * \param[in] free_key A function that accepts the key from a linked list node
 * that specifies how to free any memory associated with the key. Can be NULL
 * if the key does not need to be freed.
 * *
 * \sa cdll_delete, cdll_delete_tail, cdll_delete_node
 */
void cdll_delete_head(CDLL_NODE *sent, void (*free_key)(void *)) {
	cdll_delete_node(sent->next, free_key);
}

/**
 * \brief Deletes the node at the tail of a linked list
 *
 * \param[in,out] sent The sentinel of a linked list
 * \param[in] free_key A function that accepts the key from a linked list node
 * that specifies how to free any memory associated with the key. Can be NULL
 * if the key does not need to be freed.
 * *
 * \sa cdll_delete, cdll_delete_head, cdll_delete_node
 */
void cdll_delete_tail(CDLL_NODE *sent, void (*free_key)(void *)) {
	cdll_delete_node(sent->prev, free_key);
}

/**
 * \brief Deletes the node positions arbitrarily in a linked list 
 *
 * \param[in] n The node to be deleted
 * \param[in] free_key A function that accepts the key from a linked list node
 * that specifies how to free any memory associated with the key. Can be NULL
 * if the key does not need to be freed.
 * *
 * \sa cdll_delete, cdll_delete_head, cdll_delete_tail
 */
void cdll_delete_node(CDLL_NODE *n, void (*free_key)(void *)) {
	n->next->prev = n->prev;
	n->prev->next = n->next;

	if (free_key != NULL) {
		free_key(n->key);
	}
	free(n);
}

/**
 * \brief Joins two linked lists together
 *
 * Appends the start of n2 to the start of n1 and frees the sentinel used for
 * n2. This function is destructive to n2.
 *
 * \param[in,out] n1 The main linked list which the other will be joined to
 * \param[in] n2 The secondary linked list that will be attached to the first
 */
void cdll_join(CDLL_NODE *n1, CDLL_NODE *n2) {
	n1->prev->next = n2->next;
	n2->next->prev = n1->prev;
	n2->prev->next = n1;
	n1->prev = n2->prev;

	free(n2);
}

/**
 * \brief Extracts a node from a linked list
 *
 * Extracts a node from a linked list by removing itself from the linked list.
 * The node will still point into the list, the list will not point to the
 * node.
 *
 * \param[in,out] n The node to be extracted
 */
void cdll_extract_node(CDLL_NODE *n) {
	n->next->prev = n->prev;
	n->prev->next = n->next;
}

/**
 * \brief Prints a linked list
 *
 * \param[in] sent The sentinel of the linked list
 * \param[in] print_key A function that accepts the key from a linked list node that 
 * specifies how to print each individual node.
 */
void cdll_print(CDLL_NODE *sent, void(*print_key)(void *)) {
	CDLL_NODE *node;

	node = sent->next;
	while (node != sent) {
		if (node->key != NULL) {
			print_key(node->key);
		}
		node = node->next;
	}
}

/**
 * \brief Prints a linked list in reverse
 *
 * \param[in] sent The sentinel of the linked list
 * \param[in] print_key A function that accepts the key from a linked list node that 
 * specifies how to print each individual node.
 */
void cdll_reverse_print(CDLL_NODE *sent, void (*print_key)(void *)) {
	CDLL_NODE *node;

	node = sent->prev;
	while (node != sent) {
		if (node->key != NULL) {
			print_key(node->key);
		}
		node = node->prev;
	}
}

/**
 * \brief Frees a linked list
 *
 * \param[in] sent The sentinel of the linked list
 * \param[in] free_key A function that accepts the key from a linked list node
 * that specifies how to free any memory associated with the key. Can be NULL
 * if the key does not need to be freed.
 */
void cdll_free(CDLL_NODE *sent, void (*free_key)(void *)) {
	CDLL_NODE *n, *nt;

	if (sent == NULL) return;

	n = sent->next;
	while (n != sent) {
		nt = n;
		n = n->next;
		if (free_key != NULL) free_key(nt->key);
		free(nt);
	}

	free(sent);
}

/**
 * \brief Clones a linked list
 *
 * \param[in] sent The sentinel of the linked list to be copied
 * \param[in] copy_key A function that accepts the key from a linked list node that
 * creates and returns a copy of the key. 
 * \param[in] set_key_node A function that takes the key given and the newly
 * created node and links the node to the key. Can be NULL if the key doesn't
 * need to know which node it is in.
 *
 * \return The sentinel node of the newly cloned linked list
 */
CDLL_NODE *cdll_copy(CDLL_NODE *sent, void *(*copy_key)(void *), void (*set_key_node)(void *key, CDLL_NODE *n)) {
	CDLL_NODE *sent2, *n;

	sent2 = cdll_create();

	n = sent->prev;
	while (n != sent) {
		cdll_insert_head(sent2, copy_key(n->key), set_key_node);
		n = n->prev;
	}

	return sent2;
}

int cdll_length(CDLL_NODE *sent) {
	CDLL_NODE *n;
	int count;

	count = 0;
	n = sent->prev;
	while (n != sent) {
		count++;
		n = n->prev;
	}

	return count;
}
