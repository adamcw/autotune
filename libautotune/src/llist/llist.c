#include <stdio.h>
#include <stdlib.h>
#include "llist.h"
#include "../memory/memory.h"

/** 
 * \brief Creates a new single-linked list
 *
 * Exists to remain consistent with the dllist and cdllist libraries.
 * \return NULL
 */
LL_NODE *ll_create() {
	return NULL;
}

/**
 * \brief Inserts a new key into a linked list
 *
 * \param[in] head The head node of the linked list to be inserted into
 * \param[in] key A void pointer to the key to be inserted
 *
 * \sa ll_insert_sorted
 *
 * \return The head node of the updated linked list
 */
LL_NODE *ll_insert(LL_NODE *head, void *key) {
	LL_NODE *n;

	n = (LL_NODE *)my_malloc(sizeof(LL_NODE));
	n->key = key;
	n->next = head;
	
	return n;
}

/**
 * \brief Inserts a new key into a sorted linked list
 *
 * \param[in] head The head node of the linked list to be inserted into
 * \param[in] key A void pointer to the key to be inserted
 * \param[in] key_lt A function that takes two keys from a linked list node and
 * compares them, returns 1 when k1 is less than k2, and 0 otherwise.
 *
 * \sa ll_insert
 *
 * \return The head node of the updated linked list
 */
LL_NODE *ll_insert_sorted(LL_NODE *head, void *key, int (*key_lt)(void *k1, void *k2)) {
	LL_NODE *n, *m, *prev;

	// Create the node for insertion
	m = (LL_NODE *)my_malloc(sizeof(LL_NODE));
	m->key = key;
	m->next = NULL;

	// Search for the correct position
	prev = NULL;
	n = head;
	while (n != NULL && !key_lt(m->key, n->key)) {
		prev = n;
		n = n->next;	
	}

	// Insert Node
	m->next = n;
	if (prev != NULL) {
		prev->next = m;
	} else {
		head = m;
	}

	return head;
}

/**
 * \brief Compares two linked lists for equality
 *
 * \param[in] head1 The head node of the first linked list to be compared
 * \param[in] head2 The head node of the second linked list to be compared
 * \param key_eq A function that accepts two keys from linked list nodes,
 * compares them, and returns whether they are equal.
 *
 * \return 1 if the linked lists are equal, 0 if they differ
 */
int ll_eq(LL_NODE *head1, LL_NODE *head2, int (*key_eq)(void *k1, void *k2)) {
	LL_NODE *n1, *n2;

	n1 = head1;
	n2 = head2;

	while (n1 != NULL && n2 != NULL) {
		if (!key_eq(n1->key, n2->key)) return 0;
		n1 = n1->next;
		n2 = n2->next;
	}

	if (n1 == n2) return 1;

	return 0;
}

/**
 * \brief Clones a linked list
 *
 * \param[in] head The head node of a linked list to be copied
 * \param[in] copy_key A function that accepts the key from a linked list node that
 * creates and returns a copy of the key. 
 *
 * \return The head node of the newly cloned linked list
 */
LL_NODE *ll_copy(LL_NODE *head, void *(*copy_key)(void *)) {
	LL_NODE *head2, *n, *n2;

	head2 = ll_create();

	if (head == NULL) return head2;

	head2 = ll_insert(head2, copy_key(head->key));
	n2 = head2;
	n = head->next;

	while (n != NULL) {
		n2->next = (LL_NODE *) my_malloc(sizeof(LL_NODE));
		n2 = n2->next;
		n2->key = copy_key(n->key);
		n = n->next;
	}

	n2->next = NULL;

	return head2;
}

/**
 * \brief Deletes the first node from a linked list
 *
 * \param[in] head The head node of the linked list
 *
 * \sa ll_delete_node
 *
 * \return The new head node of the linked list
 */
LL_NODE *ll_delete(LL_NODE *head) {
	LL_NODE *n;

	n = head->next;
	free(head);

	return n;
}

/**
 * \brief Delete a node in an arbitary position in a linked list.
 * Removing the first node that matches the provided key.
 *
 * \param[in] head The head node of the linked list
 * \param[in] key The value of the node you wish to remove
 *
 * \sa ll_delete
 *
 * \return The head node of the updated linked list
 */
LL_NODE *ll_delete_node(LL_NODE *head, void *key) {
    LL_NODE *n, *prev;

    prev = NULL;
    n = head;

    while (n != NULL) {
        if (n->key == key) {
            if (prev) {
                prev->next = n->next;
            } else {
                head = n->next;
            }
            free(n);
            break;
        }
        else {
            prev = n;
            n = n->next;
        }
    }

    return head;
}
/* More efficent version, but needs testing
 *
 * LL_NODE *ll_delete_node(LL_NODE *head, void *key) {
    LL_NODE *n, **prev;

    n = head;
    prev = &head;
    while (n != NULL) {
        if (n->key == key) {
            *prev = n->next;
            free(n);
            break;
        }
        prev = &n->next;
        n = n->next;
    }
    return head;
}
*/


/**
 * \brief Prints a linked list
 *
 * \param[in] head The head node of the linked list
 * \param[in] print_key A function that accepts the key from a linked list node that 
 * specifies how to print each individual node.
 */
void ll_print(LL_NODE *head, void (*print_key)(void *)) {
	LL_NODE *temp;

	temp = head;
	while (temp != NULL) {
		print_key(temp->key);
		temp = temp->next;
	}

	printf("\n");
}

/**
 * \brief Frees a linked list
 *
 * \param[in] head The head node of the linked list
 * \param[in] free_key A function that accepts the key from a linked list node
 * that specifies how to free any memory associated with the key. Can be NULL
 * if the key does not need to be freed.
 */
void ll_free(LL_NODE *head, void (*free_key)(void *)) {
	while (head != NULL) {
		if (free_key != NULL) free_key(head->key);
		head = ll_delete(head); 
	}
}
