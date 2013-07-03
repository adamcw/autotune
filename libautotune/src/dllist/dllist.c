#include <stdio.h>
#include <stdlib.h>
#include "../memory/memory.h"
#include "dllist.h"

/**
 * \brief Creates a new double-linked list
 *
 * Exists to remain consistent with the llist and cdllist libraries.
 */
DLL_NODE *dll_create() {
	return NULL;
}

/**
 * \brief Inserts a new key into a linked list
 *
 * \param[in] head The head node of the linked list to be inserted into
 * \param[in] key A void pointer to the key to be inserted 
 * \param[in] set_key_node A function that takes the key given and the newly
 * created node and links the node to the key. Can be NULL if the key doesn't
 * need to know which node it is in.
 *
 * \return The head node of the updated linked list
 */
DLL_NODE *dll_insert(DLL_NODE *head, void *key, void (*set_key_node)(void *key, DLL_NODE *n)) {
	DLL_NODE *n;

	n = dll_create_node(key);
	if (set_key_node != NULL) {
		set_key_node(key, n);
	}
	head = dll_insert_node(head, n);
	
	return head;
}

/**
 * \brief Deletes a node from a linked list based on key
 *
 * \param[in] head The head node of the linked list to be deleted from
 * \param[in] key A void pointer to the key to be removed from the list. The
 * first instance of the key will be removed.
 * \param[in] get_key_node A function that takes a key and returns the DLL_NODE
 * associated with it.
 * \param[in] free_key A function that accepts the key from a linked list node
 * that specifies how to free any memory associated with the key. Can be NULL
 * if the key does not need to be freed.
 *
 * \sa dll_delete_head, dll_delete_node, dll_remove_node
 *
 * \return The head node of the updated linked list
 */
DLL_NODE *dll_delete(DLL_NODE *head, void *key, DLL_NODE *(*get_key_node)(void *key), void (*free_key)(void *)) {
	DLL_NODE *n;

	n = get_key_node(key);
	head = dll_delete_node(head, n, free_key);

	return head;
}

/** 
 * \brief Creates a new node for the linked list
 *
 * Creates a new node for the linked list, initialising the next and previous
 * pointers to NULL and setting the key as provided.
 *
 * \param[in] key A key to be stored in the node
 *
 * \return A newly allocated DLL_NODE
 */
DLL_NODE *dll_create_node(void *key) {
	DLL_NODE *n;

	n = (DLL_NODE *)my_malloc(sizeof(DLL_NODE));
	n->next = NULL;
	n->prev = NULL;
	n->key = key;

	return n;
}

/**
 * \brief Inserts a node into a linked list
 *
 * Inserts a node into a linked list. The node is inserted at the head of the
 * list.
 *
 * \param[in] head The head node of the linked list to be inserted into
 * \param[in] n The node that should be inserted into the list
 *
 * \return The updated head node of the linked list
 */
DLL_NODE *dll_insert_node(DLL_NODE *head, DLL_NODE *n) {
	if (head != NULL) {
		head->prev = n;
		n->next = head;
	}

	return n;
}


/**
 * \brief Clones a linked list
 *
 * \param[in] head The head node of a linked list to be copied
 * \param[in] copy_key A function that accepts the key from a linked list node that
 * creates and returns a copy of the key. 
 * \param[in] set_key_node A function that takes the key given and the newly
 * created node and links the node to the key. Can be NULL if the key doesn't
 * need to know which node it is in.
 *
 * \return The head node of the newly cloned linked list
 */
DLL_NODE *dll_copy(DLL_NODE *head, void *(*copy_key)(void *key), void (*set_key_node)(void *key, DLL_NODE *n)) {
    DLL_NODE *head2, *n, *n2, *nt;

    head2 = dll_create();

    if (head == NULL) return head2;

    head2 = dll_insert(head2, copy_key(head->key), set_key_node);
    n2 = head2;
    n = head->next;

    while (n != NULL) {
        n2->next = nt = (DLL_NODE *)my_malloc(sizeof(DLL_NODE));
        nt->prev = n2;
        n2 = nt;
        n2->key = copy_key(n->key);
        if (set_key_node != NULL) set_key_node(n2->key, n2);
        n = n->next;
    }

    n2->next = NULL;

    return head2;
}

/**
 * \brief Deletes the first node from a linked list 
 *
 * \param[in] head The head node of the linked list
 * \param[in] free_key A function that accepts the key from a linked list node
 * that specifies how to free any memory associated with the key. Can be NULL
 * if the key does not need to be freed.
 *
 * \sa dll_delete, dll_delete_node, dll_remove_node
 *
 * \return The updated head node of the linked list
 */
DLL_NODE *dll_delete_head(DLL_NODE *head, void (*free_key)(void *)) {
	head = dll_delete_node(head, head, free_key);
	return head;
}

/**
 * \brief Deletes an arbitarily positioned node from a linked list
 *
 * \param[in] head The head node of the linked list
 * \param[in] n The node to be deleted
 * \param[in] free_key A function that accepts the key from a linked list node
 * that specifies how to free any memory associated with the key. Can be NULL
 * if the key does not need to be freed.
 *
 * \sa dll_delete, dll_delete_head, dll_remove_node
 *
 * \return The updated head node of the linked list
 */
DLL_NODE *dll_delete_node(DLL_NODE *head, DLL_NODE *n, void (*free_key)(void *)) {
	if (n == NULL) return head;

	if (n->next != NULL) {
		n->next->prev = n->prev;
	}
	if (n->prev != NULL) {
		n->prev->next = n->next;
	}

	if (n == head) {
		head = n->next;
	}

	if (free_key != NULL) {
		free_key(n->key);
	}
	free(n);

	return head;
}

/**
 * \brief Removes a node from a linked list
 *
 * Removes a node from a linked list. Differs from deleting the node as the
 * memory for the node is not freed.
 *
 * \param[in] head The head node of the linked list
 * \param[in] n The node to be removed
 *
 * \sa dll_delete, dll_delete_head, dll_delete_node
 *  
 * \return The updated head node of the linked list
 */
DLL_NODE *dll_remove_node(DLL_NODE *head, DLL_NODE *n){
	if (n == NULL) return NULL;

	if (n->next != NULL) {
		n->next->prev = n->prev;
	}

	if (n->prev != NULL) {
		n->prev->next = n->next;
	}

	if (n == head) {
		head = n->next;
	}

	return head;
}

/**
 * \brief Finds a node in a linked list with the given key
 *
 * Searches linearly to find a node in a linked list with the given key.
 * Returns the first matching node.
 *
 * \param[in] head The head node of the linked list
 * \param[in] key A pointer to the key to be found
 *
 * \sa dll_find_node_fn
 *
 * \return The DLL_NODE if a match is found, NULL otherwise.
 */
DLL_NODE *dll_find_node(DLL_NODE *head, void *key){
	DLL_NODE *n;

	n = head;
	while (n != NULL) {
		if (n->key == key) {
			return n;
		}
		n = n->next;
	}

	return NULL;
}

/**
 * \brief Finds a node in a linked list with the given key, using a comparison
 * function to determine equality.
 *
 * \param[in] head The head node of the linked list
 * \param[in] key A pointer to the key to be found
 * \param[in] comparator A function that takes two keys and compares them
 * returning 1 if they are equal or 0 in they differ.
 *
 * \sa dll_find_node
 *
 * \returns The DLL_NODE if a match is found, NULL otherwise.
 */
DLL_NODE *dll_find_node_fn(DLL_NODE *head, void *key, int (*comparator)(void *key1, void *key2)){
	DLL_NODE *n;

	n = head;
	while (n != NULL) {
		if (comparator(key, n->key)) {
			return n;
		}
		n = n->next;
	}

	return NULL;
}

/**
 * \brief Prints a linked list
 *
 * \param[in] head The head node of the linked list
 * \param[in] print_key A function that accepts the key from a linked list node that 
 * specifies how to print each individual node.
 */
void dll_print(DLL_NODE *head, void (*print_key)(void *)) {
	DLL_NODE *n;

	n = head;
	while (n != NULL) {
		print_key(n->key);
		n = n->next;
		printf("\n");
	}
}

/**
 * \brief Frees a linked list
 *
 * \param[in] head The head node of the linked list
 * \param[in] free_key A function that accepts the key from a linked list node
 * that specifies how to free any memory associated with the key. Can be NULL
 * if the key does not need to be freed.
 */
void dll_free(DLL_NODE *head, void (*free_key)(void *)) {
	while (head != NULL) {
		head = dll_delete_head(head, free_key);
	}
}
