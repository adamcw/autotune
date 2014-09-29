#include <stdio.h>
#include <assert.h>
#include "bfs.h"
#include "../memory/memory.h"
#include "../match/match.h"

BFS_GRAPH *bfs_create_graph() {
	BFS_GRAPH *g;

	g = (BFS_GRAPH *)my_malloc(sizeof(BFS_GRAPH));

	g->num_vertices = 0;
	g->vertex_ll = cdll_create();
	g->vertex_h = bh_create(0);
	g->num_edges = 0;
	g->edge_ll = cdll_create();

	return g;
}

void bfs_validate(BFS_GRAPH *g) {
	CDLL_NODE *cdlln;
	cdlln = g->vertex_ll->next;
	while (cdlln != g->vertex_ll) {
		BFS_VERTEX *v = (BFS_VERTEX *)cdlln->key;
		if (v->v_num > 0 && v->num_edges == 0) {
			printf("Vertex with no edges "); bfs_print_vertex(v);
			//assert(v->num_edges != 0);
		}
		cdlln = cdlln->next;
	}
}

void bfs_free_graph(BFS_GRAPH *g) {
	cdll_free(g->vertex_ll, bfs_free_void_vertex);
	bh_free(g->vertex_h, NULL);
	cdll_free(g->edge_ll, bfs_free_void_edge);
	free(g);
}

void bfs_remove_edge(BFS_GRAPH *g, BFS_EDGE *e) {
	e->a->num_edges--;
	e->b->num_edges--;
	g->num_edges--;
	
	cdll_delete_node(e->graph_n, NULL);
	cdll_delete_node(e->a_n, NULL);
	cdll_delete_node(e->b_n, NULL);
	bfs_free_void_edge(e);
}

void bfs_remove_vertex(BFS_GRAPH *g, BFS_VERTEX *v) {
	CDLL_NODE *n;
	BFS_EDGE *e;

	// Remove the vertex from the heap if it is in there.
	if (v->i > 0) {
		bh_remove_key(g->vertex_h, v->i, bfs_vertex_lt_d, bfs_vertex_swap, bfs_set_vertex_i, NULL);
	}

	// Loop and clean up the edges
	n = v->edge_ll->next;
	while (n != v->edge_ll) {
		e = (BFS_EDGE *)n->key;
		n = n->next;
		bfs_remove_edge(g, e);
	}
	
	cdll_delete_node(v->graph_n, bfs_free_void_vertex);
	g->num_vertices--;
}

void bfs_init_search(BFS_GRAPH *g, BFS_VERTEX *v0) {
	CDLL_NODE *n;
	BFS_VERTEX *v;

	n = g->vertex_ll->next;
	while (n != g->vertex_ll) {
		v = (BFS_VERTEX *)n->key;
		v->d = -1;
		v->i = 0;
		n = n->next;
	}

	bh_reset(g->vertex_h, NULL);

	v0->d = 0;
	bfs_insert_vertex_in_h(g, v0);
}

BFS_VERTEX *bfs_remove_first_v(BHEAP *h) {
	return (BFS_VERTEX *)bh_remove_key(h, 1, bfs_vertex_lt_d, bfs_vertex_swap, bfs_set_vertex_i, NULL);
}

BFS_VERTEX *bfs_get_next_v(BFS_GRAPH *g) {
	BFS_VERTEX *v, *b;
	BFS_EDGE *e;
	CDLL_NODE *n;

	if (g->vertex_h->num_elem == 0) return NULL;

	v = bfs_remove_first_v(g->vertex_h);

	DOT *vy = (DOT *)v->key;
	
	if (vy->t < 0 || vy->i < 0 || vy->j < 0) {
		return v;
	}

	n = v->edge_ll->next;
	while (n != v->edge_ll) {
		e = (BFS_EDGE *)n->key;
		b = e->a == v ? e->b : e->a;

		if (b->i > 0) {
			if (v->d + e->wt < b->d) {
				b->d = v->d + e->wt;
				bh_bubble_up_key(g->vertex_h, b->i, bfs_vertex_lt_d, bfs_vertex_swap);
			}
		}
		else if (b->d < 0) {
			b->d = v->d + e->wt;
			bfs_insert_vertex_in_h(g, b);
		}
		n = n->next;
	}

	return v;
}

BFS_VERTEX *bfs_create_vertex(int v_num, void *key) {
	BFS_VERTEX *v;

	v = (BFS_VERTEX *)my_malloc(sizeof(BFS_VERTEX));

	v->key = key;
	v->v_num = v_num;
	v->num_edges = 0;
	v->edge_ll = cdll_create();
	v->d = -1;
	v->i = 0;
	v->graph_n = NULL;

	return v;
}

void bfs_free_vertex(BFS_VERTEX *v) {
	cdll_free(v->edge_ll, NULL);
	free(v);
}

void bfs_free_void_vertex(void *key) {
	bfs_free_vertex((BFS_VERTEX *)key);
}

int bfs_vertex_lt_d(void *k1, void *k2) {
	BFS_VERTEX *v1, *v2;

	v1 = (BFS_VERTEX *)k1;
	v2 = (BFS_VERTEX *)k2;

	return v1->d < v2->d;
}

void bfs_vertex_swap(BHEAP *h, int i, int j) {
   int temp;
   void *ptr;
   BFS_VERTEX *v1, *v2;

   ptr = h->k[i];
   h->k[i] = h->k[j];
   h->k[j] = ptr;

   v1 = (BFS_VERTEX *)h->k[i];
   v2 = (BFS_VERTEX *)h->k[j];

   temp = v1->i;
   v1->i = v2->i;
   v2->i = temp;
}

void bfs_set_vertex_i(void *key, int i) {
	((BFS_VERTEX *)key)->i = i;
}

void bfs_set_vertex_graph_n(void *key, CDLL_NODE *n) {
	((BFS_VERTEX *)key)->graph_n = n;
}

void bfs_set_edge_graph_n(void *key, CDLL_NODE *n) {
	((BFS_EDGE *)key)->graph_n = n;
}

void bfs_set_edge_a_n(void *key, CDLL_NODE *n) {
	((BFS_EDGE *)key)->a_n = n;
}

void bfs_set_edge_b_n(void *key, CDLL_NODE *n) {
	((BFS_EDGE *)key)->b_n = n;
}

void bfs_insert_vertex(BFS_GRAPH *g, BFS_VERTEX *v) {
	g->num_vertices++;
	cdll_insert_head(g->vertex_ll, v, bfs_set_vertex_graph_n);
}

void bfs_insert_vertex_in_h(BFS_GRAPH *g, BFS_VERTEX *v) {
	bh_insert_key(g->vertex_h, v, bfs_vertex_lt_d, bfs_vertex_swap, bfs_set_vertex_i, NULL);
}

BFS_EDGE *bfs_create_edge(BFS_VERTEX *a, BFS_VERTEX *b, int wt) {
	BFS_EDGE *e;

	e = (BFS_EDGE *)my_malloc(sizeof(BFS_EDGE));

	e->a = a;
	e->b = b;
	e->wt = wt;
	e->graph_n = NULL;
	e->a_n = NULL;
	e->b_n = NULL;

	return e;
}

void bfs_free_edge(BFS_EDGE *e) {
	free(e);
}

void bfs_free_void_edge(void *key) {
	bfs_free_edge((BFS_EDGE *)key);
}

void bfs_insert_edge(BFS_GRAPH *g, BFS_EDGE *e) {
	g->num_edges++;
	e->a->num_edges++;
	e->b->num_edges++;
	cdll_insert_head(g->edge_ll, e, bfs_set_edge_graph_n);
	cdll_insert_head(e->a->edge_ll, e, bfs_set_edge_a_n);
	cdll_insert_head(e->b->edge_ll, e, bfs_set_edge_b_n);
}

void bfs_print_graph(BFS_GRAPH *g) {
	int i;
	CDLL_NODE *n;

	printf("num_edges: %d\n\n", g->num_edges);

	n = g->edge_ll->next;
	for (i=0; i<g->num_edges; i++) {
		printf("%3d: ", i);
		bfs_print_edge((BFS_EDGE *)n->key);
		n = n->next;
	}

	assert(n == g->edge_ll);
}

void bfs_print_vertex(BFS_VERTEX *v) {
	int i;
	CDLL_NODE *n;
	printf("v: %d, d: %d, i: %d, ne: %d\n", v->v_num, v->d, v->i, v->num_edges);

	n = v->edge_ll->next;
	for (i = 0; i < v->num_edges; i++) {
		printf("%3d: ", i);
		bfs_print_edge((BFS_EDGE *)n->key);
		n = n->next;
	}
}

void bfs_print_void_vertex(void *key) {
	bfs_print_vertex((BFS_VERTEX *)key);
}

void bfs_print_vertex_h(BFS_GRAPH *g) {
	bh_print(g->vertex_h, bfs_print_void_vertex);
}

void bfs_print_edge(BFS_EDGE *e) {
	printf("%3d -- %3d, %d\n", e->a->v_num, e->b->v_num, e->wt);
}

