#ifndef BFS_H
#define BFS_H

#include "../cdllist/cdllist.h"
#include "../bheap/bheap.h"

typedef struct bfs_graph BFS_GRAPH;
typedef struct bfs_vertex BFS_VERTEX;
typedef struct bfs_edge BFS_EDGE;

struct bfs_graph {
	int num_vertices;
	CDLL_NODE *vertex_ll;
	BHEAP *vertex_h;

	int num_edges;
	CDLL_NODE *edge_ll;
};

struct bfs_vertex {
	CDLL_NODE *graph_n;
	
	void *key;

	int v_num;
	int num_edges;
	CDLL_NODE *edge_ll;

	int d;
	int i;
};

struct bfs_edge {
	CDLL_NODE *graph_n, *a_n, *b_n;
	
	BFS_VERTEX *a, *b;

	int wt;
};

BFS_GRAPH *bfs_create_graph();
void bfs_validate(BFS_GRAPH *g);
void bfs_free_graph(BFS_GRAPH *g);
void bfs_remove_edge(BFS_GRAPH *g, BFS_EDGE *e);
void bfs_remove_vertex(BFS_GRAPH *g, BFS_VERTEX *v);

void bfs_set_vertex_graph_n(void *key, CDLL_NODE *n);
void bfs_set_edge_graph_n(void *key, CDLL_NODE *n);

BFS_GRAPH *bfs_create_random_graph(int nv, int ne, int nwt);
void bfs_init_search(BFS_GRAPH *g, BFS_VERTEX *v0);
BFS_VERTEX *bfs_remove_first_v(BHEAP *h);
BFS_VERTEX *bfs_get_next_v(BFS_GRAPH *g);

BFS_VERTEX *bfs_create_vertex(int v_num, void *key);
void bfs_free_vertex(BFS_VERTEX *v);
void bfs_free_void_vertex(void *key);
int bfs_vertex_lt_d(void *k1, void *k2);
void bfs_vertex_swap(BHEAP *h, int i, int j);
void bfs_set_vertex_i(void *key, int i);
void bfs_insert_vertex(BFS_GRAPH *g, BFS_VERTEX *v);
void bfs_insert_vertex_in_h(BFS_GRAPH *g, BFS_VERTEX *v);

BFS_EDGE *bfs_create_edge(BFS_VERTEX *a, BFS_VERTEX *b, int wt);
void bfs_free_edge(BFS_EDGE *e);
void bfs_free_void_edge(void *key);
void bfs_insert_edge(BFS_GRAPH *g, BFS_EDGE *e);

void bfs_print_graph(BFS_GRAPH *g);
void bfs_print_vertex(BFS_VERTEX *v);
void bfs_print_void_vertex(void *key);
void bfs_print_vertex_h(BFS_GRAPH *g);
void bfs_print_edge(BFS_EDGE *e);

#endif

