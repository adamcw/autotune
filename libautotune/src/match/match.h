#ifndef __included_match_h__
#define __included_match_h__

#include "../qc/qc.h"
#include "../bfs/bfs.h"
#include "../../../blossomv/PerfectMatching.h"

/* undo definitions */

#define UNDO_CREATE_VERTEX 0
#define UNDO_CREATE_AND_INSERT_DOT 1
#define UNDO_CREATE_AND_INSERT_LINE 2
#define UNDO_UPDATE_T_BFS 3
#define UNDO_DELETE_EDGE 4
#define UNDO_MATCHED_EDGE 5
#define UNDO_UPDATE_LAST_GRAPH_HEAD 6
#define UNDO_UPDATE_LAST_GRAPH_TAIL 7
#define UNDO_CREATE_LINE 8
#define UNDO_INSERT_LINE 9
#define UNDO_CREATE_DOT 10
#define UNDO_INSERT_DOT 11
#define UNDO_LINE_MOVE 12
#define UNDO_LINE_EDIT 13
#define UNDO_LINE_WT 14
#define UNDO_DOT_MERGE 15
#define UNDO_DOT_T 16

/* general constants */

#define PRECISION 1000

#define TRUE 1
#define FALSE 0

typedef struct blossomv BLOSSOMV;
typedef struct edge EDGE;
typedef struct dot DOT;
typedef struct line LINE;
typedef struct aug_edge AUG_EDGE;
typedef struct vertex VERTEX;
typedef struct matching MATCHING;
typedef struct undo UNDO;

/**
 * \brief Stores the state of the matching problem to be passed into the
 * BlossomV module. 
 *
 * BlossomV accepts the problem in the form of two arrays, edges and weights.
 * Edges is a list where each pair of consecutive integers are the source and
 * destination ids of the edge. The weights array then forms a list of the
 * weights of each corresponding edge. Vertex IDs are never passed in, but must
 * start at 0, leading the vertices have IDs in graph space and in BlossomV
 * space.
 */
struct blossomv {
	/** The number of vertices to be matched */
	int num_vertices;

	/** The number of edges in the problem */
	int num_edges;

	/** The number of edges allocated in memory */
	int alloc_edges;
	
	/** An array of all the edges */
	int *edges;

	/** An array of the edge weights */
	int *weights;

	/** An array of all the vertices */
	VERTEX **vertices;
};

/* structure definitions */

struct edge {
	/** The source vertex of the edge */
	VERTEX *source_vertex;

	/** The destination vertex of the edge */
	VERTEX *dest_vertex;
	
	/** The reverse complement of this edge */
	EDGE *backward_edge;
};

struct dot {
	/** The i corrdinate of the dot */
	int i;
	
	/** The j coordinate of the dot */
	int j;

	/*
	 * The big_t coordinate of the dot. 
	 *
	 * This confusing syntax is due to match
	 * predating a lot of other code 
	 */
	long int t;

	/** The t coordinate of the dot */
	long int little_t;

	/** The vertex associated with the dot */
	VERTEX *v;

	/** The lines eminating from the dot */
	LL_NODE *lines;

	/** The ball associated with the dot */
	BALL *ball;

	/** The vertex in the bfs structure for this dot */
	BFS_VERTEX *bfs;

	DOT *merge;
};

struct line {
	/** The source dot of the line */
	DOT *a;

	/** The destination dot of the line */
	DOT *b;

	/** The weight of the line */
	int wt;
};

struct aug_edge {
	/** The uplink to the list of all augmented edges */
	CDLL_NODE *n;
	
	/** The source vertex of the edge */
	VERTEX *va;

	/** The destination vertex of the edge */
	VERTEX *vb;
};

struct vertex {
	/** The ID of the vertex in graph space */
	int v_num;
	
	/** The i coordinate of the vertex */
	int i;

	/** The j coordinate of the vertex */
	int j;

	/** The t coordinate of the vertex */
	long int t;

	/** The associated dot with the vertex */
	DOT *dot;

	/** The edge the vertex is matched to, NULL if no edge */
	EDGE *matched_edge; 
	
	/** The uplink from the vertex to its position in the graph */
	CDLL_NODE *graph_node;
};

struct matching {
	/** The maximum distance to search for edges when performing BFS */
	int D;

	/** The PerfectMatching object from BlossomV module */
	PerfectMatching *pm;

	/** The vertices in the graph */
	CDLL_NODE *graph;

	/** The last vertex that was at the head of the graph */
	CDLL_NODE *last_graph_head;
	
	/** The last vertex that was at the tail of the graph */
	CDLL_NODE *last_graph_tail;

	/** Whether undo is enabled on this matching */
	int undo_flag;

	/** How many undos are in the list */
	int num_undos;

	/** The undo list */
	UNDO *u; 

	/** The \ref dot%s in the matching lattice */
	CDLL_NODE *dots;

	/** A hash table of the dots in the lattice */
	HT *dot_ht;

	/** The \ref line%s in the matching lattice */
	CDLL_NODE *lines;

	/** A hash table of the lines in the lattice */

	//HT *line_ht;

	/** The graph to perform breadth first search to connect vertices to only
	 * their nearest neighbors */
	BFS_GRAPH *g;

	/** The last timestep the \ref bfs_graph was updated */
	long int t_bfs;

	CDLL_NODE *t_cdlln;

	/** The number of \ref dot%s */
	int num_dots;
	
	/** The number of \ref line%s */
	int num_lines;

	/** The number of \ref vertex */
	int num_vertices;

	/** The list of augmented edges */
	CDLL_NODE *augmented_edges;

	/** The current timestep in the matching */
	long int t;
	
	/** The number of timesteps to keep when deleting */
	int t_delete;

	int t_delay;
};

struct undo {
	/** The type of undo */
	int type;

	/** A vertex to be saved */
	VERTEX *v;

	/** An increment to be saved */
	int inc;

	/** An edge to be saved */
	EDGE *e;

	/** A second vertex to be saved */
	VERTEX *v2;

	/** A ::CDLL_NODE to be saved */
	CDLL_NODE *n;

	DOT *d1;
	DOT *d2;
	LINE *line;

	long int t;

	/** The next undo in the list */
	UNDO *next;
};

// MATCHING

MATCHING *m_create_matching(int t_delete);
void m_free_matching(MATCHING *matching);
void m_free_undo(MATCHING *matching);

void m_update_dots_and_lines_into_bfs(MATCHING *matching, BFS_GRAPH *g, int undo); 
void m_free_dots_and_lines_bfs(BFS_GRAPH *g);

BLOSSOMV *m_create_blossomv(int num_vertices);
void m_free_blossomv(BLOSSOMV *bv);

int m_add_vertex_to_blossomv(BLOSSOMV *bv, VERTEX *v, int pos);
void m_add_edge_to_blossomv(BLOSSOMV *bv, int a, int b, int wt);
void m_create_augmented_edge(MATCHING *m, VERTEX *v1, VERTEX *v2);
DOT *m_create_temporal_boundary();
void m_create_and_insert_vertex_and_edge_into_bfs(BFS_GRAPH *g, DOT *a, DOT *b, int wt);
void m_solve(MATCHING *matching, BLOSSOMV *bv);
void m_create_augmented_edges(MATCHING *matching, BLOSSOMV *bv, int undo);

void m_mwpm(MATCHING *matching, int undo);
int m_time_delete(MATCHING *m); 
AUG_EDGE *m_get_aug_edge(MATCHING *m);
void m_delete_aug_edge(AUG_EDGE *ae);

// DOT

DOT *m_create_dot(BALL *ball, VERTEX *v);
void m_insert_dot(MATCHING *m, DOT *dot);
DOT *m_create_and_insert_dot(MATCHING *m, BALL *ball); 
DOT *m_create_and_insert_dot_and_vertex(MATCHING *m, BALL *ball);
void m_free_dot(DOT *dot);
void m_free_void_dot(void *key);

// LINE

LINE *m_create_and_insert_line(MATCHING *m, DOT *da, DOT *db, float p_stick); 
LINE *m_create_and_insert_line_wt(MATCHING *m, DOT *da, DOT *db, int wt);
LINE *m_create_line_wt(MATCHING *m, DOT *da, DOT *db, int wt); 
void m_insert_line(MATCHING *m, LINE *line);

// VERTEX

bool m_is_boundary(VERTEX *v);
VERTEX *m_create_vertex(int i, int j, long int t, DOT *dot);
void m_set_vertex_graph_node(void *key, CDLL_NODE *graph_node);
void m_insert_vertex(MATCHING *matching, VERTEX *v);
void m_delete_vertex(MATCHING *matching, VERTEX *v);
void m_fast_delete_vertex(VERTEX *v);
void m_fast_delete_void_vertex(void *key);

// EDGE

EDGE *m_create_only_edge(VERTEX *v1, VERTEX *v2);
void m_delete_only_edge(EDGE *e);

// UNDO

void m_create_undo(MATCHING *matching, int type, VERTEX *v, long int inc, EDGE *e, VERTEX *v2, CDLL_NODE *n);
void m_create_line_undo(MATCHING *matching, int type, LINE *line, DOT *d1, DOT *d2, int wt, long int t);
void m_execute_and_delete_undos(MATCHING *matching);

// PRINT FUNCTIONS

void m_print_matching(MATCHING *matching);
void m_print_lattice(MATCHING *matching);
void m_print_line(LINE *line);	
void m_print_void_line(void *key);
void m_print_dot(DOT *dot);
void m_print_void_dot(void *key);
void m_print_graph(MATCHING *matching);
void m_print_void_vertex(void *k);
void m_print_augmented_edges(MATCHING *matching);
void m_print_void_edge(void *ev);

#endif
