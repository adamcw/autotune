#include <assert.h>
#include <math.h>
#include "match.h"
#include "../memory/memory.h"
#include "../bfs/bfs.h"
#include "../fast_hash/fasthash.h"

#define MATCHING_DOT_HASH_SIZE 15485863

/**
 * \brief Creates a new \ref matching 
 * 
 * \param[in] t_delete The number of timesteps to keep in memory for matching
 *
 * \return The newly created \ref matching
 */
MATCHING *m_create_matching(int t_delete) {
	MATCHING *matching;

	matching = (MATCHING *)my_malloc(sizeof(MATCHING));

	// min_horz_wt * d
	matching->D = INT_MAX;

	// Initialise the BlossomV Perfect Matching
	matching->pm = NULL;

	matching->last_graph_head = NULL;
	matching->last_graph_tail = NULL;

	matching->graph = cdll_create();
	matching->undo_flag = FALSE;
	matching->num_undos = 0;
	matching->u = NULL;

	matching->g = bfs_create_graph();
	matching->t_bfs = -1;

	// lattice
	matching->dots = cdll_create();
	matching->lines = cdll_create();
	matching->dot_ht = ht_create(MATCHING_DOT_HASH_SIZE);
	//matching->line_ht = ht_create(MATCHING_DOT_HASH_SIZE);

	// lattice and graph statistics
	matching->num_dots = 0;
	matching->num_lines = 0;
	matching->num_vertices = 0;

	matching->augmented_edges = cdll_create();
	
	// time variables
	matching->t = 0;
	matching->t_delete = t_delete;
	matching->t_delay = 0;

	return matching;
}

/**
 * \brief Frees a \ref matching
 * 
 * \param[in] matching The \ref matching to be freed
 */
void m_free_matching(MATCHING *matching) {
	CDLL_NODE *n;
	BFS_VERTEX *bfsv1;
	DOT *a;

	if (matching == NULL) {
		return;
	}
 
	// Free any straggling temporal boundaries
	n = matching->g->vertex_ll->next;
	while (n != matching->g->vertex_ll) {
		bfsv1 = (BFS_VERTEX *)n->key;
		n = n->next;
		// We only care about the vertices that are temporal boundaries
		a = (DOT *)bfsv1->key;
		if (a->v && a->v->t == -INT_MAX) {
			m_delete_vertex(matching, a->v);
			m_free_dot(a);
		}
	}

	// Free any straggling temporal boundaries
	n = matching->g->vertex_ll->next;
	while (n != matching->g->vertex_ll) {
		bfsv1 = (BFS_VERTEX *)n->key;
		n = n->next;
		// We only care about the vertices that are temporal boundaries
		a = (DOT *)bfsv1->key;
		if (a->v && a->v->t == -INT_MAX) {
			m_delete_vertex(matching, a->v);
			m_free_dot(a);
		}
	}

	cdll_free(matching->graph, m_fast_delete_void_vertex);
	m_free_undo(matching);

	cdll_free(matching->dots, m_free_void_dot);
	cdll_free(matching->lines, free);

	ht_free(matching->dot_ht, NULL);
	//ht_free(matching->line_ht, NULL);
	
	bfs_free_graph(matching->g);

	cdll_free(matching->augmented_edges, free);

	free(matching);
}

/**
 * \brief Frees the list of \ref undo
 * 
 * \param[in] matching The \ref matching from which to free all \ref undo\s 
 */
void m_free_undo(MATCHING *matching) {
	UNDO *ut;
	while (matching->u != NULL) {
		ut = matching->u;
		matching->u = matching->u->next;
		free(ut);
	}
}

/**
 * \brief Updates the \ref BFS structure to include the latest \ref dot%s and \ref line%s
 * 
 * \param[in] matching The \ref matching containing the \ref dots and \ref lines
 * \param[out] graph The \ref BFS graph to update
 * \param[in] undo Whether or not this update should be undoable  
 */
void m_update_dots_and_lines_into_bfs(MATCHING *matching, BFS_GRAPH *g, int undo) {
	int v_num;
	DOT *a, *b;
	LINE *line;
	CDLL_NODE *node;
	BFS_EDGE *e;
	BFS_VERTEX *v;
	long int latest_t;

	// Insert the dots into the BFS
	v_num = g->num_vertices;
	node = matching->dots->next;
	latest_t = matching->t_bfs;
	while (node != matching->dots) {
		a = (DOT *)node->key;

		//printf("%d, %d, %ld, %p | %ld, %ld\n", a->i, a->j, a->t, a, matching->t, matching->t_bfs);

		if (a->t < matching->t_bfs) {
			break;
		}

		if (a->t > latest_t) {
			latest_t = a->t;
		}

		// If the dot already has a BFS associated, then it's likely that the
		// BFS wasn't cleared properly from previous rounds.
		if (a->bfs != NULL) {
			node = node->next;
			//printf("--%d, %d, %ld | %ld, %ld\n", a->i, a->j, a->t, matching->t, matching->t_bfs);
			continue;

			bfs_remove_vertex(g, a->bfs);
			//printf("Dot already allocated to BFS.\n");
			//assert(a->bfs == NULL);
		}
		
		//printf("+-+ %d, %d, %ld, %p | %ld, %ld\n", a->i, a->j, a->t, a, matching->t, matching->t_bfs);

		v = bfs_create_vertex(v_num++, a);
		a->bfs = v;
		bfs_insert_vertex(g, v);

		node = node->next;
	}

	// Insert the lines into the BFS
	node = matching->lines->next;
	while (node != matching->lines) {
		line = (LINE *)node->key;
		a = line->a;
		b = line->b;

		if (line->a->t == LONG_MAX || line->b->t == LONG_MAX) {
			node = node->next;
			continue;
		}

		if ((line->a->t >= 0 && !line->a->bfs) || (line->b->t >= 0 && !line->b->bfs) ) {	
			node = node->next;
			continue;
		}

		if (a->t <= matching->t_bfs) {
			break;
		}

		// If we are next to where we have cut with t_delete, then don't
		// follow lines as they may lead into unknown memory 
		if (matching->t >= matching->t_delete && a->t < matching->t - matching->t_delete) {
			node = node->next;
			continue;
		}

		while (b->merge != NULL) {
			b = b->merge;
		}

		// Boundaries do not appear in the vertices list, and therefore have
		// not yet been allocated to the BFS until they are encountered the
		// first time from a line. 
		if (b->bfs == NULL) {
			// Find a nicer way to do this. When a ball has a line into a
			// ball that has merged into the future, the destination might not
			// exist.
			if (b->ball->type != PRIMAL_BOUNDARY && b->ball->type != DUAL_BOUNDARY) {

				node = node->next;
				continue;
			}
			assert(b->ball->type == PRIMAL_BOUNDARY || b->ball->type == DUAL_BOUNDARY);

			v = bfs_create_vertex(b->v->v_num, b);
			b->bfs = v;
			bfs_insert_vertex(g, v);
		}

		e = bfs_create_edge(a->bfs, b->bfs, line->wt);
		bfs_insert_edge(g, e);

		node = node->next;
	}

	if (undo == TRUE) {
		m_create_undo(matching, UNDO_UPDATE_T_BFS, NULL, matching->t_bfs, NULL, NULL, NULL);
	}

	//printf("Update t_bfs: %d --> %d\n", matching->t_bfs, matching->t);
	matching->t_bfs = latest_t;
}

/**
 * \brief Frees the BFS containing the dots and lines
 * 
 * \param[in] g The \ref BFS_GRAPH to be freed
 */
void m_free_dots_and_lines_bfs(BFS_GRAPH *g) {
	DOT *a;
	CDLL_NODE *node;
	BFS_VERTEX *v;
	
	node = g->vertex_ll->next;
	while (node != g->vertex_ll) {
		v = (BFS_VERTEX *)node->key;
		a = (DOT *)v->key;
		a->bfs = NULL;
		node = node->next;
	}
	bfs_free_graph(g);
}

/**
 * \brief Creates a new structure to track the state of Blossom V
 * 
 * \param[in] num_vertices The number of vertices to be matched
 *
 * \return The newly created \ref BlossomV structure
 */
BLOSSOMV *m_create_blossomv(int num_vertices) {
	BLOSSOMV *bv;

	bv = (BLOSSOMV *)my_malloc(sizeof(BLOSSOMV));
	bv->alloc_edges = 8;

	bv->num_edges = 0;
	// We need twice as many vertices because each vertex will have an
	// accompanying boundary vertex
	bv->num_vertices = 2 * num_vertices;
	bv->edges = (int *)my_malloc(sizeof(int) * 2 * bv->alloc_edges);
	bv->weights = (int *)my_malloc(sizeof(int) * bv->alloc_edges);	
	bv->vertices = (VERTEX **)my_malloc(sizeof(VERTEX *) * bv->num_vertices);

	return bv;
}

/**
 * \brief Frees a Blossom V structure
 * 
 * \param[in] bv The \ref BlossomV structure to be freed
 */
void m_free_blossomv(BLOSSOMV *bv) {
	free(bv->edges);
	free(bv->weights);
	free(bv->vertices);
	free(bv);
}

/**
 * \brief Adds a vertex to the Blossom V state
 * 
 * \param[out] bv The \ref BlossomV structure to insert the \ref vertex into
 * \param[in] v The \ref vertex to insert
 * \param[in] pos The position to insert the \ref vertex into 
 *
 * \return The BlossomV ID, the position in the BlossomV state
 */
int m_add_vertex_to_blossomv(BLOSSOMV *bv, VERTEX *v, int pos) {
	if (pos < 0) {
		pos = (bv->num_vertices / 2) - pos - 1;
	}
	//m_print_void_vertex(v);
	//printf("%d into %d\n", pos, bv->num_vertices);
	bv->vertices[pos] = v;
	return pos;
}

/**
 * \brief Adds an edge to the Blossom V state
 * 
 * \param[out] The \ref BlossomV structure to insert the \ref edge into
 * \param[in] a The BlossomV ID of the source \ref vertex 
 * \param[in] b The BlossomV ID of the destination \ref vertex 
 * \param[in] wt The weight of the edge 
 */
void m_add_edge_to_blossomv(BLOSSOMV *bv, int a, int b, int wt) {
	bv->edges[2*bv->num_edges] = a;
	bv->edges[2*bv->num_edges+1] = b;
	bv->weights[bv->num_edges] = wt;
	bv->num_edges++;

	// Resize if nessecary
	if (bv->num_edges >= bv->alloc_edges) {
		bv->alloc_edges = bv->alloc_edges << 1;
		bv->edges = (int *)my_realloc(bv->edges, sizeof(int) * 2 * bv->alloc_edges);
		bv->weights = (int *)my_realloc(bv->weights, sizeof(int) * bv->alloc_edges);
	}
}

/**
 * \brief Creates an augmented edge between two vertices
 * 
 * \param[out] matching The matching to insert the augmented edge into
 * \param[in] v1 The source vertex of the augmented edge 
 * \param[in] v2 The destination vertex of the augmented edge 
 */
void m_create_augmented_edge(MATCHING *m, VERTEX *v1, VERTEX *v2) {
	CDLL_NODE *n;
	AUG_EDGE *ae;
	
	ae = (AUG_EDGE *)my_malloc(sizeof(AUG_EDGE));
	ae->va = v1;
	ae->vb = v2;	

	n = cdll_create_node(ae);
	cdll_insert_node_head(m->augmented_edges, n);
	ae->n = n;
}

/**
 * \brief Creates a temporal boundary
 * 
 * \return The \ref dot of the new temporal boundary
 */
DOT *m_create_temporal_boundary() {
	VERTEX *v;
	DOT *dot;

	//printf("Create temporal boundary\n");
	dot = m_create_dot(NULL, NULL);
	v = m_create_vertex(-INT_MAX, -INT_MAX, -INT_MAX, dot);
	dot->v = v;
	
	return dot;
}

/**
 * \brief Creates and inserts a vertex and edge into the \ref BFS
 *
 * Creates a new \ref BFS_VERTEX for the destination \ref dot of an
 * edge. Inserting the dot into the \ref BFS and then creating an edge from the
 * source BFS_VERTEX to the newly created destination \ref BFS_VERTEX
 * 
 * \param[out] g The \ref BFS_GRAPH to insert the vertex and edge into
 * \param[in] a The source \ref dot of the edge, already has a \ref BFS_VERTEX 
 * \param[in] b The destination \ref dot of the edge, doesn't have a \ref BFS_VERTEX 
 * \param[in] wt The weight of the edge 
 */
void m_create_and_insert_vertex_and_edge_into_bfs(BFS_GRAPH *g, DOT *a, DOT *b, int wt) {
	BFS_VERTEX *v1, *v2;
	BFS_EDGE *e;

	v1 = a->bfs;

	v2 = bfs_create_vertex(g->num_vertices++, b);
	b->bfs = v2;
	bfs_insert_vertex(g, v2);

	e = bfs_create_edge(v1, v2, wt);
	bfs_insert_edge(g, e);
}

/**
 * \brief Solves the \ref matching of a given \ref BlossomV state
 * 
 * \param[in,out] matching The matching to be solved
 * \param[in] bv The BlossomV state containing the vertices and edges to be matched 
 */
void m_solve(MATCHING *matching, BLOSSOMV *bv) {
	int i;
	struct PerfectMatching::Options options;
	
	// Initialise the BlossomV Perfect Matching
	//printf("num_vertices: %d, num_edges: %d\n", bv->num_vertices, bv->num_edges);
	if (bv->num_vertices > 0 && bv->num_edges < bv->num_vertices / 2) {
		m_print_graph(matching);
		m_print_lattice(matching);
		assert(bv->num_vertices == 0 || bv->num_edges >= bv->num_vertices / 2);
	}
	
	matching->pm = new PerfectMatching(bv->num_vertices, bv->num_edges);
	options.verbose = false;
	options.fractional_jumpstart = true;
	matching->pm->options = options;

	// You can technically feed these in backwards, but this will break the
	// code. The order you read them in, then means Blossom V solves them in
	// the different order which makes reading them back out and then their
	// positions get all messed up. I started doing this then just gave up
	// because I'm not even sure if Blossom V's output is dependent on edge
	// order, or if it is on the vertex numbering. Hence, it wasn't worth the
	// continued effort to implement reverse order edges before I know if it'll
	// actually do anything.
	//printf("%d %d\n", bv->num_vertices, bv->num_edges);
	for (i = 0; i < bv->num_edges; i++) {
		matching->pm->AddEdge(bv->edges[2*i], bv->edges[2*i+1], bv->weights[i]);
	}

	matching->pm->Solve();

	//double cost = ComputePerfectMatchingCost(bv->num_vertices, bv->num_edges, bv->edges, bv->weights, matching->pm);
	//printf("%d %d cost = %.1f\n", bv->num_vertices, bv->num_edges, cost);

	//m_print_graph(matching);
	//m_print_lattice(matching);
}

/**
 * \brief Creates new augmented edges based on newly matched edges that have
 * been modified since the last matching
 * 
 * \param[in] matching The matching to obtain augmented edges from
 * \param[in] bv The \ref BlossomV state 
 * \param[in] undo Whether or not the created augmented edges should be
 * undoable 
 */
void m_create_augmented_edges(MATCHING *matching, BLOSSOMV *bv, int undo) {
	int i;
	VERTEX *v1, *v2;
	EDGE *e, *eb;
	
	// This could perhaps be implemented more efficiently. However, given the
	// large amount of time spent on the BFS, it is currently unknown if this
	// inefficiency would lead to significant enough time savings to justify
	// its implementation. Would require a way to find the edge between two
	// BlossomV NodeIds more efficiently than looping all the edges, which is
	// no more efficient than just looping them all here. If you were to create
	// a hash table, this still requires a loop of the edges to create, and
	// it's not readily apparent how you could keep these hash tables from
	// round to round given that NodeIds change every round.

	for (i = 0; i < bv->num_edges; i++) {
		if (matching->pm->GetSolution(i) == 1) {
			if (bv->weights[i] > 0) {
				v1 = bv->vertices[bv->edges[2*i]];
				v2 = bv->vertices[bv->edges[2*i+1]];

				if (v1->matched_edge && v1->matched_edge->dest_vertex == v2) {
					continue;
				}

				if ((v1->t > matching->t - matching->t_delay && v1->v_num >= 0) || (v2->t > matching->t - matching->t_delay && v2->v_num >= 0)) {
					continue;
				}

				// Delete the old matched edge
				if (v1->matched_edge) {
					assert(v1->matched_edge->source_vertex == v1);

					if (undo == TRUE) {
						// Save where both vertices were previously connected to
						m_create_undo(matching, UNDO_DELETE_EDGE, v1, 0, NULL, v1->matched_edge->dest_vertex, NULL);
					}

					// Need to undo the previous edge
					m_create_augmented_edge(matching, v1, v1->matched_edge->dest_vertex);

					// Clean up any matched edges pointing from old end points 
					v1->matched_edge->dest_vertex->matched_edge = NULL;

					// Delete the existing matched edges
					m_delete_only_edge(v1->matched_edge->backward_edge);
					m_delete_only_edge(v1->matched_edge);

					// Set the matched edges of both vertices to NULL
					v1->matched_edge = NULL;
				}

				if (v2->matched_edge) {
					if (undo == TRUE) {
						m_create_undo(matching, UNDO_DELETE_EDGE, v2, 0, NULL, v2->matched_edge->dest_vertex, NULL);
					}
				
					m_create_augmented_edge(matching, v2, v2->matched_edge->dest_vertex);

					v2->matched_edge->dest_vertex->matched_edge = NULL;

					m_delete_only_edge(v2->matched_edge->backward_edge);
					m_delete_only_edge(v2->matched_edge);
					
					v2->matched_edge = NULL;
				}

				// Create the new matched edge
				e = m_create_only_edge(v1, v2);
				eb = m_create_only_edge(v2, v1);
				e->backward_edge = eb;
				eb->backward_edge = e;

				if (undo == TRUE) {
					m_create_undo(matching, UNDO_MATCHED_EDGE, v1, 0, v1->matched_edge, NULL, NULL);
					if (v2->v_num >= 0) {
						m_create_undo(matching, UNDO_MATCHED_EDGE, v2, 0, v2->matched_edge, NULL, NULL);
					}
				}

				v1->matched_edge = e;
				if (v2->v_num >= 0) {
					v2->matched_edge = eb;
				}

				m_create_augmented_edge(matching, v1, v2);
			}
		}
	}
}

/**
 * \brief Performs all nessecary steps to perform minimum-weight perfect
 * matching
 *
 * Updates the BFS structure, finds and creates the vertices and edges that
 * require matching, builds the BlossomV state, matches using the BlossomV
 * library, then creates augmented edges from the result.
 * 
 * \param[in,out] matching The \ref matching to match
 * \param[in] undo Whether the matching should be undoable 
 */
void m_mwpm(MATCHING *matching, int undo) {
	int i, v_num, v_num1, v_num2, bdy_num;
	int found_bdy, offset;
	DOT *a, *b;
	VERTEX *v1;
	CDLL_NODE *node;

	BLOSSOMV *bv;
	BFS_VERTEX *bfsv1;

	//m_print_graph(matching);

	//printf("\nm_mwpm(): %d\n", matching->t);
	//printf("num_dots: %d\n", matching->num_dots);

	// We want to build a BFS graph of the dots and lines, such that we can
	// traverse the lattice. We will traverse the graph from each dot with a
	// vertex, and look for each other dot with a vertex at it, and construct
	// an edge between them.

	// If there are no vertices in the graph, then do nothing
	node = matching->graph->prev;
	if (node == matching->graph) {
		return;
	}

	if (matching->graph->next == matching->last_graph_head && matching->graph->prev == matching->last_graph_tail) {
		//printf("skipping!\n");
		//return;
	}

	if (undo) {
		m_create_undo(matching, UNDO_UPDATE_LAST_GRAPH_HEAD, NULL, 0, NULL, NULL, matching->last_graph_head);
		m_create_undo(matching, UNDO_UPDATE_LAST_GRAPH_TAIL, NULL, 0, NULL, NULL, matching->last_graph_tail);
	}

	matching->last_graph_head = matching->graph->next;
	matching->last_graph_tail = matching->graph->prev;
	
	// Initialise the Blossom V structure
	bv = m_create_blossomv(matching->num_vertices);

	//printf("Num Vertices: %d, %d\n", matching->num_vertices, bv->num_vertices);

	// Loop over the vertices in the matching and perform a BFS from each
	// connecting each to one another.
	offset = ((VERTEX *)node->key)->v_num; 

	while (node != matching->graph) {
		v1 = (VERTEX *)node->key;
		a = v1->dot;
		assert(a->bfs != NULL);

		// Vertices need to start at 0 for Blossom V. Since we have likely
		// trimmed with t_delete, we need to offset the v_nums to start at 0.
		v_num = a->v->v_num - offset;

		// Boundaries need to start at -1, since -0 is 0 and hence wouldn't be
		// identified as a boundary. 
		bdy_num = -(v_num + 1);

		//printf("v_num: %d, bdy_num: %d\n", v_num, bdy_num);

		v_num1 = m_add_vertex_to_blossomv(bv, a->v, v_num);

		if (a->bfs != NULL) {
			// This init adds about 2s to 100 changes without augmentation 
			bfs_init_search(matching->g, a->bfs);

			// This loop adds about another 6s.
			found_bdy = FALSE;
			bfsv1 = bfs_get_next_v(matching->g);
			while (bfsv1 != NULL) {
				b = (DOT *)bfsv1->key;
				//printf("Inner Dot: ");
				//m_print_dot(b);

				// Keep looking for edges until we have found edges that are
				// longer than the minimum path between the boundaries. After
				// this point, any edge of higher weight will never be chosen
				// as it would be 'cheaper' to connect each vertex to their
				// respective boundaries.
				if (bfsv1->d > matching->D) {
					break;
				}

				// We don't want to link a vertex to itself, or a dot without a
				// vertex
				if (b == a || b->v == NULL) {
					bfsv1 = bfs_get_next_v(matching->g);
					continue;
				}

				// We have a boundary, so we need to "create" a vertex in
				// the Blossom V graph so we can create an edge.
				if (m_is_boundary(b->v)) {
					// Already found a closer boundary for this vertex
					if (found_bdy == TRUE) {
						bfsv1 = bfs_get_next_v(matching->g);
						continue;
					}

					// We found a boundary, congrats. Now on to Carmen Sandiego.
					found_bdy = TRUE;
					v_num2 = m_add_vertex_to_blossomv(bv, b->v, bdy_num); 

					// Connect this boundary to all other inserted boundaries
					for (i = matching->num_vertices; i <= v_num2 - 1; i++) {
						// We do not want to connect normal boundaries to the
						// temporal boundaries, as this leads to matching
						// vertices to things other than the temporal boundary.
						// The temporal boundary is less of a boundary and more
						// of a fixed matching to avoid issues with BlossomV
						// rematching after t_delete.
						if (bv->vertices[i]->t == -INT_MAX) {
							continue;
						}
						m_add_edge_to_blossomv(bv, v_num2, i, 0);
					}	
				} 
				else if (b->v->t == -INT_MAX) {

					// Only create the boundary when the vertex is directly
					// connected to it, not when coming through another vertex.
					if (bfsv1->d != 0) {
						bfsv1 = bfs_get_next_v(matching->g);
						continue;
					}

					found_bdy = TRUE;
					v_num2 = m_add_vertex_to_blossomv(bv, b->v, bdy_num);
				}
				else {
					v_num2 = b->v->v_num - offset;
				}
			
				m_add_edge_to_blossomv(bv, v_num1, v_num2, bfsv1->d);

				bfsv1 = bfs_get_next_v(matching->g);
			}
		} else {
			printf("Error: %d: %d, %d, %ld\n", a->v->v_num, a->i, a->j, a->t);
			assert(a->bfs != NULL);
		}

		node = node->prev;
	}

	// Solve using Blossom V
	m_solve(matching, bv);

	// Create the augmented edges from the solution
	m_create_augmented_edges(matching, bv, undo);

	//if (matching->augmented_edges->next != matching->augmented_edges) {
		//m_print_augmented_edges(matching);
		//m_print_graph(matching);
		//m_print_lattice(matching);
	//}

	delete matching->pm;
	m_free_blossomv(bv);
}

/**
 * \brief Trims all data structures to only contain the timesteps within the
 * block of time defined by t_delete.
 * 
 * \param[in,out] m The \ref matching to trim.
 */
int m_time_delete(MATCHING *m) {
	long int t, t2, v_del;
	CDLL_NODE *n, *n2;
	DOT *dot;
	LINE *line;
	VERTEX *v, *v2;
	DOT *a, *b;
	BFS_VERTEX *bfsv1, *bfsv2;
	BFS_EDGE *bfse;
	int hash;
	long int pos[3];

	t = m->t - m->t_delete;
	t2 = t - m->t_delete;

	// Delete all temporal boundaries that will have their paired vertex
	// deleted by time delete. 
	//
	// Temporal boundaries do not exist in the graph, only in the BFS
	// structure. They could be found by looping each vertex and looking over
	// the BFS_EDGEs, however this would be less efficient.

	n = m->g->vertex_ll->next;
	while (n != m->g->vertex_ll) {
		bfsv1 = (BFS_VERTEX *)n->key;
		n = n->next;
	
		// We only care about the vertices that are temporal boundaries
		a = (DOT *)bfsv1->key;
		if (a->v && a->v->t == -INT_MAX) {
			// A temporal vertex always has only one edge
			bfse = (BFS_EDGE *)bfsv1->edge_ll->next->key;
			bfsv2 = (bfse->a == bfsv1) ? bfse->b : bfse->a;
			b = (DOT *)bfsv2->key;
			
			// If the vertex this temporal boundary is connected to will be
			// deleted, then clean up the temporal boundary. 
			if (b->v->t < t) {
				bfs_remove_vertex(m->g, bfsv1);
				m_delete_vertex(m, a->v);
				m_free_dot(a);
			}
		}
	}	

	// Delete all vertices that are less than t in the past.
	v_del = 0;
	n = m->graph->prev;
	while (n != m->graph) {
		v = (VERTEX *)n->key;
		if (v->t >= t) {
			break;
		}

		// If we are deleting a vertex with a matched edge, we need to check if
		// it is connected to a vertex that wont be deleted during this round
		// of time_delete. 
		//
		// If the edge is to a retained vertex, then we need to create a
		// pseudo-dot to connect the vertex to instead, and then insert this
		// into the BFS tree so that it can be found. 
		//
		// We insert an edge to this new boundary of 0 weight because we do not
		// expect the matching to change near the t_delete boundary. If we
		// expect this edge to be required, then we should not be using
		// time_delete in the first place because this will give us erroneous
		// results.
		if (v->matched_edge) {
			v2 = v->matched_edge->dest_vertex;
			if (v2->t >= t && v2->v_num >= 0) {
				a = v2->dot;
				b = m_create_temporal_boundary();
				m_create_and_insert_vertex_and_edge_into_bfs(m->g, a, b, 0);
			}
			
			// Delete the matched edge
			m_delete_only_edge(v->matched_edge->backward_edge);
			m_delete_only_edge(v->matched_edge);
			v->matched_edge = NULL;
			v2->matched_edge = NULL;
		}	

		n = n->prev;
		v_del++;
		m_delete_vertex(m, v);
	}

	// Delete all lines that are less than t2 in the past.
	n = m->lines->prev;
	while (n != m->lines) {
		n2 = n->prev;
		line = (LINE *)n->key;
		a = line->a;
		b = line->b;

		if (line->a->t < t2 || line->b->t < t2) {
			// We don't want to delete boundary lines until the source vertex
			// is past t_delete. Boundaries actually have low t value rather
			// than LONG MAX. This rarely has an impact, since this will only
			// delete the boundary lines close to the time boundary within a
			// small range.
			if (line->a->t >= t2 && line->b->v && line->b->v->v_num < 0) {

				n = n2;	
				continue;
			}

			cdll_delete_node(n, free);
		}
		n = n2;	
	}

	// Delete all dots that are less than t in the past.
	n = m->dots->prev;
	while (n != m->dots) {
		dot = (DOT *)n->key;
		if (dot == NULL || dot->t >= t2) {
			break;
		}
		
		if (dot->bfs != NULL) {
			bfs_remove_vertex(m->g, dot->bfs);
		}
			
		pos[0] = dot->i;
		pos[1] = dot->j;
		pos[2] = dot->little_t;
		hash = abs(fasthash32(pos, 3, 0));
		ht_delete_key(m->dot_ht, hash, dot, NULL);

		cdll_delete_node(n, m_free_void_dot);

		n = m->dots->prev;
	}

	return v_del;
}

/**
 * \brief Gets the first augmented edge
 * 
 * \param[in] m The \ref matching to get the augmented edge from
 */
AUG_EDGE *m_get_aug_edge(MATCHING *m) {
	CDLL_NODE *n;

	//print_tree(m);

	n = m->augmented_edges;
	if (n->prev == n) return NULL;

	return (AUG_EDGE *)n->prev->key;
}

/**
 * \brief Deletes an augmented edge
 * 
 * \param[in] ae The augmented edge to delete
 */
void m_delete_aug_edge(AUG_EDGE *ae) {
	cdll_delete_node(ae->n, free);
}

/**
 * \brief Creates a new dot from a \ref ball and a \ref vertex.
 *
 * If there is no BALL provided, then i = j = t = -1 and dot->ball will be
 * NULL. If no vertex is provided then dot->v will be NULL.
 * 
 * \param[in] ball The \ref ball to use for i,j,t coordinates
 * \param[in] vertex The \ref vertex associated with the dot 
 */
DOT *m_create_dot(BALL *ball, VERTEX *v) {
	DOT *dot;
	
	dot = (DOT *)my_malloc(sizeof(DOT));

	if (ball != NULL) {
		dot->i = ball->i;
		dot->j = ball->j;
		dot->t = ball->big_t;
		dot->little_t = ball->t;
	} else {
		dot->i = -1;
		dot->j = -1;
		dot->t = -1;
		dot->little_t = -1;
	}

	dot->merge = NULL;

	dot->lines = NULL;
	dot->bfs = NULL;

	dot->ball = ball;
	if (ball != NULL) {
		ball->dot = dot;
	}

	dot->v = v;
	if (v != NULL) {
		v->dot = dot;
	}

	return dot;
}

/**
 * \brief Inserts a \ref dot into a \ref matching
 * 
 * \param[out] m The \ref The \ref matching to insert the \ref dot into
 * \param[in] dot The \ref dot to insert  
 */
void m_insert_dot(MATCHING *m, DOT *dot) {
	CDLL_NODE *n;

	n = cdll_create_node(dot);
	cdll_insert_node_head(m->dots, n);
	m->num_dots++;

	if (m->undo_flag == TRUE) {
		m_create_undo(m, UNDO_CREATE_AND_INSERT_DOT, NULL, 0, NULL, NULL, n);
	}
}

/**
 * \brief Creates and inserts a new \ref dot 
 * 
 * \param[out] m The \ref matching to insert the \ref dot into
 * \param[in] ball The \ref ball to base the new \ref dot on 
 */
DOT *m_create_and_insert_dot(MATCHING *m, BALL *ball) {
	DOT *dot;
	CDLL_NODE *n;
	
	dot = m_create_dot(ball, NULL);

	n = cdll_create_node(dot);
	cdll_insert_node_head(m->dots, n);
	m->num_dots++;

	ht_insert_key(m->dot_ht, ball->hash_ijt, dot);

	if (m->undo_flag == TRUE) {
		m_create_undo(m, UNDO_CREATE_AND_INSERT_DOT, NULL, 0, NULL, NULL, n);	
	}

	return dot;
}

/**
 * \brief Creates and inserts a new \ref dot and \ref vertex
 * 
 * \param[out] m The \ref matching to insert the \ref dot and \ref vertex into
 * \param[in] ball The \ref ball to base the new \ref dot and \ref vertex on 
 */
DOT *m_create_and_insert_dot_and_vertex(MATCHING *m, BALL *ball) {
	VERTEX *v;
	DOT *dot;
	CDLL_NODE *n;
	
	dot = m_create_dot(ball, NULL);
	
	if (ball->type != PRIMAL_BOUNDARY && ball->type != DUAL_BOUNDARY) {
		n = cdll_create_node(dot);
		cdll_insert_node_head(m->dots, n);
		m->num_dots++;
		if (m->undo_flag == TRUE) {
			m_create_undo(m, UNDO_CREATE_AND_INSERT_DOT, NULL, 0, NULL, NULL, n);
		}
		ht_insert_key(m->dot_ht, ball->hash_ijt, dot);
	}

	v = m_create_vertex(ball->i, ball->j, ball->big_t, dot);
	
	if (ball->type != PRIMAL_BOUNDARY && ball->type != DUAL_BOUNDARY) {
		m_insert_vertex(m, v);	
		if (m->undo_flag == TRUE) {
			m_create_undo(m, UNDO_CREATE_VERTEX, v, 0, NULL, NULL, NULL);
		}
	}
	else {
		v->v_num = ball->i;

		// Prevents time abort as boundary sets
		// create boundary dots with t = INT_MAX
		dot->t = -1;
	}

	return dot;
}

/**
 * \brief Frees a \ref dot
 * 
 * \param[in] dot The \ref dot to be freed
 *
 * \sa m_free_void_dot
 */
void m_free_dot(DOT *dot) {
	// If the dot has a vertex and it is a boundary, then we need to free the
	// vertex as well as the dot to avoid it becoming detached.
	if (dot->v != NULL && dot->v->v_num < 0) {
		m_fast_delete_vertex(dot->v);
	}

	if (dot->ball != NULL) {
		dot->ball->dot = NULL;
	}
	
	if (dot->v != NULL) {
		dot->v->dot = NULL;
	}

	ll_free(dot->lines, NULL);
	free(dot);
}

/**
 * \brief Frees a void \ref dot
 * 
 * \param[in] key The \ref dot to be freed
 *
 * \sa m_free_dot
 */
void m_free_void_dot(void *key) {
	DOT *dot;
	dot = (DOT *)key;
	m_free_dot(dot);
}

/**
 * \brief Creates and inserts a new \ref line using a p value from a stick
 * 
 * \param[out] m The \ref matching to insert into
 * \param[in] da The source \ref dot 
 * \param[in] db The destination \ref dot
 * \param[in] p_stick The value of p from the associated stick
 *
 * \sa create_and_insert_line_wt
 */
LINE *m_create_and_insert_line(MATCHING *m, DOT *da, DOT *db, float p_stick) {
	LINE *line;
	CDLL_NODE *n;

	line = (LINE *)my_malloc(sizeof(LINE));
	line->a = da;
	line->b = db;
	line->wt = -log(p_stick)*PRECISION;
	if (line->wt < PRECISION) {
		line->wt = PRECISION;
	}
	line->wt -= line->wt%2;

	assert(da->v == NULL || da->v->v_num >= 0);
	da->lines = ll_insert(da->lines, line);		
	
	// Check that destination dot is not a boundary dot. 
	if (db->v == NULL || db->v->v_num >= 0) {
		db->lines = ll_insert(db->lines, line);
	}

	n = cdll_create_node(line);
	cdll_insert_node_head(m->lines, n);
	if (m->undo_flag == TRUE) {
		m_create_undo(m, UNDO_CREATE_AND_INSERT_LINE, NULL, 0, NULL, NULL, n);
	}
	m->num_lines++;

	return line;
}

LINE *m_create_line_wt(MATCHING *m, DOT *da, DOT *db, int wt) {
	LINE *line;
	CDLL_NODE *n;

	line = (LINE *)my_malloc(sizeof(LINE));
	line->a = da;
	line->b = db;
	line->wt = wt;

	assert(da->v == NULL || da->v->v_num >= 0);
	da->lines = ll_insert(da->lines, line);		

	if (db->v == NULL || db->v->v_num >= 0) {
		db->lines = ll_insert(db->lines, line);
	}

	if (m->undo_flag == TRUE) {
		n = cdll_create_node(line);
		m_create_undo(m, UNDO_CREATE_LINE, NULL, 0, NULL, NULL, n);
	}

	return line;
}

void m_insert_line(MATCHING *m, LINE *line) {
	CDLL_NODE *n, *x;
	LINE *l;
	assert(line->a->t != LONG_MAX);
	assert(line->b->t != LONG_MAX);

	x = m->lines->next;
	n = cdll_create_node(line);
	while (x != m->lines) {
		l = (LINE *)x->key;
		if (l->a->t <= line->a->t) {
			break;
		}
		x = x->next;
	}

	n->next = x;
	n->prev = x->prev;
	x->prev->next = n;
	x->prev = n;

	if (m->undo_flag == TRUE) {
		m_create_undo(m, UNDO_INSERT_LINE, NULL, 0, NULL, NULL, n);
	}
	m->num_lines++;
}

/**
 * \brief Creates and inserts a line using a weight
 * 
 * \param[out] m The \ref matching to insert into
 * \param[in] da The source \ref dot 
 * \param[in] db The destination \ref dot
 * \param[in] wt The weight of the \ref line 
 *
 * \return The newly created line
 *
 * \sa m_create_and_insert_line
 */
LINE *m_create_and_insert_line_wt(MATCHING *m, DOT *da, DOT *db, int wt) {
	LINE *line;

	line = m_create_line_wt(m, da, db, wt);
	m_insert_line(m, line);

	return line;
}

/**
 * \brief Checks to see if a \ref vertex is a boundary
 * 
 * \param[in] v The \ref vertex to be checked
 *
 * \return A boolean of whether the \ref vertex is a boundary
 */
bool m_is_boundary(VERTEX *v) {
	return (v->v_num < 0);
}


/**
 * \brief Creates a new \ref vertex
 * 
 * \param[in] i The i coordinate
 * \param[in] j The j coordinate
 * \param[in] t The t coordinate
 * \param[in] dot The \ref dot to associate with the \ref vertex
 *
 * \return The newly created \ref vertex
 */
VERTEX *m_create_vertex(int i, int j, long int t, DOT *dot) {
	VERTEX *v;

	//printf("Create vertex\n");
	v = (VERTEX *)my_malloc(sizeof(VERTEX));

	v->v_num = 0;
	v->i = i;
	v->j = j;
	v->t = t;
	v->dot = dot;
	
	if (dot != NULL) {
		dot->v = v;
	}

	v->matched_edge = NULL;
	v->graph_node = NULL;
	
	return v;
}

/**
 * \brief Sets the uplink from the \ref vertex to its position in the graph
 * 
 * \param[in] key The \ref vertex to set the uplink in
 * \param[in] graph_node The ::CDLL_NODE to link to 
 */
void m_set_vertex_graph_node(void *key, CDLL_NODE *graph_node) {
	((VERTEX *)key)->graph_node = graph_node;
}

/**
 * \brief Inserts a \ref vertex into the \ref matching 
 * 
 * \param[in] m The \ref matching to insert the \ref vertex into
 * \param[in] v The \ref vertex to insert 
 */
void m_insert_vertex(MATCHING *m, VERTEX *v) {
	VERTEX *v_head;

	// Get the next v_num to allocate to the vertex
	v_head = (VERTEX *)m->graph->next->key;
	if (v_head != NULL) {
		assert(v_head->v_num < INT_MAX);
		v->v_num = v_head->v_num + 1;
	}

	cdll_insert_head(m->graph, v, m_set_vertex_graph_node); 
	m->num_vertices++;
}

/**
 * \brief Deletes a \ref vertex from a \ref matching
 * 
 * \param[in] m The \ref matching the \ref vertex exists in
 * \param[in] v The \ref vertex to delete 
 */
void m_delete_vertex(MATCHING *m, VERTEX *v) {
	// Temporal boundaries are never inserted into the graph, so we only want
	// to remove the vertex from the graph when it is an actual vertex
	if (v->v_num >= 0 && v->graph_node) {
		cdll_delete_node(v->graph_node, NULL);	
		m->num_vertices--;
	}
	
	m_fast_delete_vertex(v);
}

/**
 * \brief Deletes a \ref vertex from a \ref matching quickly. Doesn't clean up
 * the graph list in matching or the number of vertices.
 * 
 * \param[in] v The \ref vertex to delete 
 *
 * \sa fast_delete_void_vertex
 */
void m_fast_delete_vertex(VERTEX *v) {
	// If the vertex still has a matched edge, we need to clean it up
	if (v->matched_edge != NULL) {
		v->matched_edge->dest_vertex->matched_edge = NULL;
		m_delete_only_edge(v->matched_edge->backward_edge);	
		m_delete_only_edge(v->matched_edge);
	}

	// Ensure there are no dots pointing to a deleted vertex
	if (v->dot != NULL) {
		v->dot->v = NULL;
	}

	free(v);
}

/**
 * \brief Deletes a \ref vertex from a \ref matching quickly. Accepts a void
 * pointer to a vertex.
 * 
 * \param[in] key The \ref vertex to be freed
 *
 * \sa fast_delete_vertex
 */
void m_fast_delete_void_vertex(void *key) {
	m_fast_delete_vertex((VERTEX *)key);
}

/**
 * \brief Creates an edge between two \ref vertex
 * 
 * \param[in] v1 The first \ref vertex to connect
 * \param[in] v2 The second \ref vertex to connect
 */
EDGE *m_create_only_edge(VERTEX *v1, VERTEX *v2) {
	EDGE *e = (EDGE *)my_malloc(sizeof(EDGE));
	e->source_vertex = v1;
	e->dest_vertex = v2;
	return e;
}

/**
 * \brief Deletes an edge from between two \ref vertex
 * 
 * \param[in] e The edge to be freed
 */
void m_delete_only_edge(EDGE *e) {
	free(e);
}

/**
 * \brief Creates a new \ref undo structure
 * 
 * \param[out] matching The \ref matching to insert the \ref undo into
 * \param[in] type The type of \ref undo
 * \param[in] v A \ref vertex (NULL if not required)
 * \param[in] inc An increment (0 if not required)
 * \param[in] e An \ref edge  (NULL if not required)
 * \param[in] v2 A \ref vertex (NULL if not required)
 * \param[in] n A ::CDLL_NODE (NULL if not required) 
 */
void m_create_undo(MATCHING *matching, int type, VERTEX *v, long int inc, EDGE *e, VERTEX *v2, CDLL_NODE *n) {
	UNDO *ut;
	ut = (UNDO *)my_malloc(sizeof(UNDO));

	ut->type = type;
	ut->v = v;
	ut->inc = inc;
	ut->e = e;
	ut->v2 = v2;
	ut->n = n;
	ut->d1 = NULL;
	ut->d2 = NULL;
	ut->line = NULL;
	ut->t = LONG_MAX;

	ut->next = matching->u;
	matching->u = ut;
}

void m_create_line_undo(MATCHING *matching, int type, LINE *line, DOT *d1, DOT *d2, int wt, long int t) {
	UNDO *ut;
	ut = (UNDO *)my_malloc(sizeof(UNDO));

	ut->type = type;
	ut->v = NULL;
	ut->inc = wt;
	ut->e = NULL;
	ut->v2 = NULL;
	ut->n = NULL;
	ut->d1 = d1;
	ut->d2 = d2;
	ut->line = line;
	ut->t = t;

	ut->next = matching->u;
	matching->u = ut;
}


/**
 * \brief Executes and deletes the \ref undo%s
 * 
 * \param[in] matching The \ref matching containing all the \ref undo%s
 */
void m_execute_and_delete_undos(MATCHING *matching) {
	CDLL_NODE *n;
	LL_NODE *ll_node;
	LINE *line;
	DOT *dot;
	EDGE *e, *eb;
	int hash;
	long int pos[3];
	
	UNDO *ut;
	ut = matching->u;

	while (ut != NULL) {
		if (ut->type == UNDO_CREATE_VERTEX) {
			// ut->v must be the head of the graph
			//printf(":undo create vertex\n");
			m_delete_vertex(matching, ut->v);
		}

		else if (ut->type == UNDO_UPDATE_LAST_GRAPH_HEAD) {
			matching->last_graph_head = ut->n;
		}
		
		else if (ut->type == UNDO_UPDATE_LAST_GRAPH_TAIL) {
			matching->last_graph_tail = ut->n;
		}

		else if (ut->type == UNDO_UPDATE_T_BFS) {
			matching->t_bfs = ut->inc;
		}

		else if (ut->type == UNDO_DELETE_EDGE) {
			e = m_create_only_edge(ut->v, ut->v2);
			eb = m_create_only_edge(ut->v2, ut->v);
			e->backward_edge = eb;
			eb->backward_edge = e;

			ut->v->matched_edge = e;
			if (ut->v2->v_num >= 0) {
				ut->v2->matched_edge = eb;
			}
		}

		else if (ut->type == UNDO_MATCHED_EDGE) {
			assert(ut->e == NULL);
			
			if (ut->v->matched_edge->dest_vertex->v_num < 0) {
				m_delete_only_edge(ut->v->matched_edge->backward_edge);
			}
			
			m_delete_only_edge(ut->v->matched_edge);	
			ut->v->matched_edge = ut->e;
		}

		// Undo create and insert dot
		else if (ut->type == UNDO_CREATE_AND_INSERT_DOT) {
			n = ut->n;

			// Traverse the lines ll in the dot and delete the back reference
			// from the line back to the dot we are deleting.
			//
			// This design choice was made due to the lines list being at most 12
			// lines in size, and hence traversing was decided to be more efficient
			// than setting up a cdll and including pointers in the line to each of
			// the lines lists of its associated dots.


			dot = (DOT *)n->key;
			ll_node = dot->lines;
			while (ll_node != NULL) {
				line = ((LINE *)ll_node->key);
				if (line->a == dot) {
					line->a = NULL;
				} else {
					line->b = NULL;
				}
				ll_node = ll_node->next;
			}

			pos[0] = dot->i;
			pos[1] = dot->j;
			pos[2] = dot->little_t;
			hash = abs(fasthash32(pos, 3, 0));
			ht_delete_key(matching->dot_ht, hash, dot, NULL);

			if (dot->bfs != NULL) {
				bfs_remove_vertex(matching->g, dot->bfs);
				dot->bfs = NULL;
			}

			// Delete the cdll_node pointing to this dot from matching->dots
			cdll_delete_node(n, NULL);

			// Free the dot and delete the lines linked list backbone
			m_free_dot(dot);
			matching->num_dots--;
		}

		// Undo create and insert line
		else if (ut->type == UNDO_CREATE_AND_INSERT_LINE) {
			n = ut->n;
			line = (LINE *)n->key;

			// Iterate the lines linked lists for each dot associated
			// with the line being deleted. Find the line being deleted,
			// and remove it from the linked list. 

			if (line->a != NULL) {
				dot = line->a;
				dot->lines = ll_delete_node(dot->lines, n->key);
			}
			if (line->b != NULL) {
				dot = line->b;
				dot->lines = ll_delete_node(dot->lines, n->key);
			}

			// Now that the back links to the line have been cleaned up,
			// we can safely remove the line and its cdll node.

			cdll_delete_node(n, NULL);
			free(line);
			matching->num_lines--;
		}

		else if (ut->type == UNDO_CREATE_LINE) {
			n = ut->n;
			line = (LINE *)n->key;

			//printf(":undo create line ");
			//m_print_line(line);

			// Iterate the lines linked lists for each dot associated
			// with the line being deleted. Find the line being deleted,
			// and remove it from the linked list. 
			
			if (line->a != NULL) {
				dot = line->a;
				if (dot->merge != NULL) {
					dot = dot->merge;
				}
				dot->lines = ll_delete_node(dot->lines, n->key);
			}

			if (line->b != NULL) {
				dot = line->b;
				if (dot->merge != NULL) {
					dot = dot->merge;
				}
				dot->lines = ll_delete_node(dot->lines, n->key);
			}

			free(line);
			free(n);
		}

		else if (ut->type == UNDO_INSERT_LINE) {
			n = ut->n;
			//printf(":undo insert line\n");
			line = (LINE *)n->key;
			cdll_delete_node(n, NULL);
			matching->num_lines--;
		}


		else if (ut->type == UNDO_CREATE_DOT) {
			n = ut->n;
			dot = (DOT *)n->key;
			//printf(":undo create dot\n");
			ll_node = dot->lines;
			while (ll_node != NULL) {
				line = ((LINE *)ll_node->key);
				if (line->a == dot) {
					line->a = NULL;
				} else {
					line->b = NULL;
				}
				ll_node = ll_node->next;
			}

			pos[0] = dot->i;
			pos[1] = dot->j;
			pos[2] = dot->little_t;
			hash = abs(fasthash32(pos, 3, 0));
			ht_delete_key(matching->dot_ht, hash, dot, NULL);

			// Free the dot and delete the lines linked list backbone
			m_free_dot(dot);
			free(n);
		}

		else if (ut->type == UNDO_INSERT_DOT) {
			n = ut->n;
			dot = (DOT *)n->key;

			//printf(":undo insert dot\n");

			// Delete the cdll_node pointing to this dot from matching->dots
			cdll_delete_node(n, NULL);

			if (dot->bfs != NULL) {
				bfs_remove_vertex(matching->g, dot->bfs);
				dot->bfs = NULL;
			}

			dot->t = ut->t;

			// Free the dot and delete the lines linked list backbone
			matching->num_dots--;
		}


		else if (ut->type == UNDO_LINE_EDIT) {
			//printf("undo line edit: ");
			ut->line->a = ut->d1;
			ut->line->b = ut->d2;
		}

		else if (ut->type == UNDO_LINE_MOVE) {
			//printf("undo move: ");
			ut->d1->lines = ll_insert(ut->d1->lines, ut->line);
			if (ut->d2 != NULL) {
				ut->d2->lines = ll_delete_node(ut->d2->lines, ut->line);
			}
		}

		else if (ut->type == UNDO_LINE_WT) {
			ut->line->wt = ut->inc;
		}

		else if (ut->type == UNDO_DOT_MERGE) {
			ut->d1->merge = ut->d2;
		}

		else if (ut->type == UNDO_DOT_T) {
			ut->d1->t = ut->t;
		}

		ut = ut->next;
	}
	
	m_free_undo(matching);
}

// PRINT FUNCTIONS

/**
 * \brief Prints a \ref matching
 * 
 * \param[in] matching The \ref matching to print
 */
void m_print_matching(MATCHING *matching) {
	if (matching != NULL) {
		if (matching->graph != NULL) {
			m_print_graph(matching);
		}
	}
}

/**
 * \brief Prints the lattice from a \ref matching
 * 
 * \param[in] matching The \ref matching to print the lattice of
 */
void m_print_lattice(MATCHING *matching) {
	CDLL_NODE *n;
	DOT *dot, *dd;
	LL_NODE *lln;
	LINE *line;

	n = matching->dots->next;
	while (n != matching->dots) {
		dot = (DOT *)n->key;
		if (dot->t < matching->t-10) break;
		if (dot->v != NULL) printf("v(%d, %d, %d, %ld) ", dot->v->v_num, dot->v->i, dot->v->j, dot->v->t);
		printf("(%d, %d, %ld) - ", dot->i, dot->j, dot->t);
		lln = dot->lines;
		while (lln != NULL) {
			line = (LINE *)lln->key;
			dd = (line->a == dot) ? line->b : line->a;
			if (dd->v != NULL) printf("v(%d, %d, %d, %ld) ", dd->v->v_num, dd->v->i, dd->v->j, dd->v->t);
			printf("(%d, %d, %ld) [%d] ", dd->i, dd->j, dd->t, line->wt);
			lln = lln->next;
		}
		printf("\n");
		n = n->next;
	}
	printf("\n");
}

/**
 * \brief Prints a \ref line
 * 
 * \param[in] line The \ref line to print
 */
void m_print_line(LINE *line) {	
	printf("line %p a (%d, %d, %ld, %ld, %p)", line, line->a->i, line->a->j, line->a->little_t, line->a->t, line->a);
	if (line->a->v != NULL) printf(" v: %d", line->a->v->v_num);
	printf(" b (%d, %d, %ld, %ld, %p)", line->b->i, line->b->j, line->b->little_t, line->b->t, line->b);
	if (line->b->v != NULL) printf(" v: %d", line->b->v->v_num);
	printf("\n");
}

void m_print_void_line(void *key) {
	m_print_line((LINE *)key);
}

/**
 * \brief Prints a \ref dot
 * 
 * \param[in] dot The \ref dot to print
 *
 * \sa m_print_void_dot
 */
void m_print_dot(DOT *dot) {
	printf("dot (%d, %d, %ld, %ld, %p, %p)", dot->i, dot->j, dot->little_t, dot->t, dot, dot->ball);
	if (dot->v != NULL) {
		printf(" v: %d", dot->v->v_num);
		if (dot->v->v_num < 0) printf(" ---bdy---");
	}
	printf("\n");
}

/**
 * \brief Prints a \ref dot 
 * 
 * \param[in] key The \ref dot to print
 *
 * \sa m_print_dot
 */
void m_print_void_dot(void *key) {
	m_print_dot((DOT *)key);
}

/**
 * \brief Prints the graph from a \ref matching
 * 
 * \param[in] matching The \ref matching to print
 */
void m_print_graph(MATCHING *matching) {
	VERTEX *v;

	if ((matching == NULL) || (matching->graph->next == matching->graph)) {
		return;
	}

	printf("num_v: %d\n", matching->num_vertices);

	v = (VERTEX *)matching->graph->next->key;
	while (v != NULL) {
		m_print_void_vertex(v);
		if (v->matched_edge) {
			m_print_void_edge(v->matched_edge);
		}
		v = (VERTEX *)v->graph_node->next->key;
	}
}

/**
 * \brief Prints a \ref vertex
 * 
 * \param[in] k The \ref vertex to print
 */
void m_print_void_vertex(void *k) {
	VERTEX *v;

	v = (VERTEX *)k;

	if (v == NULL) {
		printf("null\n");
		return;
	}

	printf("v: %d (%d, %d, %ld, %ld)\n", v->v_num, v->i, v->j, v->t, v->dot->little_t);
}

/**
 * \brief Prints the list of augmented edges
 * 
 * \param[in] matching The \ref matching to print the augmented edges from
 */
void m_print_augmented_edges(MATCHING *matching) {
	CDLL_NODE *n;
	AUG_EDGE *ae;

	n = matching->augmented_edges->next;

	while (n != matching->augmented_edges) {
		ae = (AUG_EDGE *)n->key;
		printf("augmented edge: %d (%d,%d) to %d (%d,%d)\n", ae->va->v_num, ae->va->i, ae->va->j, ae->vb->v_num, ae->vb->i, ae->vb->j);
		n = n->next;
	}
}

/**
 * \brief Prints an \ref edge
 * 
 * \param[in] ev The \ref edge to print
 */
void m_print_void_edge(void *ev) {
	EDGE *e = (EDGE *)ev;

	if (ev == NULL) {
		printf("Edge is NULL\n");
		return;
	}

	printf("edge (%d, %d)\n", e->source_vertex->v_num, e->dest_vertex->v_num); 
}
