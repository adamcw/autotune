#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include "qc.h"
#include "../memory/memory.h"
#include "../random/random.h"
#include "../fast_hash/fasthash.h"
#include <sys/time.h>

//#define DBG_MEAS_QUBIT 1
//#define DBG_CREATE_SET 1

/** The number of timesteps in each correction cycle */
#define CYCLE_T 8

#define RECIPE_ADV_OFFSET_HT_SIZE 1000117
#define RECIPE_ADV_BLOCK_HT_SIZE 1000039
#define RECIPE_ADV_DOT_HT_SIZE 1000081
#define MERGE_LINES_HT_SIZE 1000099

double t_start[TIMERS];
double t_total[TIMERS];
int t_counts[TIMERS];
const char *t_name[TIMERS];

/**
 * \brief A crude timing class for determining areas of code using large
 * amounts of time
 * 
 * This timer function is primarily usefull for finding code that becomes
 * slower over time (rather than running in constant time). Profilers can find
 * slow running functions, but narrowing it down to a specific call on a
 * specific data structure that is slowing over time can be difficult.
 *
 * \param[in] timer_num A number less than TIMERS that is the unique ID for a
 * timer
 * \param[in] action The action to start/stop/print the timers. 1 to start, 0
 * to stop, anything else will print the timers.
 * \param[in] name A string representing the name of this timer. Printed with
 * the timers to make it easier to view which timers are which
 */
void timer(int timer_num, int action, const char *name) {
	struct timeval t1970;
	double time;
	int i;
	gettimeofday(&t1970, NULL);
	time = (double)t1970.tv_sec + (double)t1970.tv_usec/1000000;
	if (action == START) {
		if (timer_num == 0) {
			for (i = 0; i < TIMERS; i++) {
				t_total[i] = 0;
				t_counts[i] = 0;
			}
		}
		t_counts[timer_num]++;
		t_name[timer_num] = name;
		t_start[timer_num] = time; 
	}
	else if (action == STOP) {
		t_total[timer_num] += time - t_start[timer_num];
		//printf("Timer %s Increment: %f\n", t_name[timer_num], time - t_start[timer_num]);
	} 
	else {
		for (i = 0; i < TIMERS; i++) {
			if (t_total[i] > 0 && t_counts[i] % 50 == 0) {
				printf("Timer: %50s\t%d\t%10.10f\t(%0.2f%%)\n", t_name[i], t_counts[i], t_total[i], (t_total[i] / t_total[0])*100); 
			}		
		}
	}
}



/**
 * \brief Initialise the quantum computuation, with an original recipe 
 *
 * \param[in] s0 The first seed of the random number generator
 * \param[in] s1 The second seed of the random number generator
 * \param[in] de_ht_size The size of the \ref de hash table
 * \param[in] stick_ht_size The size of the \ref stick hash table
 * \param[in] t_delete How many timesteps in the past to keep before deletion 
 * \param[in] recipe The \ref recipe for the quantum computation 
 *
 * \callgraph
 * \callergraph
 *
 * \return An initialised quantum computer
 */
QC *qc_create_qc(int s0, int s1, int de_ht_size, int stick_ht_size, int t_delete, RECIPE *recipe) {
	QC *qc;

	qc = (QC *)my_malloc(sizeof(QC));
	qc->recipe = recipe;
	qc->recipe_adv = NULL;
	qc_init_qc(qc, s0, s1, de_ht_size, stick_ht_size, t_delete);

	return qc;
}

/**
 * \brief Initialise the quantum computuation, with an advanced recipe
 *
 * \param[in] s0 The first seed of the random number generator
 * \param[in] s1 The second seed of the random number generator
 * \param[in] de_ht_size The size of the \ref de hash table
 * \param[in] stick_ht_size The size of the \ref stick hash table
 * \param[in] t_delete How many timesteps in the past to keep before deletion 
 * \param[in] recipe The \ref recipe_adv for the quantum computation 
 *
 * \callgraph
 * \callergraph
 *
 * \return An initialised quantum computer
 */
QC *qc_create_qc_adv(int s0, int s1, int de_ht_size, int stick_ht_size, int t_delete, RECIPE_ADV *recipe) {
	QC *qc;

	qc = (QC *)my_malloc(sizeof(QC));
	qc->recipe = NULL;
	qc->recipe_adv = recipe;
	qc_init_qc(qc, s0, s1, de_ht_size, stick_ht_size, t_delete);

	return qc;
}

/**
 * \brief Sets the default for a \ref qc
 * 
 * Initialises the quantum computer. Notable defaults include: 
 *  - perfect_gates set to FALSE
 *  - track set to TRUE
 *
 * \param[in] s0 The first seed of the random number generator
 * \param[in] s1 The second seed of the random number generator
 * \param[in] de_ht_size The size of the \ref de hash table
 * \param[in] stick_ht_size The size of the \ref stick hash table
 * \param[in] t_delete How many timesteps in the past to keep before deletion 
 * \param[in] recipe The \ref recipe_adv for the quantum computation 
 */
void qc_init_qc(QC *qc, int s0, int s1, int de_ht_size, int stick_ht_size, int t_delete) {
	randoms(&s0, &s1);
	//Init0(0);

	qc->s0 = s0;
	qc->s1 = s1;
	qc->big_t = 0;
	qc->unfinalized_big_t = 0;
	qc->num_qubits = 0;
	qc->qubit_cdll = cdll_create();
	qc->syn_heap = bh_create(0);
	qc->set_heap = bh_create(0);
	qc->finalized_balls = cdll_create();
	qc->zombie_dots = cdll_create();
	qc->num_ems = 0;
	qc->ems = NULL;
	qc->num_es = 0;
	qc->next_label = 0;
	qc->num_errors = 0;
	qc->error_cdll = cdll_create();
	qc->num_nest_errors = 0;
	qc->nest_error_cdll = cdll_create();
	qc->num_des = 0;
	qc->de_ht = ht_create(de_ht_size);
	qc->de_cdll = cdll_create();
	qc->nest_pr = qc_create_nest(stick_ht_size);
	qc->nest_du = qc_create_nest(stick_ht_size);
	qc->m_pr = m_create_matching(t_delete);
	qc->m_du = m_create_matching(t_delete);
	qc->perfect_gates = FALSE;
	qc->track = TRUE;
	qc->get_boundary = NULL;
	qc->boundaries = NULL;
	qc->logical_loss_err_pr = 0;
	qc->last_logical_loss_err_pr = 0;
	qc->logical_loss_err_du = 0;
	qc->last_logical_loss_err_du = 0;

}

/**
 * \brief Cleans up the dots and lines from any balls that exist in the nest
 * but belong to sets that have yet to be finalized. This only occurs when
 * using advanced recipes, because the dots and lines are created at the same
 * time as the set, rather than when the set is finalized.
 * 
 * \param[in] nest The \ref nest to clear of uninserted dots and lines
 */
void qc_free_uninserted_dots_and_lines(NEST *nest) {
	BALL *ball;
	DOT *dot;
	LL_NODE *ll, *ll_temp;
	LINE *line;
	CDLL_NODE *cdll;

	cdll = nest->ball_cdll->next;
	while (cdll != nest->ball_cdll) {
		ball = (BALL *)cdll->key;

		// If the ball has a value for big_t then its parent set is finalized
		// which means the ball is either converted or in finalized_balls, and
		// hence doesn't need to be cleaned up here.
		if (ball->big_t != LONG_MAX) {
			cdll = cdll->next;
			continue;
		}

		ll = ball->dot->lines;
		while (ll != NULL) {
			line = (LINE *)ll->key;
			ll_temp = ll->next;

			// Delete the line from the linked list of the destination dot.
			// This stops that dot from also trying to delete the line and
			// reading freed memory.
			dot = (line->a == ball->dot) ? line->b : line->a;
			assert(dot->merge == NULL);
			dot->lines = ll_delete_node(dot->lines, line);

			// Free the line and the node in the linked list
			free(ll);
			free(line);

			ll = ll_temp;
		}

		// Free the dot and clean up the ball's pointer
		free(ball->dot);
		ball->dot = NULL;

		cdll = cdll->next;
	}
}

/**
 * \brief Cleans up the dots and lines from any finalized balls that have yet
 * to be converted (have their dots and lines into the matching).
 *
 * Those confused will be reminded that first sets are finalized, which under
 * advanced recipes puts the ball (which already exists) into the
 * finalized_balls list. Then, after enough time has passed to consider the
 * dots and lines finalized (no chance of being merged due to loss) they are
 * themselves finalized and inserted into the matching. This function aims to
 * free lines that exist in balls that have been finalized, but the ball has
 * yet to be "converted", having its dots and lines inserted into matching.
 * 
 * \param[in] finalized_balls The finalized_balls \ref cdll_node from \ref qc
 */
void qc_free_unfinalized_dots_and_lines(CDLL_NODE *finalized_balls) {
	BALL *ball;
	DOT *dot;
	LL_NODE *ll, *ll_temp;
	LINE *line;
	CDLL_NODE *cdll;

	cdll = finalized_balls->next;
	while (cdll != finalized_balls) {
		ball = (BALL *)cdll->key;
		ll = ball->dot->lines;

		while (ll != NULL) {
			line = (LINE *)ll->key;
			ll_temp = ll->next;

			// Delete the line from the linked list of the destination dot.
			// This stops that dot from also trying to delete the line and
			// reading freed memory.
			dot = (line->a == ball->dot) ? line->b : line->a;
			assert(dot->merge == NULL);
			dot->lines = ll_delete_node(dot->lines, line);

			// Free the line and the node in the linked list
			free(ll);
			free(line);

			ll = ll_temp;
		}

		// Free the dot and clean up the ball's pointer
		free(ball->dot);
		ball->dot = NULL;

		cdll = cdll->next;
	}
}

/**
 * \brief Frees a quantum computer 
 *
 *
 * To free a copied quantum computer, see ::qc_free_qc_copy.
 *
 * \sa qc_free_qc_copy
 */
void qc_free_qc(QC *qc) {
	/*
	 * \remark Syndrome qubits are owned by the user, and are not freed here.
	 */
	CDLL_NODE *n;
	ERROR *e;

	int i;
	// Free the qubits
	cdll_free(qc->qubit_cdll, qc_free_void_qubit);

	// Syndromes are owned by the user, therefore we do no wish to free the
	// keys in the syndrome heap. 
	bh_free(qc->syn_heap, NULL); 
   
	// When using advanced recipes, balls, dots and lines are created when the
	// sets are created. As a result, the computation will end with balls at
	// the time boundary that are unfinalized and hence the dots and lines have
	// never been inserted into the matching. This loops over the nest and
	// cleans up and dots and lines in unconverted balls. This is slow, however
	// given qc_free_qc is called at the end of time, this is not an issue.
	// 
	// There are also balls that have been finalized, but are not sufficiently
	// far back in the past for their dots and lines to be finalized. That is,
	// there is still a chance for them to be merged together due to loss.
	// Hence, at the end of time, finalized balls near the current time will be
	// unconverted and their dots/lines will need to be cleaned up also.
	if (qc->recipe_adv != NULL) {
		qc_free_uninserted_dots_and_lines(qc->nest_pr);
		qc_free_uninserted_dots_and_lines(qc->nest_du);
		qc_free_unfinalized_dots_and_lines(qc->finalized_balls);
	}

	// Free the set heap
	bh_free(qc->set_heap, qc_free_void_set);

	// Free the cdll structure containing the finalized_balls, this will simply
	// be the head node when advanced recipes aren't being used.
	cdll_free(qc->finalized_balls, NULL);

	// Free the cdll structure and any zombie dots. These are dots that have
	// been merged into another dot but need to exist for the sake of
	// dot->merge.
	cdll_free(qc->zombie_dots, m_free_void_dot);

	// Free the error models
	for (i = 0; i < qc->num_ems; i++) {
		qc_free_error_model(qc->ems[i]);
	}
	free(qc->ems);

	// Free the error cdll nodes in the qc. This is not nessecarily
	// required, but it means that we can assert all errors are properly
	// removed from the simulation before being deleted.
	n = qc->error_cdll->next;
	while (n != qc->error_cdll) {
		e = (ERROR *)n->key;
		e->qc_cdlln = NULL;
		n = n->next;
	}

	cdll_free(qc->error_cdll, qc_free_void_error);
	cdll_free(qc->nest_error_cdll, qc_free_void_error);

	// Free the detection events heap and cdll
	ht_free(qc->de_ht, NULL);
	cdll_free(qc->de_cdll, qc_free_void_de);

	// Free the primal and dual nests
	qc_free_nest(qc->nest_pr);
	qc_free_nest(qc->nest_du);
	
	// Free the matchings
	m_free_matching(qc->m_pr);
	m_free_matching(qc->m_du);

	free(qc);
}

/**
 * \brief Frees a copied quantum computer and restores the previous one
 *
 * To free a non-copied quantum computer, see ::qc_free_qc.
 *
 * \param[in] qc The _copied_ quantum computer to be freed
 *
 * \sa qc_free_qc
 */
void qc_free_qc_copy(QC *qc) {
	/*
	 * \remark Syndrome qubits are owned by the user, and are not freed here.
	 * \remark We do not copy the error models, so they are not freed here
	 * \remark We do not copy the matchings, so they also are not freed here.
	 */

	// Free the qubits
	cdll_free(qc->qubit_cdll, qc_free_void_qubit);

	// Syndromes are owned by the user, therefore we do no wish to free the
	// keys in the syndrome heap. 
	bh_free(qc->syn_heap, NULL); 

	// Free the set heap
	bh_free(qc->set_heap, qc_free_void_set);
	cdll_free(qc->finalized_balls, NULL);
	cdll_free(qc->zombie_dots, NULL);

	// Free the error cdlls
	cdll_free(qc->error_cdll, qc_free_void_error);
	cdll_free(qc->nest_error_cdll, qc_free_void_error);

	// Free the detection events heap and cdll
	ht_free(qc->de_ht, NULL);
	cdll_free(qc->de_cdll, qc_free_void_de);

	// Free the primal and dual nests
	qc_free_nest(qc->nest_pr);
	qc_free_nest(qc->nest_du);

	// If there is a recipe, free it.
	if (qc->recipe) {
		my_2d_free(qc->recipe->n, (void **)qc->recipe->dotarr1);
		my_2d_free(qc->recipe->n, (void **)qc->recipe->dotarr2);
		free(qc->recipe);
	}

	if (qc->recipe_adv) {
		my_2d_free(qc->recipe_adv->n, (void **)qc->recipe_adv->block_arr);
		my_2d_free(qc->recipe_adv->n, (void **)qc->recipe_adv->t_arr);
		free(qc->recipe_adv);
	}

	free(qc);

	// Restore the random number generator to before the copy was made
	restore_rand();
}

/** 
 * \brief Increments the big_t counter and resets the time of all syndrome
 * qubits to 0. 
 *
 * \param[in,out] qc The quantum computer to increment big_t on
 */
void qc_increment_big_t(QC *qc) {
	int i;
	SYNDROME *syn;

	qc->big_t++;

	for (i = 1; i <= qc->syn_heap->num_elem; i++) {
		syn = (SYNDROME *)qc->syn_heap->k[i];
		syn->t = 0;
	}
}

/**
 * \brief Performs the minimum weight perfect matching algorithm on the \f$t = big_t\f$ time layer, for both primal and dual matchings
 *
 * \param[in,out] qc The quantum computer that contains the matchings
 * \param[in] undo 1 if the matching should be undoable, 0 otherwise.
 */
void qc_mwpm(QC *qc, int undo) {
	int d;

	if (qc->recipe) {
		d = (qc->recipe->m + 1) / 2;
	} else if (qc->recipe_adv) {
		d = (qc->recipe_adv->m + 1) / 2;
	}

	// Match the primal matching
	qc->m_pr->t = qc->big_t;
	qc->m_du->t = qc->big_t;

	if (qc->recipe) {
		qc->m_pr->D = d * qc->recipe->min_horz_wt_pr;
		qc->m_du->D = d * qc->recipe->min_horz_wt_du;
	} else if (qc->recipe_adv) {
		qc->m_pr->D = d * qc->recipe_adv->min_horz_wt_pr;
		qc->m_du->D = d * qc->recipe_adv->min_horz_wt_du;
	} else {
		qc->m_pr->D = 3 * 30000; 
		qc->m_du->D = 3 * 30000; 
	}

	// Match the dual matching
	//printf("M_PR\n");
	m_update_dots_and_lines_into_bfs(qc->m_pr, qc->m_pr->g, undo);
	m_mwpm(qc->m_pr, undo);
	//printf("M_DU\n");
	m_update_dots_and_lines_into_bfs(qc->m_du, qc->m_du->g, undo);
	m_mwpm(qc->m_du, undo);
}

/**
 * \brief Performs all of the undo steps found in the primal and dual matching
 * problems. After executing each undo step, the steps are removed.
 *
 * \param[in,out] qc The quantum computer that contains the matchings to be
 * undone.
 */
void qc_undo_mwpm(QC *qc) {
	m_execute_and_delete_undos(qc->m_pr);
	m_execute_and_delete_undos(qc->m_du);
}

/**
 * \brief Copies the state of a quantum computer into a new quantum computer. 
 *
 * \param[in] qc The quantum computer to copy
 *
 * \return A copy of the quantum computer.
 */
QC *qc_copy_qc(QC *qc) {
	QC *qc2;
	int i, j;

	// Save the current random state so that it can be undone later
	save_rand();

	qc2 = (QC *)my_malloc(sizeof(QC));
	qc2->s0 = qc->s0;
	qc2->s1 = qc->s1;
	qc2->big_t = qc->big_t;
	qc2->unfinalized_big_t = qc->unfinalized_big_t;
	qc2->num_qubits = qc->num_qubits;
	qc2->qubit_cdll = cdll_copy(qc->qubit_cdll, qc_copy_void_qubit, qc_qubit_set_qc_cdlln);
	qc2->syn_heap = bh_copy(qc->syn_heap, qc_copy_void_syndrome, NULL);
	qc2->set_heap = bh_copy(qc->set_heap, qc_copy_void_set, NULL);
	qc2->num_ems = qc->num_ems;

	// copy them. This also means they are not freed when this copy is.
	qc2->ems = qc->ems; 
	
	qc2->num_es = qc->num_es;
	qc2->next_label = qc->next_label;
	qc2->num_errors = qc->num_errors;
	qc2->error_cdll = cdll_copy(qc->error_cdll, qc_copy_void_error, qc_error_set_qc_cdlln);
	qc2->num_nest_errors = qc->num_nest_errors;
	qc2->nest_error_cdll = cdll_copy(qc->nest_error_cdll, qc_copy_void_error, NULL);
	qc2->num_des = qc->num_des;
	qc2->de_ht = ht_copy(qc->de_ht, qc_copy_void_de);
	qc2->de_cdll = cdll_copy(qc->de_cdll, qc_copy_void_de, qc_de_set_qc_cdlln);
	qc2->nest_pr = qc_copy_nest(qc->nest_pr);
	qc2->nest_du = qc_copy_nest(qc->nest_du);

	//cdll_print(qc->finalized_balls, qc_print_void_ball);
	qc2->finalized_balls = cdll_copy(qc->finalized_balls, qc_copy_null_ball, NULL);
	qc2->zombie_dots = cdll_copy(qc->zombie_dots, qc_copy_null_dot, NULL);

	// The matchings are undone at the end of correct rather than deleted, as
	// such we do not deep copy them here. This also means they are not freed
	// when this copy is.
	qc2->m_pr = qc->m_pr; 
	qc2->m_du = qc->m_du;

	// Perform a deep copy of the recipe if there is one available.
	if (qc->recipe) {
		qc2->recipe = qc_create_recipe(qc->recipe->type, qc->recipe->size, qc->recipe->n, qc->recipe->m);
		qc2->recipe->num_layers = qc->recipe->num_layers;
		free(qc2->recipe->layers);
		qc2->recipe->layers = qc->recipe->layers;
		qc2->recipe->blocks = qc->recipe->blocks;
		qc2->recipe->offsets = qc->recipe->offsets;
		qc2->recipe->repeat_id = qc->recipe->repeat_id;
		qc2->recipe->repeated_count = qc->recipe->repeated_count;
		qc2->recipe->min_horz_wt_pr = qc->recipe->min_horz_wt_pr;	
		qc2->recipe->min_horz_wt_du = qc->recipe->min_horz_wt_du;	
		
		for (i = 0; i < qc->recipe->n; ++i) {
			for (j = 0; j < qc->recipe->m; ++j) {
				qc2->recipe->dotarr1[i][j] = qc->recipe->dotarr1[i][j];
				qc2->recipe->dotarr2[i][j] = qc->recipe->dotarr2[i][j];
			}
		}
	} else {
		qc2->recipe = NULL;
	}

	if (qc->recipe_adv) {
		qc2->recipe_adv = qc_create_recipe_adv(qc->recipe_adv->type, qc->recipe_adv->n, qc->recipe_adv->m, TRUE);

		for (i = 0; i < qc->recipe_adv->n; ++i) {
			for (j = 0; j < qc->recipe_adv->m; ++j) {
				qc2->recipe_adv->block_arr[i][j] = qc->recipe_adv->block_arr[i][j];
				qc2->recipe_adv->t_arr[i][j] = qc->recipe_adv->t_arr[i][j];
			}
		}

		qc2->recipe_adv->min_horz_wt_pr = qc->recipe_adv->min_horz_wt_pr;	
		qc2->recipe_adv->min_horz_wt_du = qc->recipe_adv->min_horz_wt_du;	

		qc2->recipe_adv->in_cycle = qc->recipe_adv->in_cycle;
		qc2->recipe_adv->cycle_t0 = qc->recipe_adv->cycle_t0;
		qc2->recipe_adv->cycle_period = qc->recipe_adv->cycle_period;

		qc2->recipe_adv->block_ht = qc->recipe_adv->block_ht;
		qc2->recipe_adv->offset_ht = qc->recipe_adv->offset_ht;
		qc2->recipe_adv->block_cdll = qc->recipe_adv->block_cdll;

		qc2->recipe_adv->cycle_len = qc->recipe_adv->cycle_len;
		qc2->recipe_adv->cycle_begin = qc->recipe_adv->cycle_begin;
		qc2->recipe_adv->cycle_end = qc->recipe_adv->cycle_end;
	} else {
		qc2->recipe_adv = NULL;
	}
	
	qc2->perfect_gates = qc->perfect_gates;
	qc2->track = qc->track;
	qc2->get_boundary = qc->get_boundary;
	qc2->boundaries = qc->boundaries;

	qc2->logical_loss_err_pr = qc->logical_loss_err_pr;
	qc2->last_logical_loss_err_pr = qc->last_logical_loss_err_pr;
	qc2->logical_loss_err_du = qc->logical_loss_err_du;
	qc2->last_logical_loss_err_du = qc->last_logical_loss_err_du;

	return qc2;
}

/**
 * \ingroup qubit
 * \brief Copies a qubit 
 *
 * \param[in] q The qubit to be copies
 *
 * \sa qc_copy_void_qubit
 *
 * \returns A copy of the qubit 
 */
QUBIT *qc_copy_qubit(QUBIT *q) {
	/** 
	 * \remark q2->qc_cdlln does not need to be set explictly, it is set by cdll_copy
	 */

	QUBIT *q2;
	CDLL_NODE *n;
	ERROR *error;

	// If a copy of this qubit already exists, then return the stored copy
	if (q->copy != NULL) {
		return q->copy;
	}

	q2 = (QUBIT *)my_malloc(sizeof(QUBIT));

	// Point the original's copy to the newly created copy qubit
	q->copy = q2;

	q2->i = q->i;
	q2->j = q->j;
	q2->k = q->k;

	q2->t = q->t;
	q2->e = q->e;
	q2->num_errors = q->num_errors;

	// Copy the error cdll
	q2->error_cdll = cdll_copy(q->error_cdll, qc_copy_void_error, qc_error_set_q_cdlln);
	
	// Traverse the copied error cdll and update the qubit they belong to to be
	// the newly created copy qubit.
	n = q2->error_cdll->next;
	while (n != q2->error_cdll) {
		error = (ERROR *)n->key;
		error->q = q2;
		n = n->next;
	}

	q2->error_ht = ht_copy(q->error_ht, qc_copy_void_error);

	// Point newly created copy qubit's copy to the original qubit
	q2->copy = q;

	return q2;
}

/**
 * \ingroup qubit
 * \brief Copies a qubit. Accepts a void pointer to a qubit.
 *
 * This function is to be used when it needs to be passed generically to a
 * function.
 *
 * \param[in] key The qubit to be copied
 *
 * \sa qc_copy_qubit
 *
 * \return A void pointer to a copy of the qubit
 */
void *qc_copy_void_qubit(void *key) {
	return (void *)qc_copy_qubit((QUBIT *)key);
}

/**
 * \ingroup error
 * \brief Copies an error. 
 *
 * \param[in] error The error to be copied
 *
 * \sa qc_copy_void_error
 *
 * \return A copy of the error
 */
ERROR *qc_copy_error(ERROR *error) {
	/**
	 * \remark error2->qc_cdlln will be set by ::qc_copy_qc via ::cdll_copy when needed
	 * \remark error2->q_cdlln will be set by ::qc_copy_qubit via ::cdll_copy when needed
	 * \remark error2->q will be set by ::qc_copy_qubit when needed
	 */

	ERROR *error2;

	// If there is already a copy of this error, return that copy
	if (error->copy != NULL) {
		return error->copy;
	}

	error2 = (ERROR *)my_malloc(sizeof(ERROR));
	error->copy = error2;

	error2->big_t = error->big_t;
	error2->label = error->label;
	error2->type = error->type;
	error2->op = error->op;
	error2->p = error->p;

	error2->gate = error->gate;
	error2->i = error->i;
	error2->j = error->j;
	error2->k = error->k;
	error2->t = error->t;

	error2->qc_cdlln = NULL;
	error2->q_cdlln = NULL;
	error2->q = NULL;

	error2->stick = qc_copy_stick(error->stick);
	error2->next = qc_copy_error(error->next);
	error2->prev = qc_copy_error(error->prev);
	error2->copy = error;

	return error2;
}

/**
 * \ingroup error
 * \brief Copies an error. Accepts a void pointer to an error.
 *
 * This function is to be used when it needs to be passed generically to a
 * function.
 *
 * \param[in] key The error to be copied
 *
 * \sa qc_copy_error
 *
 * \return A void pointer to a copy of the error
 */
void *qc_copy_void_error(void *key) {
	return (void *)qc_copy_error((ERROR *)key);
}

/**
 * \ingroup stick
 * \brief Copies a stick.
 *
 * \param[in] stick The stick to be copied
 *
 * \sa qc_copy_void_stick
 *
 * \return A copy of the stick 
 */
STICK *qc_copy_stick(STICK *stick) {
	/**
	 * \remark stick2->nest_cdlln is always set by ::qc_copy_nest
	 */

	STICK *stick2;

	// If there is no stick, return
	if (stick == NULL) {
		return NULL;
	}

	// If there is already a copy of the stick, return that copy
	if (stick->copy != NULL) {
		return stick->copy;
	}

	stick2 = (STICK *)my_malloc(sizeof(STICK));
	stick->copy = stick2;

	stick2->p_stick = stick->p_stick;
	stick2->a = qc_copy_ball(stick->a);
	stick2->b = qc_copy_ball(stick->b);
	stick2->num_errors = stick->num_errors;
	stick2->error_ll = ll_copy(stick->error_ll, qc_copy_void_error);
	stick2->copy = stick;

	return stick2;
}

/**
 * \ingroup stick
 * \brief Copies a stick. Accepts a void pointer to a stick. 
 *
 * This function is to be used when it needs to be passed generically to a
 * function.
 *
 * \param[in] key The stick to be copied
 *
 * \sa qc_copy_stick
 *
 * \return A void pointer to a copy of the stick
 */
void *qc_copy_void_stick(void *key) {
	return (void *)qc_copy_stick((STICK *)key);
}

/**
 * \ingroup ball
 * \brief Copies a ball.
 *
 * \param[in] ball The ball to be copied
 *
 * \sa qc_copy_void_ball
 *
 * \return A copy of the ball
 */
BALL *qc_copy_ball(BALL *ball) {
	/**
	 * \remark Dots are not copied, so the dot in the ball copy is the same as the
	 * dot in the original ball
	 */

	BALL *ball2;

	// If there is no ball, return
	if (ball == NULL) {
		return NULL;
	}

	// If there is already a copy of the ball, return that copy
	if (ball->copy != NULL) {
		return ball->copy;
	}

	// Do not copy boundary balls
	if (ball->type == PRIMAL_BOUNDARY || ball->type == DUAL_BOUNDARY) {
		return ball;
	}

	ball2 = (BALL *)my_malloc(sizeof(BALL));
	ball->copy = ball2;

	ball2->type = ball->type;
	ball2->i = ball->i;
	ball2->j = ball->j;
	ball2->t = ball->t;
	ball2->big_t = ball->big_t;
	ball2->stick_ht = ht_copy(ball->stick_ht, qc_copy_void_stick);
	ball2->mp = ball->mp;
	ball2->hash_ijt = ball->hash_ijt;

	ball2->dot = ball->dot; 
	ball2->copy = ball;

	return ball2;
}

void *qc_copy_null_ball(void *key) {
	BALL *ball;

	ball = (BALL *)key;
	if (ball->copy == NULL) {
		printf("Attempting to use copy of ball that doesn't exist.\n");
		exit(1);
	}
	return (void *)ball->copy;
}

void *qc_copy_null_dot(void *key) {
	return key;
}

/**
 * \ingroup ball
 * \brief Copies a ball. Accepts a void pointer to a ball. 
 *
 * This function is to be used when it needs to be passed generically to a
 * function.
 *
 * \param[in] key The ball to be copied
 *
 * \sa qc_copy_ball
 *
 * \return A void pointer to a copy of the ball
 */
void *qc_copy_void_ball(void *key) {
	return (void *)qc_copy_ball((BALL *)key);
}

/**
 * \ingroup syndrome
 * \brief Copies a syndrome qubit. 
 *
 * \param[in] syn The syndrome qubit to be copied
 *
 * \sa qc_copy_void_syndrome
 *
 * \return A copy of the syndrome qubit.
 */
SYNDROME *qc_copy_syndrome(SYNDROME *syn) {
	/**
	 * \remark syn2->set is always set by ::qc_copy_syn_cdll
	 * \remark syn2->set_cdlln is always set by ::qc_copy_syn_cdll
	 */

	SYNDROME *syn2;

	// If there is already a copy of this syndrome qubit, return that copy
	if (syn->copy != NULL) {
		return syn->copy;
	}

	syn2 = (SYNDROME *)my_malloc(sizeof(SYNDROME));
	syn->copy = syn2;

	syn2->t = syn->t;
	syn2->qc_heap_i = syn->qc_heap_i;
	syn2->copy = syn;

	syn2->set = NULL;
	syn2->set_cdlln = NULL;

	return syn2;
}

/**
 * \ingroup syndrome
 * \brief Copies a syndrome qubit. Accepts a void pointer to a syndrome qubit. 
 *
 * This function is to be used when it needs to be passed generically to a
 * function.
 *
 * \param[in] key The syndrome qubit to be copied
 *
 * \sa qc_copy_syndrome
 *
 * \return A void pointer to a copy of the syndrome qubit
 */
void *qc_copy_void_syndrome(void *key) {
	return (void *)qc_copy_syndrome((SYNDROME *)key);
}

/**
 * \ingroup set
 * \brief Copies a set.
 *
 * \param[in] set The set to be copied
 *
 * \sa qc_copy_void_set
 *
 * \return A copy of the set
 */
SET *qc_copy_set(SET *set) {
	/**
	 * \remark Boundaries are not copied. Therefore, the boundary in the set is
	 * not a copy. 
	 */
	SET *set2;

	// If there is already a copy of this set, return that copy.
	if (set->copy != NULL) {
		return set->copy;
	}

	set2 = (SET *)my_malloc(sizeof(SET));
	set->copy = set2;

	set2->set_pos_cdll = cdll_copy(set->set_pos_cdll, qc_copy_void_set_pos, NULL);

	set2->type = set->type;
	set2->i = set->i;
	set2->j = set->j;
	set2->t = set->t;
	set2->big_t = set->big_t;
	set2->qc_heap_i = set->qc_heap_i;
	set2->num_meas_left = set->num_meas_left;
	set2->bdy = set->bdy;
	set2->syn_cdll = qc_copy_syn_cdll(set->syn_cdll);
	set2->num_errors = set->num_errors;
	set2->mt_cdll = cdll_copy(set->mt_cdll, qc_copy_void_measurement, NULL);
	set2->ball = qc_copy_ball(set->ball);
	set2->block = set->block;
	set2->copy = set;

	return set2;
}

/**
 * \ingroup set
 * \brief Copies a set. Accepts a void pointer to a set. 
 *
 * This function is to be used when it needs to be passed generically to a
 * function.
 *
 * \param[in] key The set to be copied
 *
 * \sa qc_copy_set
 *
 * \return A void pointer to a copy of the set
 */
void *qc_copy_void_set(void *key) {
	return (void *)qc_copy_set((SET *)key);
}

/**
 * \ingroup syndrome
 * \brief Copies a syndrome circular double-linked list 
 *
 * \param[in] syn_cdll The sentinel of the circular double-linked list to be
 * copied
 *
 * \return The sentinel to the copy of the circular double-linked list 
 */
CDLL_NODE *qc_copy_syn_cdll(CDLL_NODE *syn_cdll) {
	CDLL_NODE *syn_cdll2, *n;
	SYNDROME *syn, *syn2;

	syn_cdll2 = cdll_create();
	
	// Traverse the list, copy the syndrome and its set, build a new list
	n = syn_cdll->prev;
	while (n != syn_cdll) {
		syn = (SYNDROME *)n->key;
		syn2 = qc_copy_syndrome(syn);
		syn2->set = qc_copy_set(syn->set);
		cdll_insert_head(syn_cdll2, syn2, qc_syndrome_set_set_cdlln);
		n = n->prev;
	}

	return syn_cdll2;
}

/**
 * \ingroup measurement
 * \brief Copies a measurement
 *
 * \param[in] mt The measurement to be copied
 *
 * \sa qc_copy_void_measurement
 *
 * \return A copy of the measurement 
 */
MEASUREMENT *qc_copy_measurement(MEASUREMENT *mt) {
	MEASUREMENT *mt2;

	// If there is already a copy of this measurement, return that copy.
	if (mt->copy != NULL) {
		return mt->copy;
	}

	mt2 = (MEASUREMENT *)my_malloc(sizeof(MEASUREMENT));
	mt->copy = mt2;

	mt2->m = mt->m;
	mt2->num_parent_sets = mt->num_parent_sets;
	mt2->num_errors = mt->num_errors;
	mt2->error_cdll = cdll_copy(mt->error_cdll, qc_copy_void_error, NULL);
	mt2->copy = mt;

	return mt2;
}


/**
 * \ingroup measurement
 * \brief Copies a measurement. Accepts a void pointer to a measurement. 
 *
 * This function is to be used when it needs to be passed generically to a
 * function.
 *
 * \param[in] key The measurement to be copied
 *
 * \sa qc_copy_measurement
 *
 * \return A void pointer to a copy of the measurement
 */
void *qc_copy_void_measurement(void *key) {
	return (void *)qc_copy_measurement((MEASUREMENT *)key);
}

SET_POS *qc_copy_set_pos(SET_POS *sp) {
	SET_POS *sp2;

	sp2 = (SET_POS *)my_malloc(sizeof(SET_POS));
	sp2->i = sp->i;
	sp2->j = sp->j;
	sp2->lay = sp->lay;

	return sp2;
}

void *qc_copy_void_set_pos(void *key) {
	return (void *)qc_copy_set_pos((SET_POS *)key);
}

/**
 * \ingroup de
 * \brief Copies a detection event 
 *
 * \param[in] de The detection event to be copied
 *
 * \sa qc_copy_void_de
 *
 * \return A copy of the detection event
 */
DE *qc_copy_de(DE *de) {
	/**
	 * \remark de2->qc_cdlln always set by ::qc_copy_de_cdll
	 */
	DE *de2;

	// If there is already a copy of this detection event, return that copy.
	if (de->copy != NULL) {
		return de->copy;
	}

	de2 = qc_create_de(qc_copy_set(de->set), qc_copy_error(de->error));
	de->copy = de2;
	de2->copy = de;

	return de2;
}

/**
 * \ingroup de
 * \brief Copies a detection event. Accepts a void pointer to a detection event. 
 *
 * This function is to be used when it needs to be passed generically to a
 * function.
 *
 * \param[in] key The detection event to be copied
 *
 * \sa qc_copy_de
 *
 * \return A void pointer to a copy of the detection event
 */
void *qc_copy_void_de(void *key) {
	return (void *)qc_copy_de((DE *)key);
}

/**
 * \ingroup nest
 * \brief Copies a nest. 
 *
 * \param[in] nest The nest to be copied
 *
 * \return A copy of the nest 
 */
NEST *qc_copy_nest(NEST *nest) {
	NEST *nest2;
	BALL *ball;
	CDLL_NODE *n;
	STICK *stick2;
	int ip;

	// If there is already a copy of this nest, return that copy.
	if (nest->copy != NULL) {
		return nest->copy;
	}

	nest2 = (NEST *)my_malloc(sizeof(NEST));
	nest->copy = nest2;

	nest2->ball_cdll = cdll_copy(nest->ball_cdll, qc_copy_void_ball, qc_ball_set_nest_cdlln);

	// If the last_converted_ball_cdlln is equal to the ball_cdll then it just
	// points to the sentinel of the list, we can just set the copied
	// last_converted_ball_cdlln to be the sentinel of the newly copied
	// ball_cdll
	if (nest->last_converted_ball_cdlln == nest->ball_cdll) {
		nest2->last_converted_ball_cdlln = nest2->ball_cdll;
	}

	// Otherwise, we need to find the copied equivalent of the
	// last_converted_ball_cdlln. We can do this by following the ball of the
	// original, going to the copy, then taking its nest_cdlln, which will have
	// been set during the copying of the ball_cdlln.
	else {
		ball = (BALL *)nest->last_converted_ball_cdlln->key;
		if (ball) {
			assert(ball->copy != NULL);
			nest2->last_converted_ball_cdlln = ball->copy->nest_cdlln;
		}
	}

	// The logic here works the same as above.
	if (nest->last_blocked_ball_cdlln == nest->ball_cdll) {
		nest2->last_blocked_ball_cdlln = nest2->ball_cdll;
	}
	else {
		ball = (BALL *)nest->last_blocked_ball_cdlln->key;
		assert(ball->copy != NULL);
		nest2->last_blocked_ball_cdlln = ball->copy->nest_cdlln;
	}

	nest2->stick_cdll = cdll_copy(nest->stick_cdll, qc_copy_void_stick, qc_stick_set_nest_cdlln);

	// Copy the stick hash table by creating a new one then traversing the
	// stick cdll and inserting the elements. 
	nest2->stick_ht = ht_create(nest->stick_ht->length);
	n = nest2->stick_cdll->next;
	while (n != nest2->stick_cdll) {
		stick2 = (STICK *)n->key;
		ip = ((size_t)stick2->a + (size_t)stick2->b)%INT_MAX;
		ht_insert_key(nest2->stick_ht, ip, stick2);
		n = n->next;
	}

	nest2->copy = nest;

	return nest2;
}



// QUBIT FUNCTIONS

/**
 * \addtogroup qubit
 * @{ 
 */

/**
 * \brief Creates a new qubit
 *
 * \param[in] i The i coordinate of the \ref qubit
 * \param[in] j The j coordinate of the \ref qubit
 * \param[in] k The k coordinate of the \ref qubit
 * \param[in] ht_size The size of the \ref error \ref ht 
 *
 * \return The newly created \ref qubit
 */
QUBIT *qc_create_qubit(int i, int j, int k, int ht_size) {
	QUBIT *q;

	q = (QUBIT *)my_malloc(sizeof(QUBIT));

	q->t = 0;
	q->i = i;
	q->j = j;
	q->k = k;
	q->e = 0;
	q->num_errors = 0;
	q->error_cdll = cdll_create();
	q->error_ht = ht_create(ht_size);
	q->qc_cdlln = NULL;
	q->copy = NULL;

	return q;
}

/**
 * \brief Frees a qubit
 *
 * \param[in] q The \ref qubit to be freed 
 *
 * \sa qc_free_void_qubit 
 */
void qc_free_qubit(QUBIT *q) {
	/**
	 * \remark The errors in the qubit are owned by qc, and are hence not freed by this function. 
	 */
	CDLL_NODE *n;
	ERROR *e;

	// The errors are not being freed, but this qubit is. Ensure that no errors
	// still point back to the qubit being deleted. They should become
	// unaccessible without this qubit until freed from the list of errors,
	// however this allows us to assert that errors have been cleaned up
	// properly before being freed.
	n = q->error_cdll->next;
	while (n != q->error_cdll) {
		e = (ERROR *)n->key;
		e->q_cdlln = NULL;
		n = n->next;
	}

	cdll_free(q->error_cdll, NULL);
	ht_free(q->error_ht, NULL);
	
	// If there is a copy, NULL the void pointer pointing back to this qubit.
	if (q->copy != NULL) {
		q->copy->copy = NULL;
	}
	
	free(q);
}

/**
 * \brief Frees a \ref qubit. Accepts a void pointer to a \ref qubit.
 *
 * \param[in] key The \ref qubit to be freed 
 *
 * \sa qc_free_qubit 
 */

void qc_free_void_qubit(void *key) {
	qc_free_qubit((QUBIT *)key);
}

/**
 * \brief Inserts a \ref qubit into a \ref qc quantum computer
 *
 * \param[out] qc The \ref qc in which to insert the \ref qubit
 * \param[in] q The \ref qubit to insert into the \ref qc
 */
void qc_insert_qubit(QC *qc, QUBIT *q) {
	cdll_insert_head(qc->qubit_cdll, q, qc_qubit_set_qc_cdlln);
	qc->num_qubits++;
}

/**
 * \brief Links a \ref qubit to the \ref cdll node in the \ref qc
 *
 * \param[out] key The \ref qubit as a void pointer
 * \param[in] n The \ref cdll node to link to the \ref qubit
 */
void qc_qubit_set_qc_cdlln(void *key, CDLL_NODE *n) {
	((QUBIT *)key)->qc_cdlln = n;
}

/**
 * \brief Creates and inserts a new \ref qubit into a \ref qc
 *
 * An alias for ::qc_create_qubit then ::qc_insert_qubit
 *
 * \param[out] qc The \ref qc in which to insert the \ref qubit
 * \param[in] i The i coordinate of the \ref qubit
 * \param[in] j The j coordinate of the \ref qubit
 * \param[in] k The k coordinate of the \ref qubit
 * \param[in] ht_size The size of the \ref error \ref ht 
 *
 * \sa qc_create_qubit, qc_insert_qubit
 *
 * \return The newly created and inserted qubit
 */
QUBIT *qc_create_and_insert_qubit(QC *qc, int i, int j, int k, int ht_size) {
	QUBIT *q;

	q = qc_create_qubit(i, j, k, ht_size);
	qc_insert_qubit(qc, q);

	return q;
}

/**
 * \brief Deletes and frees a \ref qubit from a \ref qc 
 *
 * \param[out] qc The \ref qc from which the \ref qubit will be deleted
 * \param[in] q The \ref qubit to be deleted
 */
void qc_delete_qubit(QC *qc, QUBIT *q) {
	cdll_delete_node(q->qc_cdlln, qc_free_void_qubit);
	qc->num_qubits--;
}

/**
 * \brief Initialises a new qubit 
 *
 * \param[in,out] qc The \ref qc 
 * \param[in,out] q The \ref qubit to initialise
 * \param[in] gate The gate to perform the initialisation
 * \param[in] p The probability of error 
 * \param[in] multiply_es A functon that takes two es and multiplies them. 
 *
 * \sa qc_init_qubit_invisible
 */
void qc_init_qubit(QC *qc, QUBIT *q, int gate, double p,
	int (*multiply_es)(int e1, int e2)) {
	ERROR_MODEL *em;
	CDLL_NODE *n, *tn;

	assert(gate >= 0);

	// If error tracking is on, delete all errors currently in the qubit
	if (qc->track) {
		n = q->error_cdll->next;
		while (n != q->error_cdll) {
			tn = n;
			n = n->next;
			qc_delete_error(qc, (ERROR *)tn->key);
		}
	}

	q->e = 0;

	// Get the error model for the gate
	em = qc->ems[gate];

	// Set the duration of the qubit to the duration from the error model
	q->t += em->duration;

	// Simulate a random raw error
	qc_random_raw_e(qc, q, em, p, multiply_es);

	// Track all possible single qubit errors
	if (qc->track) {
		qc_add_errors(qc, q, gate, p);
	}
}

/**
 * Initialises a new qubit, but doesn't generate errors
 *
 * \param[in,out] qc The \ref qc 
 * \param[in,out] q The \ref qubit to initialise
 * \param[in] gate The gate to perform the initialisation
 * \param[in] p The probability of error 
 * \param[in] multiply_es A functon that takes two es and multiplies them. 
 *
 * \sa qc_init_qubit
 */
void qc_init_qubit_invisible(QC *qc, QUBIT *q, int gate, double p,
	int (*multiply_es)(int e1, int e2)) {
	ERROR_MODEL *em;

	assert(gate >= 0);

	q->e = 0;

	// Get the error model for the gate
	em = qc->ems[gate];

	// Simulate a random raw error
	qc_random_raw_e(qc, q, em, p, multiply_es);
}

/**
 * \brief Performs a single qubit unitary gate on a \ref qubit
 *
 * \param[in,out] qc The \ref qc 
 * \param[in,out] q The \ref qubit to perform the gate on
 * \param[in] gate The gate to use
 * \param[in] p The probability of error 
 * \param[in] transform_e A function that takes a \ref qubit and transforms the
 * e appropriate to the gate
 * \param[in] multiply_es A functon that takes two es and multiplies them. 
 * \param[in] transform_error A function that takes an \ref error and
 * transforms if appropriate to the gate
 */
void qc_single_qubit_unitary(QC *qc, QUBIT *q, int gate, double p,
	void (*transform_e)(QUBIT *q),
	int (*multiply_es)(int e1, int e2),
	void (*transform_error)(ERROR *error)) {
	ERROR_MODEL *em;
	CDLL_NODE *n;
	ERROR *error;

	assert(gate >= 0);

	em = qc->ems[gate];
	q->t += em->duration;

	// Simulate a random raw error
	transform_e(q);
	qc_random_raw_e(qc, q, em, p, multiply_es);

	// Track all possible single-qubit errors
	if (qc->track) {
		n = q->error_cdll->next;
		while (n != q->error_cdll) {
			error = (ERROR *)n->key;
			transform_error(error);
			n = n->next;
		}
		qc_add_errors(qc, q, gate, p);
	}
}

/**
 * \brief Performs a two-qubit gate on two seperate \ref qubit
 *
 * \param[in,out] qc The \ref qc 
 * \param[in,out] q1 The first \ref qubit to perform the gate on
 * \param[in,out] q2 The second \ref qubit to perform the gate on
 * \param[in] gate The gate to use
 * \param[in] p The probability of error 
 * \param[in] transform_es A function that takes two \ref qubit and transforms their 
 * e appropriate to the gate
 * \param[in] multiply_es A functon that takes two es and multiplies them. 
 * \param[in] transform_errors A function that takes two error operators and transforms them 
 */
void qc_two_qubit_unitary(QC *qc, QUBIT *q1, QUBIT *q2, int gate, double p,
	void (*transform_es)(QUBIT *q1, QUBIT *q2),
	int (*multiply_es)(int e1, int e2),
	void (*transform_errors)(int *op1, int *op2)) {
	ERROR_MODEL *em;
	int op1, op2;
	DLL_NODE *dlln;
	CDLL_NODE *temp_cdll, *cdlln;
	ERROR *error1, *error2;

	assert(gate >= 0);

	// We should never apply a two qubit unitary gate across the time dimension
	if (q1->t != q2->t) {
		printf("q1(%d, %d, %d, %ld) q2(%d, %d, %d, %ld)\n", q1->i, q1->j, q1->k, q1->t, q2->i, q2->j, q2->k, q2->t);
		assert(q1->t == q2->t);
	}

	// Get the gate error model and set the qubit durations
	em = qc->ems[gate];
	q1->t += em->duration;
	q2->t += em->duration;

	// Simulate a random raw error
	transform_es(q1, q2);
	qc_random_raw_e2(qc, q1, q2, em, p, multiply_es);

	if (qc->track) {
		// Move second qubit errors to temporary list, we want any new errors
		// to be at the head of the list, so it is easiest to join the old
		// errors to the new errors once we are finished here.
		temp_cdll = q2->error_cdll;
		q2->error_cdll = cdll_create();

		// Loop the errors on the first qubit, we want to apply errors in the
		// qubit1 -> qubit2 direction.
		cdlln = q1->error_cdll->next;
		while (cdlln != q1->error_cdll) {
			error1 = (ERROR *)cdlln->key;
			cdlln = cdlln->next;
			op1 = error1->op;

			// If the error from the first qubit is also on the second qubit,
			// then we want to apply the error
			dlln = ht_find_node_fn(q2->error_ht, error1->label, error1,
				qc_error_eq_label);
			if (dlln != NULL) {
				error2 = (ERROR *)dlln->key;
				op2 = error2->op;
			}
			else {
				error2 = NULL;
				op2 = I;
			}
	
			// Apply the two-qubit transformation
			transform_errors(&op1, &op2);

			if (op2 == I) {
				// If the error on qubit2 is now (or always was) I, then delete the error
				if (error2 != NULL) {
					qc_delete_error(qc, error2);
				}
			}
			else {
				// If there was already an error on qubit2, and it is
				// non-identity, then update the error on the second qubit
				// and insert it into the cdll of the appropriate qubit.
				if (error2 != NULL) {
					error2->op = op2;

					// We previously stored q2->error_cdll in temp_cdll, and
					// then created a new cdll in its place. We therefore need
					// to move the node from that old list to this new list.
					cdll_extract_node(error2->q_cdlln);
					cdll_insert_node_head(q2->error_cdll, error2->q_cdlln);
				}

				// If there was no error on the second qubit, and we now have a non
				// identity error, then we want to create and insert a new error.
				else {
					error2 = qc_create_error(error1->i, error1->j, error1->k, error1->t, 
						error1->big_t, error1->label, op2, error1->gate, error1->p);
					qc_insert_error_nearby(qc, q2, error2, error1);
				}
			}

			if (op1 == I) {
				// If the error on qubit1 is now (or always was) I, then delete the error
				qc_delete_error(qc, error1);
			}
			else {
				// Update the error on the first qubit
				error1->op = op1;
			}
		}

		// Loop over the errors previously on qubit2. This list now excludes
		// any errors extracted above. We want to apply errors in the qubit2 ->
		// qubit1 direction.
		cdlln = temp_cdll->next;
		while (cdlln != temp_cdll) {
			op1 = I;

			error2 = (ERROR *)cdlln->key;
			cdlln = cdlln->next;
			assert(error2 != NULL);
			op2 = error2->op;

			// Here we can assume that qubit1 doesn't have the error, if it
			// did, it would have been performed in the qubit1 -> qubit2
			// transformations earlier

			// Apply the two-qubit transformation
			transform_errors(&op1, &op2);

			// If the error on the first qubit is non-identity, create and
			// insert a new error
			if (op1 != I) {
				error1 = qc_create_error(error2->i, error2->j, error2->k, error2->t, 
					error2->big_t, error2->label, op1, error2->gate, error2->p);
				qc_insert_error_nearby(qc, q1, error1, error2);
			}

			if (op2 == I) {
				// If the error on qubit2 is now (or was always) I, then delete
				// the error
				qc_delete_error(qc, error2);
			}
			else {
				// Update the error on the second qubit
				error2->op = op2;
			}
		}

		// Join the old qubit2 errors back
		cdll_join(q2->error_cdll, temp_cdll);

		qc_add_errors2(qc, q1, q2, gate, p);
	}
}

/**
 * \brief Measures a qubit
 * 
 * \param[in] qc The \ref qc the \ref qubit is in
 * \param[in] q The \ref qubit to be measured
 * \param[in] gate The \ref gate to measure the \ref qubit
 * \param[in] basis The basis in which to measure the \ref qubit
 * \param[in] p The probability of a random error occuring 
 * \param[in] set1 The first \ref set
 * \param[in] set2 The second \ref set 
 * \param[in] multiply_es A function that takes two es and multiplies them
 * \param[in] meas_sim A function that takes a \ref qubit and simulates a
 * measurement on the qubit
 */
int qc_meas_qubit(QC *qc, QUBIT *q, int gate, int basis, double p,
	SET *set1, SET *set2,
	int (*multiply_es)(int e1, int e2),
	int (*meas_sim)(QUBIT *q)) {
	ERROR_MODEL *em;
	int m;
	SET *set;
	MEASUREMENT *mt;
	int i;
	CDLL_NODE *n, *tn;
	ERROR *error;

	#ifdef DBG_MEAS_QUBIT 
	printf("qc_meas_qubit()\n");
	qc_print_set(set1);
	qc_print_set(set2);
	#endif

	if (set1->big_t != LONG_MAX || set2->big_t != LONG_MAX) {
		qc_print_set(set1);
		qc_print_set(set2);
		assert(set1->big_t == LONG_MAX && set2->big_t == LONG_MAX);
	}

	// When using advanced recipes, sets have their dots and balls created
	// immediately, so they already have a t and big_t set at time of creation
	if (!qc->recipe_adv) {
		assert(set1->big_t == LONG_MAX);
		assert(set2->big_t == LONG_MAX);
	}

	assert(gate >= 0);
	if (qc->track) {
		assert(q->num_errors == q->error_ht->num_elem);
	}

	// Get the error model for the gate and set the duration
	em = qc->ems[gate];
	q->t += em->duration;

	// Simulate a random error
	qc_random_raw_e(qc, q, em, p, multiply_es);

	// Simulate measurement on the qubit
	m = meas_sim(q);

	// Track all possible single errors
	if (qc->track) {
		qc_add_errors(qc, q, gate, p);
		assert(q->num_errors == q->error_ht->num_elem);
	}

	// Update how many measurements are required to complete a set
	if (!qc_is_boundary_set(set1)) {
		set1->num_meas_left--;
	}
	if (!qc_is_boundary_set(set2)) {
		set2->num_meas_left--;
	}	
	
	// The qubit measured has been lost
	if (m == 0) {
		// Delete all errors from the qubit
		if (qc->track) {
			n = q->error_cdll->next;
			while (n != q->error_cdll) {
				tn = n;
				n = n->next;
				qc_delete_error(qc, (ERROR *)tn->key);
			}
			assert(q->error_ht->num_elem == 0);
			assert(q->num_errors == q->error_ht->num_elem);
			assert(n->next == n);
			assert(n->prev == n);
		}

		// If we are attempting to merge a set to itself, do nothing.
		if (set1 != set2) {
			set = qc_merge_sets(qc, set1, set2);

			// The set will be NULL if a logical loss error has occured
			// resulting in no merge operation.
			if (set != NULL && set->num_meas_left == 0) {
				qc_finalize_set(qc, set, q->t);
			}
		}
	}

	// The qubit measured has a valid measurement result
	else {
		// If there are any errors on the qubit, we want to associate them with
		// the measurement rather than the qubit. 
		if (qc->track) {
			// Filter the errors on the qubit 
			n = q->error_cdll->next;
			while (n != q->error_cdll) {
				tn = n;
				n = n->next;
				error = (ERROR *)tn->key;

				// If we are measuring in the X basis, and there is an X
				// component of the error operation, then delete the error
				if (basis == X && (error->op & Z) != Z) {
					qc_delete_error(qc, error);
				}
				
				// If we are measuring in the Z basis, and there is an Z
				// component of the error operation, then delete the error
				else if (basis == Z && (error->op & X) != X) {
					qc_delete_error(qc, error);
				}

				// The error is no longer associated with a qubit, rather with measurement
				else {
					error->q = NULL;
				}
			}

			assert(q->num_errors == q->error_ht->num_elem);

			// Empty the qubit error hash table
			for (i = 0; i < q->error_ht->length; i++) {
				dll_free(q->error_ht->table[i], NULL);
				q->error_ht->table[i] = NULL;
			}
			q->error_ht->num_elem = 0;
		}

		// Create a new measurement
		mt = qc_create_measurement(m);
		if (qc->track) {
			// Copy the errors from the qubit into the measurement
			tn = mt->error_cdll;
			mt->error_cdll = q->error_cdll;
			q->error_cdll = tn;
			mt->num_errors = q->num_errors;
			q->num_errors = 0;
		}

		// Insert measurement into the sets (if not a boundary)
		// 
		// If the sets are the same, then this measurement is interior to an
		// agglomerated set, and hence should not be considered.
		if (set1 != set2) {
			if (!qc_is_boundary_set(set1)) {
				qc_insert_measurement(set1, mt);
			}
			if (!qc_is_boundary_set(set2)) {
				qc_insert_measurement(set2, mt);
			}
		}

		// If the measurement hasn't been inserted into any sets (i.e if it joined two boundary sets), 
		// it needs to be freed now.
		if (mt->num_parent_sets == 0) {
			qc_free_measurement(mt);
		}

		// Finalize the sets if appropriate
		if (!qc_is_boundary_set(set1) && set1->num_meas_left == 0) {
			qc_finalize_set(qc, set1, q->t);
		}
		if (set1 != set2 && !qc_is_boundary_set(set2) && set2->num_meas_left == 0) {
			qc_finalize_set(qc, set2, q->t);
		}
	}

	return m;
}

/** @} */

// ERROR MODEL FUNCTIONS

/** \addtogroup error_model
 * @{
 */

/**
 * \brief Reads an error model from a file and creates a new error model 
 *
 * \param[in] qc The quantum computer that the error model will be used in
 * \param[in] fname The filename of the file containing the error model details
 *
 * \returns The error model created from the file.
 */
ERROR_MODEL *qc_create_error_model(const char *fname) {
	int i;
	long int total;															  
	FILE *in;														   
	ERROR_MODEL *em;

	// Open the file for reading
	in = fopen(fname, "r");											 
	if (in == NULL) {
		return NULL;
	}
																	   
	em = (ERROR_MODEL *)my_malloc(sizeof(ERROR_MODEL));

	// Read in the number of qubits, the scale and the number of lines in the
	// error model
	fscanf(in, "%d\n", &em->num_qubits);							
	fscanf(in, "%lf\n", &em->scale);
	fscanf(in, "%d\n", &em->num_lines);							 
											 
	// Errors should only be 1 or 2 qubit errors
	assert(em->num_qubits == 1 || em->num_qubits == 2);

	total = 0;
	if (em->num_lines > 0) {
		em->raw_em = (long int **)my_2d_calloc(em->num_lines,			   
			em->num_qubits+1, sizeof(long int));							  

		// If we are dealing with a single qubit gate, there is only one error
		// to be applied, so read in the file format accordingly.
		if (em->num_qubits == 1) {
			for (i = 0; i < em->num_lines; i++) {
				fscanf(in, "%ld %ld\n", &em->raw_em[i][0], &em->raw_em[i][1]);

				// Sum the total relative strengths so that they can be
				// converted to probabilities
				total += em->raw_em[i][0];
			}
		}
		else {
			for (i = 0; i < em->num_lines; i++) {
				fscanf(in, "%ld %ld %ld\n", &em->raw_em[i][0], &em->raw_em[i][1], &em->raw_em[i][2]);
				total += em->raw_em[i][0];
			}
		}
	}
	else {
		// No error model, must be a perfect gate
		em->raw_em = NULL;
	}
	
	// Read in the gate duration
	fscanf(in, "%d\n", &em->duration);	

	// Close the file for reading
	fclose(in);

	// If there is no error model, there is no probability of any errors
	if (em->raw_em == NULL) {
		em->raw_rel_p = NULL;
		em->sum_raw_rel_p = NULL;
		return em;
	}

	em->raw_rel_p = (double *)my_malloc(em->num_lines*sizeof(double));
	em->sum_raw_rel_p = (double *)my_malloc(em->num_lines*sizeof(double));

	// Convert the relative strengths from the file into probabilities
	em->raw_rel_p[0] = (double)em->raw_em[0][0]/total;
	em->sum_raw_rel_p[0] = em->raw_rel_p[0];
	for (i = 1; i < em->num_lines; i++) {
		em->raw_rel_p[i] = (double)em->raw_em[i][0]/total;
		em->sum_raw_rel_p[i] = em->sum_raw_rel_p[i-1] + em->raw_rel_p[i];
	}
	
	return em;
}

/**
 * \brief Frees an error model 
 *
 * \param[in] em The error model to be freed 
 */
void qc_free_error_model(ERROR_MODEL *em) {
	if (em == NULL) return;

	if (em->num_lines > 0) my_2d_free(em->num_lines, (void **)em->raw_em);
	if (em->raw_rel_p != NULL) free(em->raw_rel_p);
	if (em->sum_raw_rel_p != NULL) free(em->sum_raw_rel_p);
	free(em);
}

/**
 * \brief Inserts an error model into the quantum computer structure
 *
 * \param[out] qc The quantum computer to be inserted into
 * \param[in] em The error model to be inserted 
 *
 * \return The index of the error model inserted
 */
int qc_insert_error_model(QC *qc, ERROR_MODEL *em) {
	if (em == NULL) return -1;

	qc->num_ems++;
	qc->ems = (ERROR_MODEL **)my_realloc(qc->ems, qc->num_ems*sizeof(ERROR_MODEL *));
	qc->ems[qc->num_ems-1] = em;

	return qc->num_ems-1;
}

/**
 * \brief Perform a random error on a single qubit 
 *
 * \param[in] qc The \ref qc
 * \param[out] q The \ref qubit to perform the random error on
 * \param[in] em The \ref error_model of the random error to apply
 * \param[in] p The probability of a random error ocuring
 * \param[in] multiply_es A function that takes two es and multiplies them
 */
void qc_random_raw_e(QC *qc, QUBIT *q, ERROR_MODEL *em, double p,
	int (*multiply_es)(int e1, int e2)) {
	int i;
	unsigned int u;

	assert(em->num_qubits == 1);
	
	// If gates are perfect, or the probability of error is 0, do nothing
	if (qc->perfect_gates == TRUE || p == 0 || em->num_lines == 0 || em->scale == 0) {
		return;
	}

	// Generate a random unsigned integer 
	u = randomui();

	// Scale the probability and multiply by UINT_MAX to get the probability in
	// the range 0 to UINT_MAX.
	p *= em->scale;
	p *= UINT_MAX;
	
	// If the random number is greater than or equal to this, then no error
	// occured, do nothing
	if (u >= p) { 
		return;
	}

	// At this point, an error has occured. Loop over each possible error in
	// the error model and apply the error if it randomly occured
	qc->num_es++;
	for (i = 0; i < em->num_lines; i++) {

		// em->sum_raw_rel_p[i] is a rolling sum of the probabilities from each
		// error, such that em->sum_raw_rel_p[len-1] equals 1. This ensures
		// that at least one of the errors will be applied.
		if (u < p * em->sum_raw_rel_p[i]) {
			//printf("**RNDERR** i: %d, j: %d, k: %d, t: %ld:%ld, op: %ld\n", q->i, q->j, q->k, q->t/CYCLE_T, q->t%CYCLE_T, em->raw_em[i][1]);
			q->e = multiply_es(q->e, em->raw_em[i][1]);
		//	printf("i: %d, j: %d, k: %d, t: %ld:%ld, op: %ld\n", q->i, q->j, q->k, q->t/CYCLE_T, q->t%CYCLE_T, em->raw_em[i][1]);
			return;
		}
	}
}

/**
 * \brief Perform a random error on two qubits
 * 
 * \param[in] qc The \ref qc
 * \param[out] q1 The 1st \ref qubit to perform the random error on
 * \param[out] q2 The 2nd \ref qubit to perform the random error on
 * \param[in] em The \ref error_model of the random error to apply
 * \param[in] p The probability of a random error ocuring
 * \param[in] multiply_es A function that takes two es and multiplies them
 */
void qc_random_raw_e2(QC *qc, QUBIT *q1, QUBIT *q2, ERROR_MODEL *em, double p,
	int (*multiply_es)(int e1, int e2)) {
	int i;
	unsigned int u;

	assert(em->num_qubits == 2);
	
	// If gates are perfect, or the probability of error is 0, do nothing
	if (qc->perfect_gates == TRUE || p == 0 || em->num_lines == 0 || em->scale == 0) {
		return;
	}

	u = randomui();
	p *= em->scale;
	p *= UINT_MAX;
	
	// No error has occured, do nothing
	if (u >= p) { 
		return;
	}

	// At this point, an error has occured. Loop over each possible error in
	// the error model and apply the error if it randomly occured
	qc->num_es++;
	for (i = 0; i < em->num_lines; i++) {
		if (u < p * em->sum_raw_rel_p[i]) {
			//printf("**RNDERR2** i1: %d, j1: %d, k1: %d, i2, %d, j2: %d, k2: %d, t: %ld:%ld, op1: %ld, op2: %ld\n", q1->i, q1->j, q1->k, q2->i, q2->j, q2->k, q1->t/CYCLE_T, q1->t%CYCLE_T, em->raw_em[i][1], em->raw_em[i][2]);
			q1->e = multiply_es(q1->e, em->raw_em[i][1]);
			q2->e = multiply_es(q2->e, em->raw_em[i][2]);
		//	printf("i1: %d, j1: %d, k1: %d, i2, %d, j2: %d, k2: %d, t: %ld:%ld, op1: %ld, op2: %ld\n", q1->i, q1->j, q1->k, q2->i, q2->j, q2->k, q1->t/CYCLE_T, q1->t%CYCLE_T, em->raw_em[i][1], em->raw_em[i][2]);
			return;
		}
	}
}

/**
 * \brief Add errors for a single qubit gate 
 *
 * \param[in] qc The quantum computer
 * \param[in] q The qubit to add the errors from
 * \param[in] gate The gate that has been applied to the qubit
 * \param[in] p The probability of a error occuring
 *
 * \sa qc_add_errors2
 */
void qc_add_errors(QC *qc, QUBIT *q, int gate, double p) {
	int i, e;
	double p_e;
	ERROR *error, sentinel;
	ERROR_MODEL *em;

	// Get the correct error model for this gate
	em = qc->ems[gate];

	// Assert that we are dealing with a single qubit gate and that we have the
	// correct number of errors in the qubit and its hash table at the start
	assert(em->num_qubits == 1);
	assert(q->num_errors == q->error_ht->num_elem);

	// If our gates are perfect, or we have no chance of error, simply return
	if (qc->perfect_gates == TRUE || p == 0 || em->raw_rel_p == NULL) return;

	sentinel.next = &sentinel;
	sentinel.prev = &sentinel;

	// Multiply p by the error model scaling factor so that all the
	// probabilities are correct when multiplied.
	p *= em->scale;

	// For each possible error in the error model, add the error to the computation
	for (i = 0; i < em->num_lines; i++) {

		// The probability of any error occuring multiplied by the
		// probability of this specific error occuring
		p_e = p * em->raw_rel_p[i];

		if (p_e != 0) {
			e = em->raw_em[i][1];

			// Assert that the error is an X, Y or Z Pauli error
			assert(e != I);
			assert(e == X || e == Y || e == Z || e == L);

			error = qc_create_error(q->i, q->j, q->k, q->t, qc->big_t, qc->next_label, e, gate, p_e);
			qc_insert_error(qc, q, error, &sentinel);
			qc->next_label++;
		}
	}

	sentinel.next->prev = sentinel.prev;
	sentinel.prev->next = sentinel.next;

	// Assert that after adding the errors, there is still consistency between
	// the number of errors on the qubit and its hash table
	assert(q->num_errors == q->error_ht->num_elem);
}

/**
 * \brief Add errors for a two-qubit gate 
 *
 * \param[in] qc The quantum computer
 * \param[in] q1 The first qubit to add the errors from
 * \param[in] q2 The second qubit to add the errors from
 * \param[in] gate The gate that has been applied to the qubit
 * \param[in] p The probability of the error occuring
 *
 * \sa qc_add_errors
 */
void qc_add_errors2(QC *qc, QUBIT *q1, QUBIT *q2, int gate, double p) {
	int i, e;
	double p_e;
	ERROR *error, sentinel;
	ERROR_MODEL *em;

	// Get the correct error model for the gate
	em = qc->ems[gate];

	// Assert that we are dealing with a two qubit gate and that the number of
	// errors on each of the qubits is consistent with their respective hash
	// tables
	assert(em->num_qubits == 2);
	assert(q1->num_errors == q1->error_ht->num_elem);
	assert(q2->num_errors == q2->error_ht->num_elem);

	// If our gates are perfect, or we have no chance of error, simply return
	if (qc->perfect_gates == TRUE || p == 0 || em->raw_rel_p == NULL) return;

	sentinel.next = &sentinel;
	sentinel.prev = &sentinel;

	// Multiply p by the error model scaling factor so that all the
	// probabilities are correct when multiplied.
	p *= em->scale;

	// For each possible error in the error model, add the error to the computation
	for (i=0; i<em->num_lines; i++) {
		
		// The probability of any error occuring multiplied by the
		// probability of this specific error occuring
		p_e = p * em->raw_rel_p[i];
		
		if (p_e != 0) {
			e = em->raw_em[i][1];
			
			// Assert that the error is an X, Y or Z Pauli error
			assert(e == I || e == X || e == Y || e == Z);
			if (e != I) {
				error = qc_create_error(q1->i, q1->j, q1->k, q1->t, qc->big_t, qc->next_label, e, gate, p_e);
				qc_insert_error(qc, q1, error, &sentinel);
			}
			e = em->raw_em[i][2];
			assert(e == I || e == X || e == Y || e == Z);
			if (e != I) {
				error = qc_create_error(q2->i, q2->j, q2->k, q2->t, qc->big_t, qc->next_label, e, gate, p_e);
				qc_insert_error(qc, q2, error, &sentinel);
			}
			qc->next_label++;
		}
	}

	sentinel.next->prev = sentinel.prev;
	sentinel.prev->next = sentinel.next;

	assert(q1->num_errors == q1->error_ht->num_elem);
	assert(q2->num_errors == q2->error_ht->num_elem);
}

/** @} */





// ERROR FUNCTIONS

/**
 * \addtogroup error 
 * @{ 
 */

/**
 * \brief Creates a new \ref error
 * 
 * \param[in] i The i coordinate of the error 
 * \param[in] j The j coordinate of the error
 * \param[in] k The k coordinate of the error   
 * \param[in] t The time of the error 
 * \param[in] big_t The big_t of the error
 * \param[in] label The label of the error (the ID)
 * \param[in] op The operator of the error (I, X, Y, Z)
 * \param[in] gate The gate that caused the error
 * \param[in] p The probability that the error will occur
 *
 * \return The newly created \ref error
 */
ERROR *qc_create_error(int i, int j, int k, long int t, long int big_t, int label, int op, int gate, double p) {
	ERROR *error;

	error = (ERROR *)my_malloc(sizeof(ERROR));

	error->i = i;
	error->j = j;
	error->k = k;
	error->t = t;
	error->big_t = big_t;
	error->label = label;
	error->type = UNKNOWN;
	error->op = op;
	error->gate = gate;
	error->p = p;
	error->qc_cdlln = NULL;
	error->q_cdlln = NULL;
	error->q = NULL;
	error->stick = NULL;
	error->next = NULL;
	error->prev = NULL;
	error->copy = NULL;

	return error;
}

/**
 * \brief Frees an \ref error 
 * 
 * \param[in] error The \ref error to be freed
 *
 * \sa qc_free_void_error
 */
void qc_free_error(ERROR *error) {
	// Ensure that the error has been entirely cleaned up before being freed.
	// We do not want any links to this error left in qc or any qubits after it
	// is freed.
	if (error->qc_cdlln != NULL) {
		qc_print_error(error);
		assert(error->qc_cdlln == NULL);
	}
	if (error->q_cdlln != NULL) {
		qc_print_error(error);
		assert(error->q_cdlln == NULL);
	}

	// If there is a copy of this error, nullify the pointer from that copy to
	// this error
	if (error->copy != NULL) { 
		error->copy->copy = NULL;
	}
	
	free(error);
}

/**
 * \brief Frees an error. Accepts a void pointer to an \ref error.
 * 
 * \param[in] key The \ref error to be freed.
 *
 * \sa qc_free_error
 */
void qc_free_void_error(void *key) {
	qc_free_error((ERROR *)key);
}

/**
 * \brief NULL the qubit of an \ref error. Accepts a void pointer to an \ref error.
 * 
 * \param[in] key The \ref error to be freed
 */
void qc_null_error_q(void *key) {
	((ERROR *)key)->q = NULL;
}

/**
 * \brief Set the uplink from an \ref error to the \ref cdll in the \ref qc
 * 
 * \param[in] key The \ref error to be linked
 * \param[in] n The ::CDLL_NODE in the \ref qc to uplink to. 
 */
void qc_error_set_qc_cdlln(void *key, CDLL_NODE *n) {
	((ERROR *)key)->qc_cdlln = n;
}

/**
 * \brief Set the uplink from an \ref error to the \ref cdll in the \ref qubit
 * 
 * \param[in] key The \ref error to be linked
 * \param[in] n The ::CDLL_NODE in the \ref qubit to uplink to. 
 */
void qc_error_set_q_cdlln(void *key, CDLL_NODE *n) {
	((ERROR *)key)->q_cdlln = n;
}

/**
 * \brief Inserts an \ref error into a \ref qubit and the \ref qc 
 * 
 * \param[in] qc The \ref qc to be inserted into
 * \param[in] q The \ref qubit to be inserted into
 * \param[in] error The \ref error to be inserted
 * \param[in] sentinel The sentinel of the error list
 */
void qc_insert_error(QC *qc, QUBIT *q, ERROR *error, ERROR *sentinel) {
	assert(error != sentinel);
	assert(q->num_errors == q->error_ht->num_elem);

	// Insert into the quantum computer
	cdll_insert_head(qc->error_cdll, error, qc_error_set_qc_cdlln);
	qc->num_errors++;
	
	// Insert into the qubit
	cdll_insert_head(q->error_cdll, error, qc_error_set_q_cdlln);
	ht_insert_key(q->error_ht, error->label, error);
	error->q = q;
	q->num_errors++;

	error->next = sentinel->next;
	error->prev = sentinel;
	sentinel->next->prev = error;
	sentinel->next = error;

	assert(q->num_errors == q->error_ht->num_elem);
}

/**
 * \brief Inserts an \ref error into a \ref qubit and the \ref qc. The \ref
 * error is inserted next the \ref qc next to an existing \ref error. 
 *
 * This function is slightly more efficient for putting related errors into the
 * \ref qc.
 * 
 * \param[in] qc The \ref qc to be inserted into
 * \param[in] q The \ref qubit to be inserted into
 * \param[in] error The \ref error to be inserted
 * \param[in] sentinel The sentinel of the error list
 */
void qc_insert_error_nearby(QC *qc, QUBIT *q, ERROR *error, ERROR *sentinel) {
	assert(error != sentinel);
	assert(sentinel->qc_cdlln != NULL);
	assert(q->num_errors == q->error_ht->num_elem);

	// Insert into the quantum computer
	cdll_insert_head(sentinel->qc_cdlln, error, qc_error_set_qc_cdlln);
	qc->num_errors++;

	// Insert into the qubit
	cdll_insert_head(q->error_cdll, error, qc_error_set_q_cdlln);
	ht_insert_key(q->error_ht, error->label, error);
	error->q = q;
	q->num_errors++;

	error->next = sentinel->next;
	error->prev = sentinel;
	sentinel->next->prev = error;
	sentinel->next = error;

	assert(q->num_errors == q->error_ht->num_elem);
}

/**
 * \brief Inserts an \ref error into a \ref ht if the error does not exist in
 * the \ref ht. Otherwise, if the \ref error is already in the \ref ht, remove it.
 * 
 * \param[in] ht The \ref ht to insert/delete from
 * \param[in] error The \ref error to insert/delete
 */
void qc_mod_insert_error_in_ht(HT *ht, ERROR *error) {
	int i;
	DLL_NODE *n;
	ERROR *err2;

	// Determine the bin in the hash table by modulus on the ID by the size
	i = error->label % ht->length;
	n = ht->table[i];

	// If there is something in the dll, search the dll for the error
	while (n != NULL) {
		err2 = (ERROR *)n->key;

		// If the error is already in the hash table, remove it and returns
		if (err2->label == error->label) {
			ht->table[i] = dll_delete_node(ht->table[i], n, NULL);
			ht->num_elem--;
			return;
		}
		
		n = n->next;
	}

	// If the error was not found, insert it
	ht->table[i] = dll_insert(ht->table[i], error, NULL);
	ht->num_elem++;
}

/**
 * \brief Deletes and frees an \ref error
 *
 * Will delete and free an \ref error from both the associated \ref qc and \ref qubit.
 * 
 * \param[in] qc The \ref qc to have the error deleted from
 * \param[in] error The \ref error to be deleted 
 */
void qc_delete_error(QC *qc, ERROR *error) {
	// Remove the error from the group of errors 
	error->next->prev = error->prev;
	error->prev->next = error->next;

	// Delete the error from the quantum computer
	assert(error->qc_cdlln);
	cdll_delete_node(error->qc_cdlln, NULL);
	qc->num_errors--;

	error->qc_cdlln = NULL;

	// Delete the error from the qubit
	assert(error->q_cdlln);
	cdll_delete_node(error->q_cdlln, NULL);
	ht_delete_key(error->q->error_ht, error->label, error, NULL);
	error->q->num_errors--;

	error->q_cdlln = NULL;

	// Free the error
	qc_free_error(error);
}

/**
 * \brief Compares the labels of two \link error errors\endlink for equality
 * 
 * \param[in] key1 The first \ref error to be compared
 * \param[in] key2 The second \ref error to be compared
 *
 * \return 1 if the two \link error errors\endlink have the same label, 0 otherwise.
 */
int qc_error_eq_label(void *key1, void *key2) {
	return ((ERROR *)key1)->label == ((ERROR *)key2)->label;
}

/**
 * \brief Moves an \ref error into a \ref nest
 * 
 * \param[in] qc The \ref qc containing the \ref error
 * \param[in] error The \ref error to move 
 */
void qc_move_error_to_nest(QC *qc, ERROR *error) {
	CDLL_NODE *n;

	assert(error != NULL);

	// Remove the error from the quantum computer
	n = error->qc_cdlln;
	n->next->prev = n->prev;
	n->prev->next = n->next;
	qc->num_errors--;

	// Remove the error from the qubit that created it
	if (error->q != NULL) {
		assert(error->q->num_errors == error->q->error_ht->num_elem);
		assert(error->q_cdlln != NULL);

		cdll_delete_node(error->q_cdlln, NULL);
		error->q_cdlln = NULL;
		ht_delete_key(error->q->error_ht, error->label, error, NULL);
		error->q->num_errors--;
		assert(error->q->num_errors == error->q->error_ht->num_elem);
		error->q = NULL;
	}

	// Add the error to the quantum computer nest list
	n->next = qc->nest_error_cdll->next;
	n->prev = qc->nest_error_cdll;
	qc->nest_error_cdll->next->prev = n;
	qc->nest_error_cdll->next = n;
	qc->num_nest_errors++;

	// Prevents additional detection events occurring
	error->qc_cdlln = NULL;
}

/** @} */













// DE FUNCTIONS

/**
 * \addtogroup de 
 * @{
 */

/**
 * \brief Creates a new \ref de
 * 
 * \param[in] set The \ref set for the \ref de
 * \param[in] error The \ref error for the \ref de 
 *
 * \return The newly created \ref de
 */
DE *qc_create_de(SET *set, ERROR *error) {
	DE *de;

	de = (DE *)my_malloc(sizeof(DE));

	de->set = set;
	de->error = error;

	// Done for faster access (no need to dereference every time)
	de->label = error->label; 
	
	de->qc_cdlln = NULL;
	de->copy = NULL;

	return de;
}

/**
 * \brief Frees a \ref de
 * 
 * \param[in] de The \ref de to be freed
 */
void qc_free_de(DE *de) {
	// If there is a copy of the detection event, nullify the pointer from the
	// copy to this detection event
	if (de->copy != NULL)  {
		de->copy->copy = NULL;
	}
	
	free(de);
}

/**
 * \brief Frees a \ref de. Accepts a void pointer to a \ref de.
 * 
 * \param[in] key The \ref de to be freed
 */
void qc_free_void_de(void *key) {
	qc_free_de((DE *)key);
}

/**
 * \brief Set the uplink from an \ref de to the \ref cdll in the \ref qc
 * 
 * \param[in] key The \ref de to be linked
 * \param[in] n The ::CDLL_NODE in the \ref qc to uplink to. 
 */
void qc_de_set_qc_cdlln(void *key, CDLL_NODE *n) {
	((DE *)key)->qc_cdlln = n;
}

/**
 * \brief Inserts a \ref de into a \ref qc.
 * 
 * \param[in] qc The \ref qc to insert the \ref de into
 * \param[in] de The \ref de to be inserted. 
 */
void qc_insert_de(QC *qc, DE *de) {
	ht_insert_key(qc->de_ht, de->label, de);
	cdll_insert_head(qc->de_cdll, de, qc_de_set_qc_cdlln);
	qc->num_des++;
}

/**
 * \brief Deletes and frees a \ref de from a \ref qc.
 * 
 * \param[in] qc The \ref qc to delete the \ref de from
 * \param[in] de The \ref de to be deleted 
 */
void qc_delete_de(QC *qc, DE *de) {
	ht_delete_key(qc->de_ht, de->label, de, NULL);
	cdll_delete_node(de->qc_cdlln, qc_free_void_de);
	qc->num_des--;
}

/**
 * \brief Compares the labels of two \link de Detection Events\endlink for equality
 * 
 * \param[in] key1 The first \ref de to be compared
 * \param[in] key2 The second \ref de to be compared
 *
 * \return 1 if the two \link de Detection Events\endlink have the same label, 0 otherwise.
 */
int qc_de_eq_label(void *key1, void *key2) {
	return ((DE *)key1)->label == ((DE *)key2)->label;
}

/**
 * \brief Locate a \ref de by label and type
 * 
 * \param[in] qc The \ref qc to search for the \ref de in
 * \param[in] label The label of the \ref de to find
 * \param[in] type The type of the \ref de to find
 *
 * \return The \ref de if found, otherwise NULL
 */
DE *qc_find_de(QC *qc, int label, int type) {
	HT *ht;
	int i;
	DLL_NODE *n;
	DE *de;

	ht = qc->de_ht;
	i = label % ht->length;
	n = ht->table[i];

	// Search the hash table for the detection event
	while (n != NULL) {
		de = (DE *)n->key;
		if (de->label == label && de->set->type == type) {
			return de;
		}
		n = n->next;
	}

	return NULL;
}

/**
 * \brief Finalizes a \ref de
 *
 * Inserts the \ref error from the \ref de into the \ref stick that links
 * the \link de Detection Event's\endlink \ref set to its boundary.
 * 
 * \param[in] qc The \ref qc containing the \ref de
 * \param[in] de The \ref de to finalize 
 */
void qc_finalize_de(QC *qc, DE *de) {
	SET *set, *bdy;
	BALL *a, *b;
	STICK *stick;

	assert(de != NULL);
	set = de->set;
	bdy = set->bdy;

	// A finalised detection event should always have a boundary, error.
	if (bdy == NULL) {
		printf("qc->big_t: %ld\n", qc->big_t);
		qc_print_set(de->set);
		qc_print_error(de->error);
		printf("next_label: %ld, FAIL_LABEL: %d\n", qc->next_label, FAIL_LABEL);
		assert(bdy != NULL);
	}

	// Both the detection event set, and the boundary, should have balls
	a = set->ball;
	b = bdy->ball;
	assert(a != NULL && b != NULL);

	// Find the stick between the set's ball and the boundary's ball, if it
	// doesn't exist, then create it.
	stick = qc_find_stick(qc, a, b);
	if (stick == NULL) {
		stick = qc_create_stick(a, b);
		if (set->type == PRIMAL) {
			qc_insert_stick(qc->nest_pr, stick);
		}
		else {
			qc_insert_stick(qc->nest_du, stick);
		}
	}

	// Track the error in the stick
	stick->error_ll = ll_insert(stick->error_ll, de->error);
	de->error->stick = stick;
	stick->num_errors++;
}

/** @} */




// SYNDROME FUNCTIONS

/**
 * \addtogroup syndrome
 * @{
 */

/**
 * \brief Creates a new \ref syndrome
 * 
 * \return The newly created Syndrome
 */
SYNDROME *qc_create_syndrome() {
	SYNDROME *syn;

	syn = (SYNDROME *)my_malloc(sizeof(SYNDROME));

	syn->t = 0;
	syn->qc_heap_i = 0;
	syn->set = NULL;
	syn->set_cdlln = NULL;
	syn->copy = NULL;

	return syn;
}

/**
 * \brief Frees a \ref syndrome
 * 
 * \param[in] syn The \ref syndrome to be freed
 */
void qc_free_syndrome(SYNDROME *syn) {
	// If there is no syndrome, do nothing
	if (syn == NULL) {
		return;
	}
	
	// If there is a copy of this syndrome, nullify the pointer from the copy
	// back to this syndrome.
	if (syn->copy != NULL) {
		syn->copy->copy = NULL;
	}

	free(syn);
}

/**
 * \brief Inserts a \ref syndrome into a \ref qc
 * 
 * \param[in] qc The \ref qc to insert the \ref syndrome into
 * \param[in] syn The \ref syndrome to insert 
 */
void qc_insert_syndrome(QC *qc, SYNDROME *syn) {
	bh_insert_key(qc->syn_heap, syn, qc_syndrome_lt_t, qc_syndrome_swap,
		qc_syndrome_set_qc_heap_i, NULL);
}

/**
 * \brief Removes a \ref syndrome from a \ref qc
 *
 * \param[in] qc The \ref qc to remove the syndrome for
 * \param[in] syn The \ref syndrome to remove 
 */
void qc_uninsert_syndrome(QC *qc, SYNDROME *syn) {
	/**
	 * \remark This function does not free the \ref syndrome
	 */
	bh_remove_key(qc->syn_heap, syn->qc_heap_i, qc_syndrome_lt_t, qc_syndrome_swap,
		qc_syndrome_set_qc_heap_i, NULL);
}

/**
 * \brief Compares two \link syndrome Syndrome Qubits\endlink to see if the
 * first is less than the second, based on time.
 * 
 * \param[in] key1 The first \ref syndrome to be compared
 * \param[in] key2 The second \ref syndrome to be compared
 *
 * \return 1 if key1's t is less than key2's t, 0 otherwise.
 */
int qc_syndrome_lt_t(void *key1, void *key2) {
	return ((SYNDROME *)key1)->t < ((SYNDROME *)key2)->t;
}

/**
 * \brief Swaps two \link syndrome Syndrome Qubits\endlink inside a heap
 * 
 * \param[in] h The heap that both \link syndrome Syndrome Qubits\endlink are
 * in
 * \param[in] i1 The index in the heap of the first \ref syndrome
 * \param[in] i2 The index in the heap of the second \ref syndrome 
 */
void qc_syndrome_swap(BHEAP *h, int i1, int i2) {
	void *key;
	SYNDROME *syn1, *syn2;
	int i;

	key = h->k[i1];
	h->k[i1] = h->k[i2];
	h->k[i2] = key;

	syn1 = (SYNDROME *)h->k[i1];
	syn2 = (SYNDROME *)h->k[i2];

	i = syn1->qc_heap_i;
	syn1->qc_heap_i = syn2->qc_heap_i;
	syn2->qc_heap_i = i;
}

/**
 * \brief Sets the uplink from a \ref syndrome to the \ref heap in a \ref qc.
 * 
 * \param[in] key The \ref syndrome to be linked
 * \param[in] i The index in the \ref heap in the \ref qc for the \ref syndrome 
 */
void qc_syndrome_set_qc_heap_i(void *key, int i) {
	((SYNDROME *)key)->qc_heap_i = i;
}

/**
 * \brief Associates a \ref set with a \ref syndrome
 * 
 * \param[in] set The \ref set to be associated
 * \param[in] syn The \ref syndrome to be associated 
 */
void qc_associate_syndrome(SET *set, SYNDROME *syn) {
	// If there is no set, no syndrome or the set is a boundary, do nothing
	if (set == NULL || syn == NULL || set->type == PRIMAL_BOUNDARY || set->type == DUAL_BOUNDARY) {
		return;
	}

	cdll_insert_head(set->syn_cdll, syn, qc_syndrome_set_set_cdlln);
	syn->set = set;
}

/**
 * \brief Unassociates a \ref set from a \ref syndrome
 * 
 * \param[in] syn The \ref syndrome to be unassociated
 */
void qc_unassociate_syndrome(SYNDROME *syn) {
	// If there is no syndrome, or no associated set, do nothing
	if (syn == NULL || syn->set_cdlln == NULL) { 
		return;
	}

	cdll_delete_node(syn->set_cdlln, NULL);
	syn->set = NULL;
	syn->set_cdlln = NULL;
}

/**
 * \brief Sets the uplink from a \ref syndrome to the \ref cdll in a \ref set 
 * 
 * \param[in] key The \ref syndrome to be linked
 * \param[in] n The ::CDLL_NODE in the \ref set \ref cdll 
 */
void qc_syndrome_set_set_cdlln(void *key, CDLL_NODE *n) {
	((SYNDROME *)key)->set_cdlln = n;
}

/** @} */












// SET FUNCTIONS

/** \addtogroup set
 * @{
 */

/**
 * \brief Creates a new \ref set
 * 
 * \param[in] type The type of \ref set (PRIMAL, DUAL, PRIMAL_BOUNDARY, DUAL_BOUNDARY)
 * \param[in] i The i coordinate of the \ref set
 * \param[in] j The j coordinate of the \ref set
 * \param[in] num_meas_left The number of measurements left for the \ref set
 * \param[in] bdy The nearest boundary to the \ref set
 *
 * \return The newly created \ref set
 */
SET *qc_create_set(int type, int i, int j, int num_meas_left, SET *bdy) {
	/**
	 * \remark If type is a boundary, a \ref ball will also be created and linked to the \ref set.
	 */
	SET *set;

	set = (SET *)my_malloc(sizeof(SET));

	set->set_pos_cdll = cdll_create();

	set->type = type;
	set->i = i;
	set->j = j;
	set->big_t = LONG_MAX;
	set->t = LONG_MAX;
	set->qc_heap_i = 0;
	set->num_meas_left = num_meas_left;
	set->bdy = bdy;
	set->syn_cdll = cdll_create();
	set->num_errors = 0;
	set->mt_cdll = cdll_create();
	if (type == PRIMAL_BOUNDARY || type == DUAL_BOUNDARY) {
		set->ball = qc_create_ball(set, 0);
		set->ball->mp = 1;
	}
	else {
		set->ball = NULL;
	}
	set->copy = NULL;
	set->block = NULL;

	#ifdef DBG_CREATE_SET
	printf("Create Set: ");
	qc_print_set(set);
	#endif

	return set;
}

SET *qc_create_set_adv(QC *qc, int type, int i, int j, int num_meas_left, SET *bdy) {
	/**
	 * \remark If type is a boundary, a \ref ball will also be created and linked to the \ref set.
	 */
	SET *set;
	BLOCK *block;
	NEST *nest;
	MATCHING *matching;
	//DOT *dot;

	set = (SET *)my_malloc(sizeof(SET));

	set->type = type;
	set->i = i;
	set->j = j;
	set->big_t = LONG_MAX;
	set->t = LONG_MAX;
	set->qc_heap_i = 0;
	set->num_meas_left = num_meas_left;
	set->bdy = bdy;
	set->block = NULL;
	set->set_pos_cdll = cdll_create();
	set->syn_cdll = cdll_create();
	set->num_errors = 0;
	set->mt_cdll = cdll_create();
	set->ball = NULL;

	if (type == PRIMAL_BOUNDARY || type == DUAL_BOUNDARY) {
		set->ball = qc_create_ball(set, 0);
		set->ball->mp = 1;
	}
	else {
		if (qc->recipe_adv && qc->recipe_adv->cycle_len > 0) {
			//qc_print_set(set);
			set->block = qc_get_block_for_set(qc, set);

			block = (BLOCK *)qc->recipe_adv->block_arr[set->i][set->j]->prev->key;

			if (set->block->t > qc->recipe_adv->t_arr[set->i][set->j]) {
				qc->recipe_adv->t_arr[set->i][set->j] = set->block->t;
			} else {
				qc->recipe_adv->t_arr[set->i][set->j] += set->block->t - block->t;
			}

			//qc_print_block(set->block);
			set->t = qc->recipe_adv->t_arr[set->i][set->j];
			//printf("--%d %d %ld\n", set->i, set->j, set->t);
			set->ball = qc_create_ball_adv(set, STICK_HT_SIZE);
			//qc_print_ball(set->ball);

			matching = (set->type == PRIMAL) ? qc->m_pr : qc->m_du;
			nest = (set->type == PRIMAL) ? qc->nest_pr : qc->nest_du;

			qc_insert_ball(nest, set->ball);
			qc_convert_block_to_dot_and_lines(qc, matching, set->ball, set->block);
			//m_print_dot(set->ball->dot);

			//dot = m_create_and_insert_dot(matching, set->ball);
		}
	}
	set->copy = NULL;

	#ifdef DBG_CREATE_SET
	printf("Create Set: ");
	qc_print_set(set);
	#endif

	return set;
}

SET_POS *qc_create_set_pos(int i, int j, int lay) {
	SET_POS *set_pos;
	
	set_pos = (SET_POS *)my_malloc(sizeof(SET_POS));
	set_pos->i = i;
	set_pos->j = j;
	set_pos->lay = lay;

	return set_pos;
}

void qc_swap_set_layer(SET *set, int lay) {
	CDLL_NODE *n;
	SET_POS *sp;

	n = set->set_pos_cdll->next;
	while (n != set->set_pos_cdll) {
		sp = (SET_POS *)n->key;
		n = n->next;
		if (sp->lay != lay) {
			sp->lay = !sp->lay;
		} else {
			//cdll_delete_node(n->prev, free);
		}
	}
}

/**
 * \brief Frees a \ref set.
 * 
 * Will free a \ref set and linked \link measurement Measurements\endlink where
 * this is the last \ref set to which they belong.
 * If the \ref set is a boundary, it will free the \ref ball also.
 *
 * \param[in] set The \ref set to be freed.
 *
 * \sa qc_free_void_set
 */
void qc_free_set(SET *set) {
	/**
	 * \remark Sets DO NOT own \link ball Balls\endlink nests do.
	 */
	CDLL_NODE *n, *tn;
	MEASUREMENT *mt;

	// At the moment we do not want to accidentally delete boundary sets. To
	// ensure this, boundaries cannot be deleted with free_set and must be
	// manually freed using qc_free_bdy. Once loss development is over, these
	// functions could be remerged somewhat safely to ensure backwards
	// compatibility with the previous qc library.
	if (qc_is_boundary_set(set)) {
		printf("Trying to free a boundary before the end of execution.\n");
		assert(!qc_is_boundary_set(set));	
	}
	
	if (set->syn_cdll != NULL) {
		cdll_free(set->syn_cdll, NULL);
	}

	if (set->mt_cdll != NULL) {
		// Free the measurements in the set
		n = set->mt_cdll->next;
		while (n != set->mt_cdll) {
			tn = n;
			n = n->next;
			mt = (MEASUREMENT *)tn->key;
			free(tn);

			// If the measurement belongs in multiple sets, then decrement the
			// counter
			if (mt->num_parent_sets > 1) {
				mt->num_parent_sets--;
			}

			// If this is the last set the measurement is in, then free the
			// measurement.
			else {
				qc_free_measurement(mt);
			}
		}
		free(set->mt_cdll);
	}

	if (set->set_pos_cdll != NULL) {
		cdll_free(set->set_pos_cdll, free);
	}

	// If there is a copy of this set, nullify the pointer from the copy to
	// this set.
	if (set->copy != NULL) {
		set->copy->copy = NULL;
	}

	free(set);
}

/**
 * \brief Frees a boundary \ref set
 * 
 * Will free a boundary set. Seperated from qc_free_set as to avoid accidental
 * deletion of boundary sets.
 *
 * \param[in] set The boundary \ref set to be freed.
 */
void qc_free_bdy(SET *set) {
	if (set->syn_cdll != NULL) {
		cdll_free(set->syn_cdll, NULL);
	}
	if (set->ball != NULL) {
		qc_free_ball(set->ball);
	}
	if (set->set_pos_cdll != NULL) {
		cdll_free(set->set_pos_cdll, free);
	}
	free(set->mt_cdll);
	free(set);
}

/**
 * \brief Frees a \ref set. Accepts a void pointer to a \ref set.
 * 
 * \param[in] key The \ref set to be freed.
 *
 * \sa qc_free_set
 */
void qc_free_void_set(void *key) {
	qc_free_set((SET *)key);
}

/**
 * \brief Inserts a \ref set into a \ref qc
 * 
 * \param[in] qc The \ref qc to insert the \ref set into
 * \param[in] set The \ref set to be inserted
 */
void qc_insert_set(QC *qc, SET *set) {
	bh_insert_key(qc->set_heap, set, qc_set_lt_big_t, qc_set_swap,
		qc_set_set_qc_heap_i, NULL);
}

/**
 * \brief Uninserts a \ref set from \ref qc
 * 
 * \param[in] qc The \ref qc to insert the \ref set into
 * \param[in] set The \ref set to be uninserted
 */
void qc_uninsert_set(QC *qc, SET *set) {
	bh_remove_key(qc->set_heap, set->qc_heap_i, qc_set_lt_big_t, qc_set_swap,
		qc_set_set_qc_heap_i, NULL);
}

/**
 * \brief Compares two \link set Sets\endlink by their values of big_t and
 * returns if the first is less than the second.
 * 
 * \param[in] key1 The first set to be compared
 * \param[in] key2 The second set to be compared
 *
 * \return 1 if key1 is less than key2, 0 otherwise.
 */
int qc_set_lt_big_t(void *key1, void *key2) {
	return ((SET *)key1)->big_t < ((SET *)key2)->big_t;
}

/**
 * \brief Swaps the position of two \link set Sets\endlink in a \ref heap.
 * 
 * \param[in] h The \ref heap that both \link set Sets\endlink belong to
 * \param[in] i1 The index of the first \ref set to be swapped in the \ref heap 
 * \param[in] i2 The index of the second \ref set to be swapped in the \ref heap 
 */
void qc_set_swap(BHEAP *h, int i1, int i2) {
	void *key;
	SET *set1, *set2;
	int i;

	key = h->k[i1];
	h->k[i1] = h->k[i2];
	h->k[i2] = key;

	set1 = (SET *)h->k[i1];
	set2 = (SET *)h->k[i2];

	i = set1->qc_heap_i;
	set1->qc_heap_i = set2->qc_heap_i;
	set2->qc_heap_i = i;
}

/**
 * \brief Sets a link from a \ref set to the location of said \ref set to
 * the heap inside a \ref qc
 * 
 * \param[in] key The \ref set to be linked
 * \param[in] i The index location of the \ref set in the \ref heap. 
 */
void qc_set_set_qc_heap_i(void *key, int i) {
	((SET *)key)->qc_heap_i = i;
}

/**
 * \brief Returns where a \ref set is a boundary or not
 * 
 * \param[in] set The \ref set to test
 *
 * \return 1 if the set is a boundary, 0 otherwise.
 */
int qc_is_boundary_set(SET *set) {
	return (set->type == PRIMAL_BOUNDARY || set->type == DUAL_BOUNDARY);
}

/**
 * \brief Merges two dots when neither of them are the boundary
 * 
 * \param[in] m The \ref matching to store the undo steps in
 * \param[in] dot1 The \ref dot being kept after the merge
 * \param[in] dot2 The \ref dot being removed after the merge
 * \param[in] lln_save The \ref ll_node of all the lines coming from dot2
 * before the merge 
 */
void qc_merge_dots_without_boundary(MATCHING *m, DOT *dot1, DOT *dot2, LL_NODE *lln_save) {
	LL_NODE *lln, *kept_lines; 
	DLL_NODE *dlln;
	LINE *line, *line2;
	DOT *d1, *d2, *dot_keep;
	HT *line_ht;
	int hash, found;
	double p_new;

	kept_lines = ll_create();
	line_ht = ht_create(MERGE_LINES_HT_SIZE);

	// Loop over all the lines from the dot we are keeping and remove
	// any internal lines (lines that go from this dot to itself now
	// that the other dot is gone). 

	lln = dot1->lines;
	while (lln != NULL) {
		line = (LINE *)lln->key;	
		d1 = line->a;
		d2 = line->b;

		// If either dot at the end of the lines is a merged dot, it
		// wont have any lines, so recurse down the mergings until we
		// get to the dot containing the relevant lines
		while (d1->merge != NULL) {
			d1 = d1->merge;
		}
		while (d2->merge != NULL) {
			d2 = d2->merge;
		}

		// Since we have just merged two sets, at minimum, the line
		// connecting the two sets will now appear to be a line to
		// itself as traversing down the dot from the removed set will
		// result in the dot that was kept. 
		// 
		// If the line is between itself, it is internal, so we
		// do not want to keep the line. Create an UNDO so we can get
		// the line back later, and then continue.
		if (d1 == d2) {
			if (m->undo_flag == TRUE) {
				m_create_line_undo(m, UNDO_LINE_MOVE, line, dot1, NULL, 0, 0);
			} 
			lln = lln->next;
			continue;
		}

		// We want to hash the line based on its destination, this lets
		// us easily determine which lines we already have later.
		if (d1 == dot1) {
			hash = qc_hash_ijt(d2->i, d2->j, d2->little_t);
		} else {
			hash = qc_hash_ijt(d1->i, d1->j, d1->little_t);
		}
		ht_insert_key(line_ht, hash, line);

		kept_lines = ll_insert(kept_lines, line);
		lln = lln->next;
	}

	// Loop over all the lines in the dot that is being removed.
	// Each line is checked to see if its destination is already in the
	// hash. If so, this means both dot1 and dot2 are connected share a
	// destination. In this instance, the two lines are merged
	// together. If the destination is new, then it is simply inserted.
	
	lln = lln_save; 
	while (lln != NULL) {
		line = (LINE *)lln->key;	
		d1 = line->a;
		d2 = line->b;

		// This line is between the dot (dot2) we are removing and
		// another dot which we are going to keep. This dot can also be
		// thought of as the "destination" of the line.
		dot_keep = (d1 == dot2) ? d2 : d1;

		// Hash the line based on the location of the destination 
		assert(dot_keep->merge == NULL);
		hash = qc_hash_ijt(dot_keep->i, dot_keep->j, dot_keep->little_t);

		// Search the hash table created when looping over the lines in
		// dot1 and look to see if there is already line to the destination 
		found = FALSE;
		dlln = ht_hash_lookup(line_ht, hash); 

		while (dlln != NULL) {
			line2 = (LINE *)dlln->key;
			if (line2->a == line->a || line2->b == line->b) {
				found = TRUE;
				break;
			}
			dlln = dlln->next;
		}

		// dot_keep is already completely merged, but the dX corresponding 
		// to dot2 isn't currently merged. So traverse on down!
		while (d1->merge != NULL) {
			d1 = d1->merge;
		}
		while (d2->merge != NULL) {
			d2 = d2->merge;
		}

		// If the line is internal, then flag it as moved to no where
		// and continue looking at more lines 
		if (d1 == d2) {
			if (m->undo_flag == TRUE) {
				m_create_line_undo(m, UNDO_LINE_MOVE, line, dot2, NULL, 0, 0);
			} else {
				free(line);
			}
			lln = lln->next;
			continue;
		}

		d1 = line->a;
		d2 = line->b;

		// If the line already exists in the hash table, then we want
		// to merge the probabilities of the two lines into one line
		if (found) {
			// One end of the line has to be our dot 
			assert(line->a == dot2 || line->b == dot2);

			// Ensure that no line has been left un updated after a
			// constituent dot has been merged
			assert(line->a->merge == NULL || line->a == dot2);
			assert(line->b->merge == NULL || line->b == dot2);
			assert(line2->a->merge == NULL);
			assert(line2->b->merge == NULL);

			// Create an undo for the line at its current weight
			if (m->undo_flag == TRUE) {
				m_create_line_undo(m, UNDO_LINE_WT, line2, NULL, NULL, line2->wt, 0);
			}

			// Convert the weights back to probabilities, then add them
			p_new = exp(-line->wt/PRECISION) + exp(-line2->wt/PRECISION);

			// Convert the new probability back to a weight
			line2->wt = -log(p_new) * PRECISION;
			if (line2->wt < PRECISION) {
				line2->wt = PRECISION;
			}
			line2->wt -= line2->wt % 2;

			// We are removing this line from both end points. Hence,
			// we need to be able to undo and add it back to both ends.
			// The end in dot_keep needs to be explictly remove, the
			// end in dot2 is removed by nullifying the entire lines list.

			if (m->undo_flag == TRUE) {
				// dot2 has its lines nullified, so we need an undo to
				// put this line back given that it is being removed
				// and isn't being moved into another dot.
				m_create_line_undo(m, UNDO_LINE_MOVE, line, dot2, NULL, 0, 0);

				// Both line and line2 have the same destination at this
				// point. We want to remove line from the dot we are
				// keeping, as we have updated line2 to have the merged
				// weight of the both lines.
				//
				// We don't want to do this if the dot being kept is a boundary
				// as the line was never in there in the first place.
				if (!dot_keep->v || dot_keep->v->v_num >= 0) { 
					m_create_line_undo(m, UNDO_LINE_MOVE, line, dot_keep, NULL, 0, 0);
				}
			}
			dot_keep->lines = ll_delete_node(dot_keep->lines, line);
			if (m->undo_flag != TRUE) {
				assert(d1 == dot2 || d2 == dot2);
				assert(d1 == dot_keep || d2 == dot_keep);
				free(line);
			}

			// If we have merged lines, then we have nothing else to do
			// here as line2 already has the correct end points since
			// it is from dot1 (the dot being merged into). It's also
			// already in the hash table and the kept_lines list.
			lln = lln->next;
			continue;
		}

		// dot_keep is already completely merged, but the dX corresponding 
		// to dot2 isn't currently merged. So traverse on down!
		while (d1->merge != NULL) {
			d1 = d1->merge;
		}
		while (d2->merge != NULL) {
			d2 = d2->merge;
		}

		// One end of the line has to be our dot 
		assert(line->a == dot2 || line->b == dot2);

		// Ensure that no line has been left un updated after a
		// constituent dot has been merged
		assert(line->a->merge == NULL || line->a == dot2);
		assert(line->b->merge == NULL || line->b == dot2);

		// Insert the line into the hash table
		ht_insert_key(line_ht, hash, line);
			
		// Update the line so that the end originally pointing to dot2,
		// now points to dot1
		if (m->undo_flag == TRUE) {
			m_create_line_undo(m, UNDO_LINE_EDIT, line, line->a, line->b, 0, 0);
		}
		if (line->a == dot2) {
			line->a = dot1;
		} else {
			line->b = dot1;
		}

		// We are moving this line from dot2 to dot1, so we need to
		// keep an undo so that we can move the line back later
		if (m->undo_flag == TRUE) {
			m_create_line_undo(m, UNDO_LINE_MOVE, line, dot2, dot1, 0, 0);
		}
		kept_lines = ll_insert(kept_lines, line);

		lln = lln->next;
	}
	
	ht_free(line_ht, NULL);

	ll_free(dot1->lines, NULL);
	ll_free(lln_save, NULL);

	// Update dot1 with the lines we've kept!
	dot1->lines = kept_lines;
}

/**
 * \brief Merges two dots when one of them is the boundary
 * 
 * \param[in] m The \ref matching to store the undo steps of the merge in
 * \param[in] dot1 The \ref dot of the boundary to be merged into
 * \param[in] dot2 The \ref dot being merged into the boundary
 * \param[in] lln_save The \ref ll_node of lines that existed in dot2 before
 * being merged 
 */
int qc_merge_dots_with_boundary(MATCHING *m, DOT *dot1, DOT *dot2, LL_NODE *lln_save) {
	LL_NODE *lln, *lln2, *llnk;
	LINE *line, *line2;
	DOT *d1, *d2, *dot_keep;
	double p_new;
	int found;
	int logical_error;

	//m_print_dot(dot1);
	//m_print_dot(dot2);

	logical_error = FALSE;

	// When the dot being kept (dot1) is a boundary, we do not have the luxury
	// of a list of lines to compare between to determine when to merge. 
	lln = lln_save; 
	while (lln != NULL) {
		line = (LINE *)lln->key;	
		d1 = line->a;
		d2 = line->b;

		// dot1 in this instance is a boundary, so all lines will be moved from
		// dot2 to nowhere, regardless of what is done elsewhere.
		if (m->undo_flag == TRUE) {
			m_create_line_undo(m, UNDO_LINE_MOVE, line, dot2, NULL, 0, 0);
		}

		// This line is between the dot (dot2) we are removing and
		// another dot which we are going to keep. This dot can also be
		// thought of as the "destination" of the line.
		dot_keep = (d1 == dot2) ? d2 : d1;
		assert(dot_keep->merge == NULL);
	
		// If this line is simply between the two dots being merged, do
		// nothing, as we just remove internal edges (done already as
		// dot2->lines was nullified in qc_merge_sets)
		if ((d1 == dot1 && d2 == dot2) || (d2 == dot1 && d1 == dot2)) {
			assert(dot_keep->lines == NULL);
			if (m->undo_flag != TRUE) {
				free(line);
			}
			lln = lln->next;
			continue;
		}

		// We should never be merging lines into a boundary
		assert(dot_keep->i >= 0);

		// Loop the lines on the destination dot, looking to see if this
		// destination is also connected to the same boundary as dot2. If it
		// is, then the connection from dot_keep to the boundary and dot2 to
		// the boundary need to be merged.
		found = FALSE;
		lln2 = dot_keep->lines;
		while (lln2 != NULL) {
			line2 = (LINE *)lln2->key;

			assert(line->a->merge == NULL || line->a == dot2);
			assert(line->b->merge == NULL || line->b == dot2);
			assert(line2->a->merge == NULL || line2->a == dot2);
			assert(line2->b->merge == NULL || line2->b == dot2);

			llnk = lln2->next;


			// Keep searching until we find a line with the boundary set (dot1)
			// at one of the ends 
			if (line2->a != dot1 && line2->b != dot1) {
				lln2 = llnk;
				continue;
			}

			// Create an undo for the line at its current weight
			if (m->undo_flag == TRUE) {
				m_create_line_undo(m, UNDO_LINE_WT, line2, NULL, NULL, line2->wt, 0);
			}

			// Convert the weights back to probabilities, then add them
			p_new = exp(-line->wt/PRECISION) + exp(-line2->wt/PRECISION);

			// Convert the new probability back to a weight
			line2->wt = -log(p_new) * PRECISION;
			if (line2->wt < PRECISION) {
				line2->wt = PRECISION;
			}
			line2->wt -= line2->wt % 2;

			// We have merged the line between dot_keep and dot2 into the line
			// between dot_keep and dot1, hence we now want to remove line from
			// dot_keep, else there would be two lines to the boundary.
			if (m->undo_flag == TRUE) {
				m_create_line_undo(m, UNDO_LINE_MOVE, line, dot_keep, NULL, 0, 0);
			}	
			dot_keep->lines = ll_delete_node(dot_keep->lines, line);

			found = TRUE;
			break;
		}

		// If a line was found and merged, then we don't need to update the
		// line as it already has the correct end points
		if (found) {
			if (m->undo_flag != TRUE) {
				free(line);
			}

			lln = lln->next;
			continue;
		}

		assert(line->a == dot2 || line->b == dot2);
		assert(line->a->merge == NULL || line->a == dot2);
		assert(line->b->merge == NULL || line->b == dot2);

		// Update the end points of the line so that references to dot2 are now
		// pointing to dot1
		if (m->undo_flag == TRUE) {
			m_create_line_undo(m, UNDO_LINE_EDIT, line, line->a, line->b, 0, 0);
		}
		if (line->a == dot2) {
			line->a = dot1;
		} else {
			line->b = dot1;
		}

		// If the line now has a boundary in line->a, then swap the
		// order of the two dots. This is simply because matching
		// determines when it has found all the lines base on the
		// t value of line->a
		if (line->a->i < 0) {
			dot_keep = line->a;
			line->a = line->b;
			line->b = dot_keep;
		}

		lln = lln->next;
	}
	ll_free(lln_save, NULL);

	return (!logical_error);
}

/**
 * \brief Merges two \link set Sets\endlink together
 *
 * Will return the merge result in the first set and free the memory associated
 * with the second set.
 * 
 * \param[in] qc The \ref qc containing both sets
 * \param[in] set1 The first \ref set to be merged 
 * \param[in] set2 The second \ref set to be merged 
 *
 * \return The first \ref set provided, now merged with the second \ref set.
 */
SET *qc_merge_sets(QC *qc, SET *set1, SET *set2) {
	SET *setr, *setk;
	SYNDROME *syn;
	CDLL_NODE *node;
	LL_NODE *lln_save, *lln;
	DOT *dot1, *dot2, *dot_keep;
	SET_POS *sp;
	MATCHING *m;
	MEASUREMENT *mt;
	LINE *line;

	/*
	printf("### qc_merge_sets ###\n");
	qc_print_set(set1);
	qc_print_set(set2);
	//*/

	// Assert that the sets being merged are the same type
	assert(
		((set1->type == PRIMAL || set1->type == PRIMAL_BOUNDARY) && 
			(set2->type == PRIMAL || set2->type == PRIMAL_BOUNDARY)) || 
		((set1->type == DUAL || set1->type == DUAL_BOUNDARY) && 
			(set2->type == DUAL || set2->type == DUAL_BOUNDARY)));

	// We cannot merge sets if one of them is NULL
	assert(set1 != NULL);
	assert(set2 != NULL);

	// qc_meas_qubit should ensure that we never reach this case.	
	assert(set1 != set2);

	// Sets shouldn't have a big_t value yet if they are being merged
	assert(set1->big_t == LONG_MAX);
	assert(set2->big_t == LONG_MAX);

	setk = (qc_is_boundary_set(set1) || set1->t > set2->t) ? set1 : set2;
	setr = (setk == set1) ? set2 : set1;

	/*
	printf("Keep Set: "); qc_print_set(setk);
	printf("Remove Set: "); qc_print_set(setr);
	//*/


	//qc_print_set(set1);
	//qc_print_set(set2);
	//printf("%p %p\n", set1->bdy, set2->bdy);


	// If we are merging a boundary set with another boundary set, or we are
	// merging a boundary set with a set that is adjacent to the other boundary
	// (effectively about to merge two boundaries) then we want to log a
	// logical error and abort the merging process
	if ((qc_is_boundary_set(set1) && qc_is_boundary_set(set2)) || 
		(qc_is_boundary_set(setk) && (setr->bdy != NULL && setr->bdy != setk))) {

		if (setk->type == PRIMAL || setk->type == PRIMAL_BOUNDARY) {
			qc->logical_loss_err_pr += 1;
		} else {
			qc->logical_loss_err_du += 1;
		}

		//printf("\nM Logical Loss Error\n");
		return NULL;
	}
	
	// If we are using an advanced recipe, then we have dots already in
	// existance at the time of merging. We therefore need to merge the dots in
	// order to be able to continue.
	if (qc->recipe_adv) {	

		// If we are merging into a boundary set, there is a chance that the
		// set being lost is connected to both boundaries. If we are merging
		// this set, then this is a logical error. 
		//
		// Example:
		// v1 is connected to bdy1, v2 is connected to bdy2, v1 and v2 are
		// connected. v1 is merged with bdy1, then v2 is merged with bdy2. When
		// v2 is merged with bdy2, it is connected to both bdy1 and bdy2, so
		// this action will merge the boundaries, but is NOT caught by the
		// previous check for logical loss errors.

		if (qc_is_boundary_set(setk)) {
			// Loop the lines of the dot being removed
			lln = setr->ball->dot->lines; 
			while (lln != NULL) {
				line = (LINE *)lln->key;	
				
				// Find the destination
				dot_keep = (line->a == setr->ball->dot) ? line->b : line->a;

				// If the dot is a boundary, but is not the boundary we expect,
				// then we have multiple boundaries, and should increase the
				// loss error count and return immediately.
				if (dot_keep->ball->type == PRIMAL_BOUNDARY || dot_keep->ball->type == DUAL_BOUNDARY) {
					if (dot_keep != setk->ball->dot) {
						if (set1->type == PRIMAL || set1->type == PRIMAL_BOUNDARY) {
							qc->logical_loss_err_pr += 1;
						} else {
							qc->logical_loss_err_du += 1;
						}

						//printf("\nM2 Logical Loss Error\n");
						return NULL;
					}
				}
				lln = lln->next;
			}
		}
		
		// Sets shouldn't be finalized if they are being merged. 
		assert(setr->ball->big_t == LONG_MAX);
		assert(setk->ball->big_t == LONG_MAX);

		// If the removed set is NOT a boundary, we expect it to have no t value.
		if (setr->type == PRIMAL || setr->type == DUAL) {
			assert(setr->ball->dot->t == LONG_MAX);
		}
		// If the kept set is NOT a boundary, we expect it to have no t value.
		if (setk->type == PRIMAL || setk->type == DUAL) {
			assert(setk->ball->dot->t == LONG_MAX);
		}
		
		// Determine the type of the set being kept, and pick the correct matching
		m = (setk->type == PRIMAL || setk->type == PRIMAL_BOUNDARY) ? qc->m_pr : qc->m_du;

		dot1 = setk->ball->dot;
		dot2 = setr->ball->dot;

		// The dot being kept should never already have been merged into
		// another dot.
		assert(dot1->merge == NULL);

		// We are about to merge the dot2 into dot1, so we need to be able to
		// undo this process 
		if (m->undo_flag == TRUE) {
			m_create_line_undo(m, UNDO_DOT_MERGE, NULL, dot2, dot2->merge, 0, 0);
		}

		// Save the list of lines from dot2, before NULLing the list and
		// setting dot2 to be merged into dot1.  
		lln_save = dot2->lines;
		dot2->lines = NULL;
		dot2->merge = dot1;

		if (!qc_is_boundary_set(setk)) {
			qc_merge_dots_without_boundary(m, dot1, dot2, lln_save);
		} else {
			qc_merge_dots_with_boundary(m, dot1, dot2, lln_save);
		}

		if (m->undo_flag == FALSE) {
			// If this process is not going to be undone, this ball will never
			// be used again and hence can be freed.
			if (setr->ball->type == PRIMAL) {
				qc_delete_ball(qc->nest_pr, setr->ball);
			} else {
				assert(setr->ball->type == DUAL);
				qc_delete_ball(qc->nest_du, setr->ball);
			}
		}

		if (m->undo_flag == FALSE) {
			// After this loop ends, there is no longer any lines pointing to
			// this dot. This dot continues to exist as the lattice uses it to
			// find what it has been merged into. That is, this dot still
			// resides where it does, however dot->merge will be followed to
			// find the dot1 that was just merged with. There is, however, no
			// way to determine where these `zombie` dots when it comes to time
			// delete, and hence they need to be tracked in this list. 
			dot2->t = qc->big_t;
			cdll_insert_head(qc->zombie_dots, dot2, NULL);
			//printf("Zombie Dots: %d\n", cdll_length(qc->zombie_dots));
		}
	}
	

	/***********************/

	// Remove set2 from the quantum computer
	qc_uninsert_set(qc, setr);

	// Add the measurements from the set being removed to the set being kept.
	// Given that boundaries do not track remaining measurements, only do this
	// if the set being kept is not a boundary.
	if (!qc_is_boundary_set(setk)) {
		setk->num_meas_left += setr->num_meas_left;

		// If the set being kept wasn't next to the boundary, give it the boundary
		// of the set being removed.
		if (setk->bdy == NULL) {
			setk->bdy = setr->bdy;
		}
	}

	// Merge the set_pos cdlls into the set being kept. This can/is used from
	// ex to update tracking of set locations.
	cdll_join(setk->set_pos_cdll, setr->set_pos_cdll);

	// Move the syndromes into the kept set. We do this even with boundaries as
	// syndromes are recycled and they need to be kept somewhere.
	cdll_join(setk->syn_cdll, setr->syn_cdll);
	
	if (!qc_is_boundary_set(setk)) {
		// Join the measurement cdlls from the set being removed into the set
		// being kept. 
		cdll_join(setk->mt_cdll, setr->mt_cdll);

		// Nullify the pointers for the syndrome and measurement cdlls so that the
		// set can be removed cleanly.
		setr->mt_cdll = NULL;
	} else {
		node = setr->mt_cdll->next;
		assert(setk->ball->mp == 1);
		setk->ball->mp = 1;
		while (node != setr->mt_cdll) {
			mt = (MEASUREMENT *)node->key;
			setk->ball->mp *= mt->m;
			node = node->next;
		}
	}

	// Loop the new set_pos cdll and make sure they all point to the kept set
	node = setk->set_pos_cdll->next;
	while (node != setk->set_pos_cdll) {
		sp = (SET_POS *)node->key;
		sp->set = setk;	
		node = node->next;
	}

	// Update all the syndromes in the newly merged set to point to the set we
	// are keeping.
	node = setk->syn_cdll->next;
	while (node != setk->syn_cdll) {
		syn = (SYNDROME *)node->key;
		if (syn->set != setk) {
			syn->set = setk;
		}

		// If merging into the boundary, advance the syndrome.
		if (setk->type == PRIMAL_BOUNDARY || setk->type == DUAL_BOUNDARY) {
			syn->t = 1;
			bh_bubble_down_key(qc->syn_heap, syn->qc_heap_i,
				qc_syndrome_lt_t, qc_syndrome_swap);
		}

		node = node->next;
	}

	setr->syn_cdll = NULL;
	setr->set_pos_cdll = NULL;

	qc_free_set(setr);

	//qc_print_set(setk);
	//if (setk->copy) {
	//	qc_print_set(setk->copy);
	//}
	return setk;
}

/**
 * \brief Finalizes a \ref set 
 * 
 * \param[in] qc The \ref qc the \ref set is in
 * \param[in] set The \ref set to finalize
 * \param[in] t The time to finalize 
 */
void qc_finalize_set(QC *qc, SET *set, long int t) {
	CDLL_NODE *n, *n2;
	SYNDROME *syn;
	MEASUREMENT *mt;
	ERROR *err;
	DLL_NODE *dlln;
	DE *de;
	STICK *stick;
	int i;

	if (qc_is_boundary_set(set)) {
		n = set->syn_cdll->next;
		while (n != set->syn_cdll) {
			syn = (SYNDROME *)n->key;
			if (syn->set == set) {
				syn->t = t;
				bh_bubble_down_key(qc->syn_heap, syn->qc_heap_i,
					qc_syndrome_lt_t, qc_syndrome_swap);
			}
			n = n->next;
		}

		// Get the syndrome qubit with the lowest t
		syn = (SYNDROME *)qc->syn_heap->k[1];

		// If the lowest syndrome qubit t is greater than 0, then all syndrome
		// qubits have been measured. Increment big_t and finalize nests.
		//
		// With asynchronous circuits, errors can create sticks with both balls
		// having big_t two greater
		if (syn->t > 0) {
			qc_increment_big_t(qc);
			qc_finalize_nests(qc, qc->big_t - 3);
		}

		return;
	}

	// Initialised to NULL to avoid using unitialised (quiet compiler warning) 
	HT *err_ht = NULL;

	// Update the big_t in the set to be the current big_t, the bubble its
	// location in the set heap inside qc
	if (set->t == LONG_MAX) {
		set->t = t;
	}
	if (set->t != t) {
		printf("set->t: %ld, t: %ld, set->big_t: %ld\n", set->t, t, set->big_t);
		assert(set->t == t);
	}
	set->big_t = qc->big_t;
	bh_bubble_up_key(qc->set_heap, set->qc_heap_i,
		qc_set_lt_big_t, qc_set_swap);

	// If we are using an advanced recipe, these values were set at the time of
	// creation of the set
	if (!qc->recipe_adv) {
		// Create a ball for the set and insert it into the correct nest
		set->ball = qc_create_ball_adv(set, STICK_HT_SIZE);

		if (set->type == PRIMAL) {
			qc_insert_ball(qc->nest_pr, set->ball);
		}
		else {
			qc_insert_ball(qc->nest_du, set->ball);
		}
	}

	// Create a hash table for all the errors in the set
	if (qc->track) {
		err_ht = ht_create(set->num_errors * 10);
	}

	// Set the measurement product
	//qc_print_set(set);
	set->ball->mp = 1;

	// Loop the measurements in the set
	n = set->mt_cdll->next;
	while (n != set->mt_cdll) {
		mt = (MEASUREMENT *)n->key;
		
		// Multiply the set's ball's measurement product by the measurement result. 
		set->ball->mp *= mt->m;

		// Loop the errors in the measurement, if they are in the quantum
		// computer, then insert them into the newly created error hash table. 
		if (qc->track) {
			n2 = mt->error_cdll->next;
			while (n2 != mt->error_cdll) {
				err = (ERROR *)n2->key;
				if (err->qc_cdlln != NULL) {
					qc_mod_insert_error_in_ht(err_ht, err);
				}
				n2 = n2->next;
			}
		}
		n = n->next;
	}

	if (qc->track) {
		// Loop the error hash table
		for (i = 0; i < err_ht->length; i++) {
			dlln = err_ht->table[i];
			while (dlln != NULL) {
				err = (ERROR *)dlln->key;
				dlln = dlln->next;
				err->type = set->type;

				de = qc_find_de(qc, err->label, err->type);

				// If there is already a detection event
				if (de != NULL) {
					assert(de->label == err->label);
					assert(set->ball->type == err->type);
					
					// Attempt to find the stick between this set and the
					// detection event's
					stick = qc_find_stick(qc, set->ball, de->set->ball);
					
					// No stick? Create and insert one.
					if (stick == NULL) {
						stick = qc_create_stick(set->ball, de->set->ball);
						if (set->type == PRIMAL) {
							qc_insert_stick(qc->nest_pr, stick);
						}
						else {
							qc_insert_stick(qc->nest_du, stick);
						}
					}

					// Delete the detection event
					qc_delete_de(qc, de);

					// Insert the error into the stick
					stick->error_ll = ll_insert(stick->error_ll, err);
					err->stick = stick;
					stick->num_errors++;
				}

				// No detection event? Create and insert one.
				else {
					de = qc_create_de(set, err);
					qc_insert_de(qc, de);
				}
			}
		}

		ht_free(err_ht, NULL);
	} 
	
	// If tracking is off, and an advanced recipe is available, then we want to
	// mark the ball as finalized by setting its big_t and then inserting it
	// into qc->finalized_balls to be converted later.
	else if (qc->recipe_adv) {
		if (set->ball->type != PRIMAL_BOUNDARY && set->ball->type != DUAL_BOUNDARY) {
			set->ball->big_t = set->big_t;
			assert(set->ball->t == set->t);
		}

		cdll_insert_head(qc->finalized_balls, set->ball, NULL);
	}

	// Update each syndrome qubit, setting their time to be the given t, then
	// update their positions in the quantum computer's syndrome qubit heap
	n = set->syn_cdll->next;
	while (n != set->syn_cdll) {
		syn = (SYNDROME *)n->key;
		if (syn->set == set) {
			syn->t = t;
			bh_bubble_down_key(qc->syn_heap, syn->qc_heap_i,
				qc_syndrome_lt_t, qc_syndrome_swap);
		}
		n = n->next;
	}

	// Get the syndrome qubit with the lowest t
	syn = (SYNDROME *)qc->syn_heap->k[1];

	// If the lowest syndrome qubit t is greater than 0, then all syndrome
	// qubits have been measured. Increment big_t and finalize nests.
	//
	// With asynchronous circuits, errors can create sticks with both balls
	// having big_t two greater
	if (syn->t > 0) {
		qc_increment_big_t(qc);
		qc_finalize_nests(qc, qc->big_t - 3);
	}
}

/** @} */






// MEASUREMENT FUNCTIONS

/** \addtogroup measurement
 * @{
 */

/**
 * \brief Creates a new \ref measurement
 * 
 * \param[in] m The stochastic measurement result (-1, 0, 1)
 *
 * \return The newly created \ref measurement.
 */
MEASUREMENT *qc_create_measurement(int m) {
	MEASUREMENT *mt;

	mt = (MEASUREMENT *)my_malloc(sizeof(MEASUREMENT));

	mt->m = m;
	mt->num_parent_sets = 0;
	mt->num_errors = 0;
	mt->error_cdll = cdll_create();
	mt->copy = NULL;
	
	return mt;
}

/**
 * \brief Frees a \ref measurement
 * 
 * \param[in] mt The \ref measurement to be freed
 */
void qc_free_measurement(MEASUREMENT *mt) {
	CDLL_NODE *n, *nt;
	ERROR *error;

	// Free the error cdll without freeing the errors. Nullify the uplinks from
	// the errors to the cdll being freed.
	n = mt->error_cdll->next;
	while (n != mt->error_cdll) {
		nt = n;
		error = (ERROR *)nt->key;
		error->q_cdlln = NULL;
		n = n->next;
		free(nt);
	}
	free(mt->error_cdll);

	// If there is a copy of this measurement, nullify the pointer from the
	// copy to this measurement.
	if (mt->copy != NULL) {
		mt->copy->copy = NULL;
	}
	free(mt);
}

/**
 * \brief Frees a \ref measurement. Accepts a void pointer to a \ref
 * measurement.
 * 
 * \param[in] mt The \ref measurement to be freed
 */
void qc_free_void_measurement(void *key) {
	qc_free_measurement((MEASUREMENT *)key);
}

/**
 * \brief Inserts a \ref measurement into a \ref set
 * 
 * \param[in] set The \ref set to insert the \ref measurement into
 * \param[in] mt The \ref measurement to be inserted
 */
void qc_insert_measurement(SET *set, MEASUREMENT *mt) {
	cdll_insert_head(set->mt_cdll, mt, NULL);

	// Increment the number of parents this measurement has
	mt->num_parent_sets++;

	// Add any measurement errors to the set error counter
	set->num_errors += mt->num_errors;
}

/** @} */






// NEST FUNCTIONS

/**
 * \addtogroup nest
 * @{
 */

/**
 * \brief Creates a new \ref nest
 * 
 * \param[in] stick_ht_stick The size of the \ref stick \ref ht 
 *
 * \return The newly created \ref nest
 */
NEST *qc_create_nest(int stick_ht_size) {
	NEST *nest;

	nest = (NEST *)my_malloc(sizeof(NEST));

	nest->ball_cdll = cdll_create();
	nest->last_blocked_ball_cdlln = nest->ball_cdll;
	nest->last_converted_ball_cdlln = nest->ball_cdll;
	nest->stick_cdll = cdll_create();
	nest->stick_ht = ht_create(stick_ht_size);
	nest->copy = NULL;
	
	return nest;
}

/**
 * \brief Frees a \ref nest
 * 
 * \param[in] nest The \ref nest to be freed
 */
void qc_free_nest(NEST *nest) {
	/*
	int clen;
	CDLL_NODE *n;

	clen = 0;
	n = nest->ball_cdll->next;
	while (n != nest->ball_cdll) {
		clen++;
		n = n->next;
	}
	printf("len: %d\n", clen);
	*/

	cdll_free(nest->ball_cdll, qc_free_void_ball);
	cdll_free(nest->stick_cdll, qc_free_void_stick);
	ht_free(nest->stick_ht, NULL);

	// If there is a copy of this nest, nullify the pointer from the copy to
	// this nest.
	if (nest->copy != NULL) {
		nest->copy->copy = NULL;
	}
	free(nest);
}

/**
 * \brief Finalizes a \ref nest locally
 * 
 * \param[in] ball The \ref ball to finalize locally around
 */
void qc_finalize_local_nest(BALL *ball) {
	HT *ht;
	int i;
	DLL_NODE *dlln;
	STICK *stick;

	// Loop the sticks in the ball and finalize them
	ht = ball->stick_ht;
	for (i = 0; i < ht->length; i++) {
		dlln = ht->table[i];
		while (dlln != NULL) {
			stick = (STICK *)dlln->key;
			qc_finalize_stick(stick);
			dlln = dlln->next;
		}
	}
}

/**
 * \brief Finalizes the \link nest Nests in a \ref qc 
 * 
 * \param[in] qc The \ref qc with the \link nest Nests\endlink
 * \param[in] big_t The big_t of the finalization 
 */
void qc_finalize_nests(QC *qc, long int big_t) {
	CDLL_NODE *n;
	DE *de;
	SET *set;
	ERROR *err;

	// Moves all errors that are less than or equal to big_t from the
	// unfinalised list into the nests.
	n = qc->error_cdll->prev;
	err = (ERROR *)n->key;
	while (n != qc->error_cdll && err->big_t <= big_t) {
		qc_move_error_to_nest(qc, err);
		n = qc->error_cdll->prev;
		err = (ERROR *)n->key;
	}

	// Finalizes the detection events that are less than or equal to big_t
	n = qc->de_cdll->prev;
	de = (DE *)n->key;
	while (n != qc->de_cdll && de->set->big_t <= big_t) {
		n = n->prev;
		qc_finalize_de(qc, de);
		qc_delete_de(qc, de);
		de = (DE *)n->key;
	}

	// Extract each set from the heap, stopping when big_t in the set is
	// greater than big_t. 
	set = (SET *)bh_get_top(qc->set_heap);
	while (set != NULL && set->big_t <= big_t) {
		// Finalise the around ball of the set
		if (qc->track) {
			qc_finalize_local_nest(set->ball);
		}

		// Remove and free the set
		qc_uninsert_set(qc, set);
		qc_free_set(set);
		set = (SET *)bh_get_top(qc->set_heap);
	}

	// Update where we have finalized up to
	qc->unfinalized_big_t = big_t + 1;
}

/**
 * \brief Trim the \link nest Nests\endlink 
 * 
 * \param[in] qc The \ref qc containing the \link nest Nests\endlink
 * \param[in] big_t The big_t to trim the nests up to 
 */
void qc_trim_nests(QC *qc, long int big_t) {
	/**
	 * \remark This can be called with big_t = qc->unfinalized_big_t - 2
	 */
	CDLL_NODE *s, *n, *nt;
	BALL *b;
	DOT *d;
	ERROR *error;
	long int t;

	// Delete all balls from the primary nest where b->big_t is less than the
	// provided big_t
	s = qc->nest_pr->ball_cdll;
	n = s->prev;
	b = (BALL *)n->key;
	while (n != s && b->big_t <= big_t) {
		n = n->prev;
		qc_delete_ball(qc->nest_pr, b);
		b = (BALL *)n->key;
	}

	// Delete all balls from the dual nest where b->big_t is less than the
	// provided big_t
	s = qc->nest_du->ball_cdll;
	n = s->prev;
	b = (BALL *)n->key;
	while (n != s && b->big_t <= big_t) {
		n = n->prev;
		qc_delete_ball(qc->nest_du, b);
		b = (BALL *)n->key;
	}

	t = qc->m_pr->t - qc->m_pr->t_delete;
	s = qc->zombie_dots;
	n = s->prev;
	d = (DOT *)n->key;
	while (n != s && d->t <= t) {
		n = n->prev;

		ht_delete_key(qc->m_pr->dot_ht, qc_hash_ijt(d->i, d->j, d->little_t), d, NULL);
		ht_delete_key(qc->m_du->dot_ht, qc_hash_ijt(d->i, d->j, d->little_t), d, NULL);
		cdll_delete_node(n->next, m_free_void_dot);

		d = (DOT *)n->key;
	}

	// With asynchronous circuits, an error can lead to a stick with both balls
	// having big_t two greater
	
	// Loop all errors in the nests
	s = qc->nest_error_cdll;
	n = s->prev;
	error = (ERROR *)n->key;
	while (n != s && error->big_t < big_t-1) {
		nt = n;
		n = n->prev;

		assert(error->qc_cdlln == NULL);
		error->next->prev = error->prev;
		error->prev->next = error->next;

		// Delete the error from the nest cdll
		cdll_delete_node(nt, NULL);
		qc->num_nest_errors--;

		// If the error is linked to a qubit, delete the error from the qubit
		if (error->q != NULL) {
			assert(error->q->num_errors == error->q->error_ht->num_elem);
			assert(error->q_cdlln != NULL);
			cdll_delete_node(error->q_cdlln, NULL);
			ht_delete_key(error->q->error_ht, error->label, error, NULL);
			error->q_cdlln = NULL;

			error->q->num_errors--;
			assert(error->q->num_errors == error->q->error_ht->num_elem);
		}

		// If the error isn't linked to a qubit, but it has a link to a
		// q_cdlln, then delete the error from the measurement. The q_cdlln is
		// moved to the measurement and reused, technically the name becomes
		// invalid as it is no longer related to a qubit, but this is done to
		// avoid needing yet another variable.
		else if (error->q_cdlln != NULL) {
			cdll_delete_node(error->q_cdlln, NULL);
			error->q_cdlln = NULL;
		}

		qc_free_error(error);

		error = (ERROR *)n->key;
	}
}

/**
 * \brief Converts \link nest Nests\endlink  
 * 
 * \param[in] qc The \ref qc to convert
 * \param[in] undo Whether or not undo is enabled 
 */
void qc_convert_nests(QC *qc, int undo) {
	// If there is a recipe, convert the nests based on the recipe type
	if (qc->recipe != NULL) {
		if (qc->recipe->type == RECIPE_INFINITE) {
			qc_convert_nest_infinite(qc, undo);
		} else {
			printf("Finite case not yet supported");
			exit(0);
		}
	}
	else if (qc->recipe_adv != NULL) {
		CDLL_NODE *cdlln, *cdlln2, *cdlln_temp, *n, *valid_cdll;
		BALL *ball;
		MATCHING *matching;
		LL_NODE *lln;
		LINE *line, *line2;
		DOT *a, *b, *other;
		int valid;

		/*
		 * This code became highly non-trivial due to a rare edge case where it
		 * was possible that a ball that should be finalized however all the
		 * balls around it belonged to sets that would be finalized in the
		 * future. This would cause a dot to be inserted into the matching with
		 * no viable lines, and hence could not be matched. This code was
		 * updated to ensure that all balls actually had a viable line to
		 * follow before being created.
		 */

		/*
		 * Step 1: Loop all the finalized balls and set them to have a dot->t
		 * value. This essentially marks them as "possibly finalized".
		 */

		//printf("%p %p\n", qc->finalized_balls->key, qc->finalized_balls->prev->key);
		
		//cdll_print(qc->finalized_balls, qc_print_void_ball);

		valid_cdll = cdll_create();

		cdlln = qc->finalized_balls->prev;
		while (cdlln != qc->finalized_balls) {
			ball = (BALL *)cdlln->key;

			if (ball->big_t >= qc->unfinalized_big_t) {
				cdlln = cdlln->prev;
				break;
			}

			matching = (ball->type == PRIMAL) ? qc->m_pr : qc->m_du;

			if (ball->type != PRIMAL_BOUNDARY && ball->type != DUAL_BOUNDARY) {
				ball->dot->t = ball->big_t;
	
				if (ball->dot->little_t != ball->t) {
					m_print_dot(ball->dot);
					assert(ball->dot->little_t == ball->t);
				}

				// Keep track of the potentially valid balls
				cdll_insert_head(valid_cdll, ball, NULL);
			}

			cdlln = cdlln->prev;
		}

		/*
		 * Step 2: Identify which dots have no valid lines because they only
		 * connect to dots that aren't finalized. 
		 */
			
		cdlln = valid_cdll->prev;
		while (cdlln != valid_cdll) {
			ball = (BALL *)cdlln->key;
			matching = (ball->type == PRIMAL) ? qc->m_pr : qc->m_du;

			valid = FALSE;
			lln = ball->dot->lines;
			while (lln != NULL) {
				line2 = (LINE *)lln->key;
				other = (line2->a == ball->dot) ? line2->b : line2->a;
				if (other->t != LONG_MAX) {
					valid = TRUE;
					break;
				}
				lln = lln->next;
			}

			if (!valid) {
				//printf("Can't convert invalid ball.\n");
				//m_print_dot(ball->dot);
				//ll_print(ball->dot->lines, m_print_void_line);

				ball->dot->t = LONG_MAX;

				cdlln_temp = cdlln->prev;
				cdll_delete_node(cdlln, NULL);
				cdlln = cdlln_temp;
				continue;
			}

			cdlln = cdlln->prev;
		}

		cdlln = valid_cdll->prev;
		while (cdlln != valid_cdll) {
			ball = (BALL *)cdlln->key;
			ball->dot->t = LONG_MAX;
			cdlln = cdlln->prev;
		}

		cdlln = valid_cdll->prev;
		while (cdlln != valid_cdll) {
			ball = (BALL *)cdlln->key;
			matching = (ball->type == PRIMAL) ? qc->m_pr : qc->m_du;

			// Insert the dot into the matching
			n = cdll_create_node(ball->dot);
			cdll_insert_node_head(matching->dots, n);
			matching->num_dots++;

			if (matching->undo_flag == TRUE) {
				m_create_undo(matching, UNDO_INSERT_DOT, NULL, ball->dot->t, NULL, NULL, n);

				// We set the undo for the dot time here, as we don't want to
				// create it until after we've confirmed the dot is indeed valid
				m_create_line_undo(matching, UNDO_DOT_T, NULL, ball->dot, NULL, 0, ball->dot->t);
			}

			ball->dot->t = ball->big_t;

			/*
			 * Because all dots are made "valid" ahead of time, both dots will
			 * want to create the edge between the two dots. This used to be
			 * avoided by handling the dots in order, however we can no longer
			 * do that due to needing to work out all the potentially final
			 * dots before we can determine if a dot truly has no lines to
			 * follow. 
			 */

			lln = ball->dot->lines;
			while (lln != NULL) {

				line = (LINE *)lln->key;
				a = line->a;
				b = line->b;

				if (a == ball->dot && b->t > ball->big_t) {
					//printf("SKIP LINE TO FUTURE ");	m_print_line(line);
					lln = lln->next;
					continue;
				}
				if (b == ball->dot && a->t > ball->big_t) {
					//printf("SKIP LINE TO FUTURE ");	m_print_line(line);
					lln = lln->next;
					continue;
				}

				m_insert_line(matching, line);
				lln = lln->next;
			}

			if (ball->mp < 0) {
				VERTEX *v = m_create_vertex(ball->i, ball->j, ball->big_t, ball->dot);
				if (ball->type != PRIMAL_BOUNDARY && ball->type != DUAL_BOUNDARY) {
					m_insert_vertex(matching, v);	
					if (matching->undo_flag == TRUE) {
						m_create_undo(matching, UNDO_CREATE_VERTEX, v, 0, NULL, NULL, NULL);
					}
				}
				else {
					v->v_num = ball->i;

					// Prevents time abort as boundary sets
					// create boundary dots with t = INT_MAX
					ball->dot->t = -1;
				}
			}

			// Clean up this ball from the finalized_balls list
			cdlln2 = qc->finalized_balls->prev;
			while (cdlln2 != qc->finalized_balls) {
				cdlln_temp = cdlln2->prev;
				if (cdlln->key == cdlln2->key) {
					cdll_delete_node(cdlln2, NULL);
					break;
				}
				cdlln2 = cdlln_temp;
			}

			cdlln = cdlln->prev;
		}

		cdll_free(valid_cdll, NULL);

		//printf("%p %p\n", qc->finalized_balls->key, qc->finalized_balls->prev->key);
	}
	else {
		// Convert the nests without recipe guidance
		qc_convert_nest(qc, qc->nest_pr, qc->m_pr, undo);
		qc_convert_nest(qc, qc->nest_du, qc->m_du, undo);
	}
}

/**
 * \brief Converts a \ref nest
 * 
 * \param[in] qc The \ref qc containing the \ref nest
 * \param[in] nest The \ref nest to be converted
 * \param[in] m The \ref matching for the \ref nest 
 * \param[in] undo Whether or not undo is enabled 
 */
void qc_convert_nest(QC *qc, NEST *nest, MATCHING *m, int undo) {
	CDLL_NODE *n;
	BALL *ball, *ball2;
	DOT *dot, *dot2;
	HT *ht;
	int i;
	DLL_NODE *n2;
	STICK *s;

	m->undo_flag = undo;

	// Loop all balls in the nest since the last converted ball until the last
	// of the finalized balls
	n = nest->last_converted_ball_cdlln->prev;
	ball = (BALL *)n->key;
	while (n != nest->ball_cdll && ball->big_t < qc->unfinalized_big_t) {
		// If the ball has no dot, create and insert one
		if (ball->dot == NULL) {
			// If the measurement product is negative, then create a vertex
			if (ball->mp < 0) {
				dot = m_create_and_insert_dot_and_vertex(m, ball);
			}
			else {
				dot = m_create_and_insert_dot(m, ball);
			}
			ball->dot = dot;
		}			

		n = n->prev;
		ball = (BALL *)n->key;
	}

	n = nest->last_converted_ball_cdlln->prev;
	ball = (BALL *)n->key;
	while (n != nest->ball_cdll && ball->big_t < qc->unfinalized_big_t) {
		dot = ball->dot;

		// For each stick in the ball
		ht = ball->stick_ht;
		for (i = 0; i < ht->length; i++) {
			n2 = ht->table[i];
			while (n2 != NULL) {
				s = (STICK *)n2->key;

				// Find the destination ball for the stick
				ball2 = (s->a == ball) ? s->b : s->a;
				
				// If there is no dot on the destination ball, create and
				// insert it.
				if (ball2->dot == NULL) {
					// If the destination doesn't already have a dot, then it's
					// unfinalized, so we don't want to create it yet. However,
					// boundaries wont have dots, so we do want to create a dot
					// for them that this ball is connected to.
					if (ball2->type != PRIMAL_BOUNDARY && ball2->type != DUAL_BOUNDARY) {
						n2 = n2->next;
						continue;
					}

					// If the measurement product is negative, then create a vertex
					if (ball2->mp < 0) {
						dot2 = m_create_and_insert_dot_and_vertex(m, ball2);
					}
					else {
						dot2 = m_create_and_insert_dot(m, ball2);
					}
					ball2->dot = dot2;
				}
				else {
					dot2 = ball2->dot;
				}

				// Create a line between the two dots with the weight of the
				// stick
				m_create_and_insert_line(m, dot, dot2, s->p_stick);
				n2 = n2->next;
			}
		}

		n = n->prev;
		ball = (BALL *)n->key;
	}
	
	nest->last_converted_ball_cdlln = n->next;
}

/** @} */






// BALL FUNCTIONS

/**
 * \addtogroup ball
 * @{
 */

/**
 * \brief Creates a new \ref ball
 * 
 * \param[in] set The \ref set to create the \ref ball from
 *
 * \return The newly created \ref ball
 */
BALL *qc_create_ball(SET *set, int stick_ht_size) {
	return qc_create_ball_raw(set->type, set->i, set->j, set->big_t, stick_ht_size);
}

/**
 * \brief Creates a new \ref ball, suitable for advanced recipe usage
 * 
 * \param[in] set The \ref set to convert to a ball
 * \param[in] stick_ht_size The size of the hash table to use for the \ref stick%s. 
 */
BALL *qc_create_ball_adv(SET *set, int stick_ht_size) {
	return qc_create_ball_raw_adv(set->type, set->i, set->j, set->t, set->big_t, stick_ht_size);
}

/**
 * \brief Creates a new \ref ball from raw parameters
 * 
 * \param[in] type The type of \ref ball
 * \param[in] i The i coordinate
 * \param[in] j The j coordinate
 * \param[in] big_t The big_t coordinate
 * \param[in] stick_ht_size The size of the hash table to use for the \ref stick%s. 
 *
 * \return The newly created \ref ball
 */
BALL *qc_create_ball_raw(int type, int i, int j, long int big_t, int stick_ht_size) {
	BALL *ball;

	ball = (BALL *)my_malloc(sizeof(BALL));

	ball->type = type;
	ball->i = i;
	ball->j = j;
	ball->t = 0;
	ball->big_t = big_t;
	ball->stick_ht = ht_create(stick_ht_size);
	ball->nest_cdlln = NULL;
	ball->mp = 0;
	ball->dot = NULL;
	ball->hash_ijt = INT_MAX;

	ball->copy = NULL;

	return ball;
}

/**
 * \brief Creates a new \ref ball from raw parameters, suitable for an advanced recipe
 * 
 * \param[in] type The type of \ref ball
 * \param[in] i The i coordinate
 * \param[in] j The j coordinate
 * \param[in] t The t coordinate 
 * \param[in] big_t The big_t coordinate
 * \param[in] stick_ht_size The size of the hash table to use for the \ref stick%s. 
 *
 * \return The newly created \ref ball
 */
BALL *qc_create_ball_raw_adv(int type, int i, int j, long int t, long int big_t, int stick_ht_size) {
	BALL *ball;

	ball = (BALL *)my_malloc(sizeof(BALL));

	ball->type = type;
	ball->i = i;
	ball->j = j;
	ball->t = t;
	ball->big_t = big_t;
	ball->stick_ht = ht_create(stick_ht_size);
	ball->nest_cdlln = NULL;
	ball->mp = 0;
	ball->dot = NULL;

	ball->hash_ijt = qc_hash_ball(ball);
	
	ball->copy = NULL;

	return ball;
}

/**
 * \brief Frees a \ref ball
 * 
 * If the ball is a boundary, the dot will also be freed.
 *
 * \param[in] ball The \ref ball to be freed
 */
void qc_free_ball(BALL *ball) {
	// If the ball is a boundary, then also free the dot
	if (ball->type == PRIMAL_BOUNDARY || ball->type == DUAL_BOUNDARY) {
		if (ball->dot != NULL) {
			m_free_dot(ball->dot);
		}
	}
	ht_free(ball->stick_ht, NULL);
	
	// If there is a copy of this ball, nullify the pointer from the copy to
	// this ball
	if (ball->copy != NULL) {
		ball->copy->copy = NULL;
	}

	// If there is a dot and the dot is pointing to this ball, then nullify the
	// pointer. Second condition is because this ball may be a copy
	if (ball->dot != NULL && ball->dot->ball == ball) { 
		ball->dot->ball = NULL;
	}

	free(ball);
}

/**
 * \brief Frees a \ref ball. Accepts a void pointer to a \ref ball.
 * 
 * If the ball is a boundary, the dot will also be freed.
 *
 * \param[in] key The \ref ball to be freed
 */
void qc_free_void_ball(void *key) {
	//printf("FREE VOID BALL\n");
	//qc_print_ball((BALL *)key);
	qc_free_ball((BALL *)key);
}

/**
 * \brief Inserts a \ref ball into a \ref nest
 * 
 * \param[in] nest The \ref nest to insert the \ref ball into
 * \param[in] ball The \ref ball to insert 
 */
void qc_insert_ball(NEST *nest, BALL *ball) {
	cdll_insert_head(nest->ball_cdll, ball, qc_ball_set_nest_cdlln);
}

/**
 * \brief Deletes and frees a \ref ball from a \ref nest
 * 
 * \param[in] nest The \ref nest that the \ref ball is in
 * \param[in] ball The \ref ball to be deleted 
 */
void qc_delete_ball(NEST *nest, BALL *ball) {
	HT *ht;
	int i;
	DLL_NODE *n;
	STICK *s;

	// If there is no ball, do nothing
	if (ball == NULL) {
		return;
	}

	// Delete all sticks from the ball
	ht = ball->stick_ht;
	for (i = 0; i < ht->length; i++) {
		n = ht->table[i];
		while (n != NULL) {
			s = (STICK *)n->key;
			n = n->next;
			qc_delete_stick(nest, s);
		}
	}
	ht_free(ht, NULL);
	
	// Remove from the nest cdll
	cdll_delete_node(ball->nest_cdlln, NULL);

	// If there is a dot, nullify the pointer from the dot to this ball
	if (ball->dot != NULL) {
		ball->dot->ball = NULL;
	}

	free(ball);
}

/**
 * \brief Sets an link from a \ref ball to the \ref cdlln of \link ball
 * Balls\endlink in the \ref nest. 
 * 
 * \param[in] key The \ref ball to be linked
 * \param[in] n The ::CDLL_NODE of the ball in the \ref cdlln in the \ref nest. 
 */
void qc_ball_set_nest_cdlln(void *key, CDLL_NODE *n) {
	((BALL *)key)->nest_cdlln = n;
}

/** @} */




// STICK FUNCTIONS

/**
 * \addtogroup stick
 * @{
 */

/**
 * \brief Creates a new \ref stick between two \link ball Balls\endlink
 * 
 * \param[in] a The first \ref ball to make a stick between
 * \param[in] b The second \ref ball to make a stick between
 *
 * \return The newly created \ref stick
 */
STICK *qc_create_stick(BALL *a, BALL *b) {
	STICK *stick;

	stick = (STICK *)my_malloc(sizeof(STICK));

	stick->p_stick = 0;
	stick->a = a;
	stick->b = b;
	stick->num_errors = 0;
	stick->error_ll = NULL;
	stick->nest_cdlln = NULL;
	stick->copy = NULL;

	return stick;
}

/**
 * \brief Frees a \ref stick
 * 
 * \param[in] stick The \ref stick to be freed
 *
 * \sa qc_free_void_stick
 */
void qc_free_stick(STICK *stick) {
	ll_free(stick->error_ll, NULL);
	
	// If there is a copy of this stick, nullify the pointer from the copy to
	// this stick.
	if (stick->copy != NULL) {
		stick->copy->copy = NULL;
	}
	
	free(stick);
}

/**
 * \brief Frees a \ref stick. Accepts a void pointer to a \ref stick.
 * 
 * \param[in] key The \ref stick to be freed
 *
 * \sa qc_free_stick
 */
void qc_free_void_stick(void *key) {
	qc_free_stick((STICK *)key);
}

/**
 * \brief Finds a \ref stick that connects two \link ball Balls\endlink
 * 
 * \param[in] qc The \ref qc that contains the \link ball Balls\endlink
 * \param[in] b1 The first \ref ball that the \ref stick connects 
 * \param[in] b2 The second \ref ball that the \ref stick connects 
 *
 * \return The \ref stick if it is found, NULL otherwise.
 */
STICK *qc_find_stick(QC *qc, BALL *b1, BALL *b2) {
	int type;
	int ip;
	HT *ht;
	DLL_NODE *n;
	STICK *s;

	// Detemine which nest to use, and get the stick hash table
	type = b1->type;
	if (type == PRIMAL) {
		ht = qc->nest_pr->stick_ht;
	}
	else {
		ht = qc->nest_du->stick_ht;
	}

	// Search the hash table for the stick
	assert(ht->length > 0);
	ip = ((size_t)b1 + (size_t)b2)%INT_MAX;
	n = ht->table[ip%ht->length];
	while (n != NULL) {
		s = (STICK *)n->key;

		// If the stick starts and ends with the balls provided, return the stick
		if ((s->a == b1 && s->b == b2) || (s->a == b2 && s->b == b1)) {
			return s;
		}
		n = n->next;
	}

	return NULL;
}

/**
 * \brief Sets a link from a \ref stick to the \ref cdll in a \ref nest
 * 
 * \param[in] key The \ref stick to link
 * \param[in] n The ::CDLL_NODE of the \ref stick in the \ref nest 
 */
void qc_stick_set_nest_cdlln(void *key, CDLL_NODE *n) {
	((STICK *)key)->nest_cdlln = n;
}

/**
 * \brief Inserts a \ref stick into a \ref nest
 * 
 * \param[in] nest The \ref nest to insert the \ref stick into
 * \param[in] stick The \ref stick to insert 
 */
void qc_insert_stick(NEST *nest, STICK *stick) {
	int ip;
	BALL *a, *b;

	a = stick->a;
	b = stick->b;
	ip = ((size_t)a + (size_t)b)%INT_MAX;

	// Insert the stick into the balls at either end
	if (a->stick_ht->length > 0) {
		ht_insert_key(a->stick_ht, ip, stick);
	}
	if (b->stick_ht->length > 0) {
		ht_insert_key(b->stick_ht, ip, stick);
	}

	// Insert the stick into the nest cdll and hash table
	cdll_insert_head(nest->stick_cdll, stick, qc_stick_set_nest_cdlln);
	ht_insert_key(nest->stick_ht, ip, stick);
}

/**
 * \brief Deletes and frees a \ref stick from a \nest
 * 
 * \param[in] nest The \ref nest to delete the \ref stick from
 * \param[in] stick The \ref stick to delete 
 */
void qc_delete_stick(NEST *nest, STICK *stick) {
	int ip;
	BALL *a, *b;

	a = stick->a;
	b = stick->b;
	ip = ((size_t)a + (size_t)b)%INT_MAX;
	
	// Delete the stick from the balls at either end
	if (a->stick_ht->length > 0) {
		ht_delete_key(a->stick_ht, ip, stick, NULL);
	}
	if (b->stick_ht->length > 0) {
		ht_delete_key(b->stick_ht, ip, stick, NULL);
	}

	// Delete the stick from the cdll and hash table, free the stick
	ht_delete_key(nest->stick_ht, ip, stick, NULL);
	cdll_delete_node(stick->nest_cdlln, qc_free_void_stick);
}

/**
 * \brief Compares two \link ball Balls\endlink for equality
 * 
 * \param[in] key1 The first \ref ball to compare
 * \param[in] key2 The second \ref ball to compare
 */
int qc_stick_eq_balls(void *key1, void *key2) {
	STICK *s1, *s2;

	s1 = (STICK *)key1;
	s2 = (STICK *)key2;

	return ((s1->a == s2->a && s1->b == s2->b) || (s1->a == s2->b && s1->b == s2->a));
}

/**
 * \brief Finalizes a \ref stick
 * 
 * \param[in] stick The \ref stick to be finalized
 */
void qc_finalize_stick(STICK *stick) {
	LL_NODE *n;
	double *p_arr, sum;
	int i;
	ERROR *err, **err_arr;

	// If the stick has a negative probability, error
	if (stick->p_stick < 0) {
		printf("p_stick = %g < 0\n", stick->p_stick);
		assert(0);
	}

	// If the stick already has a p_stick, do nothing
	if (stick->p_stick != 0) {
		return;
	}

	p_arr = (double *)my_calloc(stick->num_errors, sizeof(double));
	err_arr = (ERROR **)my_calloc(stick->num_errors, sizeof(ERROR *));

	// Create an array of the errors in the stick
	n = stick->error_ll;
	i = 0;
	while (n != NULL) {
		err_arr[i] = (ERROR *)n->key;
		n = n->next;
		i++;
	}
	assert(i == stick->num_errors);

	// Simple calculation of stick probability
	// Sum the probabilities of each of the errors in the stick
	//
	// This version is simpler, but also the most cost efficient
	sum = 0;
	n = stick->error_ll;
	while (n != NULL) {
		err = (ERROR *)n->key;
		sum += err->p;
		n = n->next;
	}

	stick->p_stick = sum;

	free(p_arr);
	free(err_arr);
}

/** @} **/







// RECIPE FUNCTIONS

/**
 * \addtogroup recipe
 * @{
 */

/**
 * \brief Creates a \ref recipe
 * 
 * \param[in] type The type of \ref recipe
 * \param[in] size The size of \ref recipe (number of layers)
 * \param[in] n The size of the n dimension of the computation 
 * \param[in] m The size of the m dimension of the computation 
 *
 * \return The newly created recipe
 */
RECIPE *qc_create_recipe(int type, int size, int n, int m) {
	RECIPE *recipe;

	recipe = (RECIPE *)my_malloc(sizeof(RECIPE));
	recipe->type = type;
	recipe->size = size;

	recipe->n = n;
	recipe->m = m;
	
	recipe->min_horz_wt_pr = INT_MAX;
	recipe->min_horz_wt_du = INT_MAX;

	recipe->num_layers = 0;
	recipe->layers = (LAYER **)my_malloc(size * sizeof(LAYER *));
	recipe->blocks = NULL;
	recipe->offsets = NULL;

	recipe->dotarr1 = (DOT ***)my_2d_calloc(n, m, sizeof(DOT *));
	recipe->dotarr2 = (DOT ***)my_2d_calloc(n, m, sizeof(DOT *));
	
	recipe->repeat_id = -1;
	recipe->repeated_count = 0;

	return recipe;
}

/**
 * \brief Creates a \ref recipe_adv
 *
 * \param[in] type The type of \ref recipe_adv 
 * 
 * \return The newly created \ref recipe_adv
 */
RECIPE_ADV *qc_create_recipe_adv(int type, int n, int m, int copy) {
	RECIPE_ADV *recipe;
	int i, j;

	recipe = (RECIPE_ADV *)my_malloc(sizeof(RECIPE_ADV));
	recipe->n = n;
	recipe->m = m;
	recipe->type = type;
	recipe->in_cycle = FALSE;
	recipe->cycle_t0 = LONG_MAX;
	recipe->cycle_period = 0;
	recipe->cycle_len = 0;
	recipe->cycle_begin = NULL;
	recipe->cycle_end = NULL;

	recipe->min_horz_wt_pr = INT_MAX;
	recipe->min_horz_wt_du = INT_MAX;

	recipe->t_arr = (long int **)my_2d_calloc(n, m, sizeof(long int));
	recipe->block_arr = (CDLL_NODE ***)my_2d_calloc(n, m, sizeof(CDLL_NODE *));

	if (!copy) {
		recipe->block_cdll = cdll_create();
		recipe->block_ht = ht_create(RECIPE_ADV_BLOCK_HT_SIZE);
		recipe->offset_ht = ht_create(RECIPE_ADV_OFFSET_HT_SIZE);
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				recipe->block_arr[i][j] = cdll_create();
			}
		}
	}

	return recipe;
}

/**
 * \brief Frees a \ref recipe
 * 
 * \param[in] recipe The recipe to be freed
 */
void qc_free_recipe(RECIPE *recipe) {
	int i;

	// If there is no recipe, do nothing
	if (recipe == NULL) {
		return;
	}
	// Free each of the layers
	for (i = 0; i < recipe->size; i++) {
		qc_free_layer(recipe->layers[i], recipe->n);
	}
	free(recipe->layers);

	// Free the dot arrays
	my_2d_free(recipe->n, (void **)recipe->dotarr1);
	my_2d_free(recipe->n, (void **)recipe->dotarr2);

	// Free the blocks and offsets
	ll_free(recipe->blocks, qc_free_void_block);
	ll_free(recipe->offsets, free);
	free(recipe);
}

/**
 * \brief Frees an advanced recipe
 * 
 * \param[in] recipe The \ref recipe_adv to free
 */
void qc_free_recipe_adv(RECIPE_ADV *recipe) {
	int i, j;
	// If there is no recipe, do nothing
	if (recipe == NULL) {
		return;
	}

	for (i = 0; i < recipe->n; i++) {
		for (j = 0; j < recipe->m; j++) {
			cdll_free(recipe->block_arr[i][j], NULL);
		}
	}

	my_2d_free(recipe->n, (void **)recipe->block_arr);
	my_2d_free(recipe->n, (void **)recipe->t_arr);

	// Free the blocks and offsets
	ht_free(recipe->block_ht, qc_free_void_block);
	ht_free(recipe->offset_ht, free);
	cdll_free(recipe->block_cdll, NULL);
	free(recipe);
}

/**
 * \brief Resets an advanced recipe
 * 
 * \param[in] recipe The \ref recipe_adv to be reset
 */
void qc_reset_recipe_adv(RECIPE_ADV *recipe) {
	int i;
	int j;
	CDLL_NODE *cdlln;

	// Need to reset the state that determines if we are in the cycle
	recipe->in_cycle = FALSE;

	// Need to reset the time array that tracks where each set is
	my_2d_free(recipe->n, (void **)recipe->t_arr);
	recipe->t_arr = (long int **)my_2d_calloc(recipe->n, recipe->m, sizeof(long int));

	// Need to cycle each of the block cdlls to their head nodes so that they
	// start in the correct positions.
	for (i = 0; i < recipe->n; i++) {
		for (j = 0; j < recipe->m; j++) {
			cdlln = (CDLL_NODE *)recipe->block_arr[i][j];
			while (cdlln->key != NULL) {
				cdlln = cdlln->next;
			}
			recipe->block_arr[i][j] = cdlln;
		}
	}
}

/** @} */





// LAYER FUNCTIONS

/**
 * \addtogroup layer
 * @{
 */

/**
 * \brief Creates a \ref layer
 * 
 * \param[in] n The size of n dimension of the \ref layer
 * \param[in] m The size of m dimension of the \ref layer
 *
 * \return The newly created \ref layer
 */
LAYER *qc_create_layer(int n, int m) {
	LAYER *layer;

	layer = (LAYER *)my_malloc(sizeof(LAYER));
	layer->num_blocks_pr = 0;
	layer->num_blocks_du = 0;
	layer->blocks = (BLOCK ***)my_2d_calloc(n, m, sizeof(BLOCK *));

	return layer;
}

/**
 * \brief Adds a \ref layer into a \ref recipe
 * 
 * \param[in] recipe The \ref recipe to insert the \ref layer into
 * \param[in] layer The \ref layer to insert 
 */
void qc_add_layer(RECIPE *recipe, LAYER *layer) {
	recipe->num_layers++;

	// Resize the recipe to make room for the new layer if necessary
	if (recipe->size < recipe->num_layers) {
		recipe->size = recipe->num_layers;
		recipe->layers = (LAYER **)my_realloc(recipe->layers, recipe->size*sizeof(LAYER *));
	}

	recipe->layers[recipe->num_layers - 1] = layer;
}

/**
 * \brief Compares two \link layer Layers\endlink for equality
 * 
 * \param[in] l1 The first \ref layer to compare	
 * \param[in] l2 The second \ref layer to compare  
 * \param[in] n The size of the n dimension of the \ref layer 
 * \param[in] m The size of the m dimension of the \ref layer 
 *
 * \return 1 if the \link layer Layers\endlink are identical, 0 otherwise
 */
int qc_compare_layers(LAYER *l1, LAYER *l2, int n, int m) {
	int i, j;

	// Compare two layers and detemrine if they are the same
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			if (l1->blocks[i][j] != l2->blocks[i][j]) {
				return FALSE;
			}
		}
	}
	return TRUE;
}

/**
 * \brief Determines if a given \ref layer is empty
 * 
 * Empty is defined as having no blocks, or having blocks which all have no
 * offsets.
 *
 * \param[in] layer The \ref layer to check for emptiness
 * \param[in] n The size of the n dimension of the \ref layer 
 * \param[in] m The size of the m dimension of the \ref layer 
 *
 * \return 1 if the \link layer is empty, 0 otherwise
 */
int qc_is_empty_layer(LAYER *layer, int n, int m) {
	int i, j, not_empty;

	not_empty = FALSE;
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			// If there is a block at a given location, check if there are offsets
			if (layer->blocks[i][j] && layer->blocks[i][j]->offsets != NULL) {
				not_empty = TRUE;
				break;
			}
		}
		if (not_empty) {
			break;
		}
	}
	return !not_empty;
}

/**
 * \brief Creates and builds a new \ref layer
 * 
 * \param[in] qc The \ref qc to build the \ref layer from
 * \param[in] recipe The \ref recipe to build the \ref layer in
 * \param[in] n The size of the n dimension of the \ref layer 
 * \param[in] m The size of the m dimension of the \ref layer  
 * \param[in] layer_big_t The big_t of the \ref layer
 *
 * \return The newly created and build \ref layer
 */
LAYER *qc_build_layer(QC *qc, RECIPE *recipe, int n, int m, long int layer_big_t) {
	LAYER *layer;

	layer = qc_create_layer(n, m);
	qc_create_and_insert_blocks(qc, recipe, layer, layer_big_t, PRIMAL);
	qc_create_and_insert_blocks(qc, recipe, layer, layer_big_t, DUAL);

	return layer;
}

/**
 * \brief Frees a \ref layer
 * 
 * \param[in] layer The \ref layer to free
 * \param[in] n The size of the n dimension of the \ref layer 
 */
void qc_free_layer(LAYER *layer, int n) {
	my_2d_free(n, (void **)layer->blocks);
	free(layer);
}

/** @} */



// BLOCK FUNCTIONS

/**
 * \addtogroup block
 * @{
 */

/**
 * \brief Creates a new \ref block
 * 
 * \param[in,out] recipe The \ref recipe to create the \ref block in
 * \param[in] ball The \ref ball to create the \ref block from 
 * \param[in] type The type of \ref block to create 
 *
 * \return The newly created \ref block
 */
BLOCK *qc_create_block(RECIPE *recipe, BALL *ball, int type) {
	int i, d_i, d_j, wt;
	long int d_t;
	HT *ht;
	STICK *stick;
	BLOCK *block;
	DLL_NODE *dlln;
	BALL *a, *b;
	OFFSET *offset, *found;

	block = (BLOCK *)my_malloc(sizeof(BLOCK));
	block->id = -1;
	block->type = type;
	block->offsets = NULL;

	assert(ball != NULL);

	block->i = INT_MAX; // Unused
	block->j = INT_MAX; // Unused
	block->t = LONG_MAX; // Unused
	block->big_t = LONG_MAX; // Unused
	block->hash_ijt = INT_MAX; // Unused

	// Loop all of the sticks in the ball
	ht = ball->stick_ht;
	for (i = 0; i < ht->length; i++) {
		dlln = ht->table[i];
		while (dlln != NULL) {
			stick = (STICK *)dlln->key;
			
			a = stick->a;
			b = stick->b;
		
			// Convert the stick to a weight
			wt = -log(stick->p_stick)*PRECISION;
			if (wt < PRECISION) {
				wt = PRECISION;
			}
			wt -= wt % 2;
			
			//
			// Create an offset corresponding to the stick
			//
			
			// If we are connected to the boundary we 
			// want to use a spacial coordinate for
			// the offset, as to signify which boundary
			// we are connected to.
			if (a->type == PRIMAL_BOUNDARY || a->type == DUAL_BOUNDARY) {
				offset = qc_create_offset(a->i, a->j, a->big_t, wt, a->type);
			}
			else if (b->type == PRIMAL_BOUNDARY || b->type == DUAL_BOUNDARY) {
				offset = qc_create_offset(b->i, b->j, b->big_t, wt, b->type);
			}

			// Calculate the delta coordinates based on
			// taking the value of the coordinates of the 
			// considered ball away from the ball at the
			// other end of the stick. That is, a - b, if 
			// the current ball is b
			else {
				d_i = (a == ball) ? b->i - a->i : a->i - b->i;
				d_j = (a == ball) ? b->j - a->j : a->j - b->j;
				d_t = (a == ball) ? b->big_t - a->big_t : a->big_t - b->big_t;

				offset = qc_create_offset(d_i, d_j, d_t, wt, a->type);
			}
			
			// Add the offset to the recipe
			recipe->offsets = qc_add_offset(recipe->offsets, offset);

			// See if the offset is in the recipe, if the found one is not the
			// one we just inserted, then the offset was not added, so we can
			// free it
			found = qc_offset_lookup(recipe->offsets, offset);
			if (found != offset) {
				free(offset);
			}

			// If an offset was found, insert it into the block
			if (found != NULL) {
				block->offsets = ll_insert_sorted(block->offsets, found, qc_offset_lt);
			}

			dlln = dlln->next;
		}
	}

	return block;
}

/**
 * \brief Creates a block and inserts it into an advanced recipe
 * 
 * \param[in,out] recipe The \ref recipe_adv to insert the \ref block into.
 * \param[in] ball The \ref ball to convert into a \ref block
 *
 * \return The newly created \ref block
 */
BLOCK *qc_create_block_adv(RECIPE_ADV *recipe, BALL *ball) {
	int i, d_i, d_j, wt;
	long int d_t;
	HT *ht;
	STICK *stick;
	BLOCK *block;
	DLL_NODE *dlln;
	BALL *a, *b;
	OFFSET *offset, *found;

	block = (BLOCK *)my_malloc(sizeof(BLOCK));
	block->id = -1;
	block->offsets = NULL;

	assert(ball != NULL);

	block->i = ball->i;
	block->type = ball->type;
	block->j = ball->j;
	block->t = ball->t;
	block->big_t = ball->big_t;
	block->hash_ijt = ball->hash_ijt;

	// Loop all of the sticks in the ball
	ht = ball->stick_ht;
	for (i = 0; i < ht->length; i++) {
		dlln = ht->table[i];
		while (dlln != NULL) {
			stick = (STICK *)dlln->key;

			a = stick->a;
			b = stick->b;

			// Convert the stick to a weight
			wt = -log(stick->p_stick)*PRECISION;
			if (wt < PRECISION) {
				wt = PRECISION;
			}
			wt -= wt % 2;

			//printf("[%d] ", wt);
			//printf("line a (%d, %d, %ld, %p)", a->i, a->j, a->big_t, a);
			//printf(" b (%d, %d, %ld, %p)\n", b->i, b->j, b->big_t, b);

			//
			// Create an offset corresponding to the stick
			//

			// If we are connected to the boundary we
			// want to use a spacial coordinate for
			// the offset, as to signify which boundary
			// we are connected to.
			if (a->type == PRIMAL_BOUNDARY || a->type == DUAL_BOUNDARY) {
				offset = qc_create_offset_adv(a->i, a->j, a->t, wt, a->type);
			}
			else if (b->type == PRIMAL_BOUNDARY || b->type == DUAL_BOUNDARY) {
				offset = qc_create_offset_adv(b->i, b->j, b->t, wt, b->type);
			}

			// Calculate the delta coordinates based on
			// taking the value of the coordinates of the
			// considered ball away from the ball at the
			// other end of the stick. That is, a - b, if
			// the current ball is b
			else {
				d_i = (a == ball) ? b->i - a->i : a->i - b->i;
				d_j = (a == ball) ? b->j - a->j : a->j - b->j;
				d_t = (a == ball) ? b->t - a->t : a->t - b->t;

				offset = qc_create_offset_adv(d_i, d_j, d_t, wt, a->type);
			}

			found = qc_offset_lookup_adv(recipe->offset_ht, offset->hash_ijt, offset);
			if (!found) {
				ht_insert_key(recipe->offset_ht, offset->hash_ijt, offset);
				found = offset;
			} else {
				// If we found the offset already, we don't need the one we
				// just made.
				free(offset);
			}

			/*** FILTER OFFSETS TO THE FUTURE ***/
			if (found->type != PRIMAL_BOUNDARY && found->type != DUAL_BOUNDARY && found->t > 0) {
				dlln = dlln->next;
				continue;
			}

			block->offsets = ll_insert_sorted(block->offsets, found, qc_offset_lt_adv);

			dlln = dlln->next;
		}
	}

	return block;
}

/**
 * \brief Create and insert \link block Blocks\endlink into a \ref layer
 * 
 * \param[in] qc The \qc that will be used with the \ref recipe
 * \param[in] recipe The \ref recipe of the \ref layer
 * \param[in] layer The \ref layer to create and insert into
 * \param[in] layer_big_t The big_t of the \ref layer
 * \param[in] type The type of the \ref layer (PRIMAL/DUAL) 
 */
void qc_create_and_insert_blocks(QC *qc, RECIPE *recipe, LAYER *layer, long int layer_big_t, int type) {
	NEST *nest;
	BALL *ball;
	BLOCK *block, *found;
	CDLL_NODE *current_ball;

	nest = (type == PRIMAL) ? qc->nest_pr : qc->nest_du;

	current_ball = nest->last_blocked_ball_cdlln->prev;

	// Loop over the latest layer of finalised balls.
	// Create any new blocks, and create a new layer
	// to store the blocks that have occured.

	ball = (BALL *)current_ball->key;

	// This assertion fails if there are no balls that haven't already been
	// blocked. i.e you're trying to create a layer but there are no balls to
	// put in it.
	assert(ball);

	while (ball->big_t == layer_big_t) {

		// Create the block and add it to the recipe,
		// add_block ignores the block if it already exists.
		block = qc_create_block(recipe, ball, type);
		recipe->blocks = qc_add_block(recipe->blocks, block);
	
		// Find the block in the blocks list, if the block found isn't the one
		// we inserted, then free the block (it wasn't inserted)
		found = qc_block_lookup(recipe->blocks, block);
		if (found != block) {
			qc_free_block(block);
		}

		// Insert the block found into the layer
		layer->blocks[ball->i][ball->j] = found;

		// Keep track of how many blocks are in the layer
		// helpful for knowing how many balls are created 
		// BLOCK *block;
		// for a given layer when converting the layer to
		// dots and lines.
		if (type == PRIMAL) {
			++layer->num_blocks_pr;
		} else {
			++layer->num_blocks_du;
		}

		nest->last_blocked_ball_cdlln = current_ball;
		
		current_ball = current_ball->prev;
		ball = (BALL *)current_ball->key;

		// Stop looping if we run out of balls and reach the sentinel node.
		if (!ball) {
			break;
		}
	}
}

/**
 * \brief Creates blocks inside an advanced recipe from any available finalized
 * balls in the provided qc. When a new ball is inserted, it also checks if
 * this new ball forms a repeating cycle.
 * 
 * \param[in] qc The \ref qc containing the finalized balls
 * \param[out] recipe The \ref recipe_adv to insert the blocks into
 * \param[in] type The type of nest to search for balls to be converted. 
 *
 * \returns True if the blocks inserted have formed a repeating cycle, and
 * hence indicates the recipe has stabilized. False otherwise.
 */
int qc_create_and_insert_blocks_adv(QC *qc, RECIPE_ADV *recipe, int type, long int thresh_big_t, int thresh_cycle_len) {
	int cycle, count;
	long int big_t;
	NEST *nest;
	BALL *ball;
	BLOCK *block, *block2;
	OFFSET *offset;
	LL_NODE *ll;
	CDLL_NODE *current_ball, *n1, *n2, *node, *old_node;

	nest = (type == PRIMAL) ? qc->nest_pr : qc->nest_du;

	// Get the next available ball to be blockified
	current_ball = nest->last_blocked_ball_cdlln->prev;
	ball = (BALL *)current_ball->key;
	assert(ball);

	// Get big_t from the first ball we will be checking.
	big_t = ball->big_t;
	while (ball && ball->big_t == big_t) {
		// Create and insert the block into the recipe in in the hash table
		// (based on i,j,t coords) and the cdll.
		block = qc_create_block_adv(recipe, ball);
		ht_insert_key(recipe->block_ht, block->hash_ijt, block);
		cdll_insert_tail(recipe->block_cdll, block, NULL);
		cdll_insert_tail(recipe->block_arr[ball->i][ball->j], block, NULL);

		// We want to find all other blocks identical to this one (in position
		// and content). We do this because we want to traverse backward from
		// each of these duplicates and compare each block along the way to see
		// if we have a sequence that has repeated.
		//
		// A sequence has repeated if we travel from the new node to the old
		// node and all of the blocks match. Think of this as comparing
		// old_node to new_node, then old_node - 1 to new_node - 1, then
		// old_node - 2 to new_node - 2, and so forth until new_node - x is the
		// same position of old_node. We therefore have a cycle of length
		// new_node - old_node.
		cycle = FALSE;
		node = recipe->block_cdll->prev->prev;
		while (node != recipe->block_cdll) {
			n1 = recipe->block_cdll->prev;
			n2 = node;
	
			// Save where we started the comparison so we know when we have
			// traversed the length of the section to be checked 
			old_node = n2;

			// Keep looping from the backwards from old_node and new_node and
			// comparing the blocks along the way for equality. If either hits
			// the head node of the sequence (key == NULL) then we can't have a
			// repeating sequence.
			count = 0;
			while (n1->key != NULL && n2->key != NULL && qc_block_eq_adv((BLOCK *)n1->key, (BLOCK *)n2->key)) {
				count++;
				n1 = n1->prev;
				n2 = n2->prev;

				// If we have detected a cyclic sequence, then break. This is
				// defined as a sequence of at least thresh_len nodes in
				// length, and that we have balls of at least thresh_big_t
				// considered before looking for a sequence. The sequence, as
				// described above, exists if all nodes match when traversing
				// backwards from old_node and new_node until new_node
				// traversing pointer is equal to the starting position of the
				// old node traversing pointer.
				if (n1 == old_node) {
					if (big_t >= thresh_big_t && count >= thresh_cycle_len) {
						cycle = TRUE;
					}
					// We want to stop searching at this point whether a cycle
					// is found or not as our search domain has ended.
					break;
				}
			}

			// If we have a cycle, we want to stop and break out of the loop.
			if (cycle == TRUE) {
				break;
			}
			node = node->prev;
		}

		nest->last_blocked_ball_cdlln = current_ball;

		// If we have a cycle, we want to stop and break out of the loop.
		if (cycle == TRUE) {
			break;
		}
		current_ball = current_ball->prev;
		ball = (BALL *)current_ball->key;
	}

	if (cycle) {
		// Calculate the period of the recipe, that is, the time between the
		// end of the recipe, and it's corresponding block the cycle before.
		block = (BLOCK *)recipe->block_cdll->prev->key;
		block2 = (BLOCK *)old_node->key;
		recipe->cycle_period = block->t - block2->t;

		// Calculate the t0 of the recipe. That is, the lowest time of any
		// block in the cycling part of the recipe.
		node = recipe->block_cdll->prev;
		while (node != old_node) {
			block = (BLOCK *)node->key;
			if (block->t < recipe->cycle_t0) {
				recipe->cycle_t0 = block->t;
			}
			node = node->prev;
		}

		// This is used to determine when to start using the cyclic region
		recipe->cycle_end = recipe->block_cdll->prev;

		// These values are mostly useful for printing the cyclic region of the
		// recipe and aren't actually required.
		recipe->cycle_begin = old_node->next;
		recipe->cycle_len = count;

		node = recipe->cycle_begin;
		while (node != recipe->block_cdll) {
			block = (BLOCK *)node->key;
			if (block == NULL) {
				continue;
			}

			ll = block->offsets;
			while (ll != NULL) {
				offset = (OFFSET *)ll->key;
				if (offset->t != 0) {
					ll = ll->next;
					continue;
				}

				if (block->type == PRIMAL && offset->j == 0 && offset->wt < recipe->min_horz_wt_pr) {
					recipe->min_horz_wt_pr = offset->wt;
				}
				if (block->type == DUAL && offset->i == 0 && offset->wt < recipe->min_horz_wt_du) {
					recipe->min_horz_wt_du = offset->wt;
				}

				ll = ll->next;
			}
			node = node->next;
		}

		//qc_print_recipe_adv(recipe);	
	}

	return cycle;
}

/**
 * \brief Adds a \ref block to a \ref llist of \link block Blocks\endlink. If
 * the \ref block already exists, then do nothing.
 * 
 * \param[in] blocks_head The head node of the \ref llist
 * \param[in] block The \ref block to look for
 *
 * \return The found \ref block if the \ref block already existed, otherwise the newly
 * inserted \ref block.
 */
LL_NODE *qc_add_block(LL_NODE *blocks_head, BLOCK *block) {
	BLOCK *found;

	found = qc_block_lookup(blocks_head, block);

	// If the block wasn't found, then insert the block
	if (found == NULL) {
		block->id = (blocks_head == NULL) ? 0 : ((BLOCK *)blocks_head->key)->id + 1;
		blocks_head = ll_insert(blocks_head, block);
	} 

	return blocks_head;
}

/**
 * \brief Look for a \ref block in a \ref llist of \link block Blocks\endlink.
 * 
 * \param[in] blocks_head The head node of the \ref llist
 * \param[in] block The \ref block to look for
 *
 * \return The found \ref block if the \ref block was found, NULL otherwise
 */
BLOCK *qc_block_lookup(LL_NODE *blocks_head, BLOCK *block) {
	LL_NODE *node;
	
	node = blocks_head;
	while (node != NULL) {
		if (qc_block_eq((BLOCK *)node->key, block)) {
			return (BLOCK *)node->key;
		}
		node = node->next;
	}

	return NULL;
}

/**
 * \brief Compares two \link block Blocks\endlink for equality
 * 
 * \param[in] b1 The first \ref block to compare
 * \param[in] b2 The second \ref block to compare
 *
 * \return 1 if b1 is equal to b2, 0 otherwise
 */
int qc_block_eq(BLOCK *b1, BLOCK *b2) {
	return ll_eq( b1->offsets, b2->offsets, qc_offset_eq );
}

/**
 * \brief Compares two \link block Blocks\endlink for equality
 *
 * Differs from the regular equaltiy as advanced blocks are equal if their i &
 * j coordinates match, as well as their offsets, rather than just their
 * offsets.
 * 
 * \param[in] b1 The first \ref block to compare
 * \param[in] b2 The second \ref block to compare
 *
 * \return 1 if b1 is equal to b2, 0 otherwise
 */
int qc_block_eq_adv(BLOCK *b1, BLOCK *b2) {
	return b1->i == b2->i && b1->j == b2->j && ll_eq( b1->offsets, b2->offsets, qc_offset_eq_adv );
}

/**
 * \brief Frees a \ref block
 * 
 * \param[in] block The \ref block to be freed
 *
 * \sa qc_free_void_block
 */
void qc_free_block(BLOCK *block) {
	ll_free(block->offsets, NULL);
	free(block);
}

/**
 * \brief Frees a \ref block. Accepts a void pointer to a \ref block.
 * 
 * \param[in] block The \ref block to be freed
 *
 * \sa qc_free_block
 */
void qc_free_void_block(void *block) {
	qc_free_block((BLOCK *)block);
}

/** @} */







// OFFSET FUNCTIONS

/**
 * \addtogroup offset
 * @{
 */

/**
 * \brief Create a new \ref offset
 * 
 * \param[in] i The i coodinate of the \ref offset
 * \param[in] j The j coordinate of the \ref offset
 * \param[in] big_t The big_t of the \ref offset
 * \param[in] wt The weight of the \ref offset
 * \param[in] type The type of the \ref offset
 *
 * \return The newly created \ref offset
 */
OFFSET *qc_create_offset(int i, int j, long int big_t, double wt, int type) {
	OFFSET *offset;

	offset = (OFFSET *)my_malloc(sizeof(OFFSET));
	offset->id = -1;
	offset->i = i;
	offset->j = j;
	offset->t = LONG_MAX;
	offset->big_t = big_t;
	offset->wt = wt;
	offset->type = type;
	offset->hash_ijt = INT_MAX;

	return offset;
}

/**
 * \brief Create a new \ref offset, for an advanced recipe
 * 
 * \param[in] i The i coodinate of the \ref offset
 * \param[in] j The j coordinate of the \ref offset
 * \param[in] t The t of the \ref offset
 * \param[in] big_t The big_t of the \ref offset
 * \param[in] wt The weight of the \ref offset
 * \param[in] type The type of the \ref offset
 *
 * \return The newly created \ref offset
 */
OFFSET *qc_create_offset_adv(int i, int j, long int t, double wt, int type) {
	OFFSET *offset;

	offset = (OFFSET *)my_malloc(sizeof(OFFSET));
	offset->id = -1;
	offset->i = i;
	offset->j = j;
	offset->t = t;
	offset->big_t = LONG_MAX;
	offset->wt = wt;
	offset->type = type;

	offset->hash_ijt = qc_hash_ijt(i, j, t);

	return offset;
}

/**
 * \brief Adds an \ref offset to a \ref llist of \link offset Offsets\endlink
 * 
 * \param[in] offset_head The head node of a \ref llist containing \link offset
 * Offsets\endlink
 * \param[in] offset The \ref offset to add 
 *
 * \return The head of the \ref llist
 */
LL_NODE *qc_add_offset(LL_NODE *offset_head, OFFSET *offset) {
	OFFSET *found;

	// If the offset isn't a boundary, and has a positive big_t, do nothing.
	if (offset->type != PRIMAL_BOUNDARY && offset->type != DUAL_BOUNDARY && offset->big_t > 0) {
		return offset_head;
	}

	// Look for the offset in the list
	found = qc_offset_lookup(offset_head, offset);

	// If the offset does not exist, look for its inverse
	if (found == NULL) {
		offset->i *= -1;
		offset->j *= -1;
		offset->big_t *= -1;
	
		found = qc_offset_lookup(offset_head, offset);
	
		// If the inverse is not found either, then create the new offset
		if (found == NULL) {	
			offset->id = (offset_head == NULL) ? 0 : ((OFFSET *)offset_head->key)->id + 1;
			offset_head = ll_insert(offset_head, offset);
		} 
	
		offset->i *= -1;
		offset->j *= -1;
		offset->big_t *= -1;
	} 

	return offset_head;
}

/**
 * \brief Finds an \ref offset in a \ref block
 * 
 * \param[in] block The \ref block to search for the \ref offset in
 * \param[in] i The i coordinate of the \ref offset
 * \param[in] j The j coordinate of the \ref offset
 * \param[in] big_t The big_t of the \ref offset
 * \param[in] wt The wt of the \ref offset
 *
 * \return The \ref offset if found, NULL otherwise.
 */
OFFSET *qc_find_offset(BLOCK *block, int i, int j, long int big_t, double wt) {
	LL_NODE *n;
	OFFSET *o;

	n = block->offsets;
	while (n != NULL) {
		o = (OFFSET *)n->key;
		if (o->i == i && o->j == j && o->big_t == big_t && o->wt == wt) {
			return o;
		}
		n = n->next;
	}

	return NULL;
}

/**
 * \brief Looks for an \ref offset in a \ref llist
 * 
 * \param[in] offset_head The head node of the \ref llist
 * \param[in] offset The \ref offset to locate
 *
 * \return The offset if found, NULL otherwise.
 */
OFFSET *qc_offset_lookup(LL_NODE *offset_head, OFFSET *offset) {
	LL_NODE *node;
	
	node = offset_head;
	while (node != NULL) {
		if ( qc_offset_eq((OFFSET *)node->key, offset) ) {
			return (OFFSET *)node->key;
		}
		node = node->next;
	}

	return NULL;
}

/**
 * \brief Looks for an \ref offset in a \ref ht as used in advanced recipes
 * 
 * \param[in] ht The \ref ht to look up
 * \param[in] hash  The hash of the key \ref offset to look for
 * \param[in] offset The offset to look for
 *
 * \return The offset if found, otherwise NULL.
 */
OFFSET *qc_offset_lookup_adv(HT *ht, int hash, OFFSET *offset) {
	DLL_NODE *node;

	node = ht_hash_lookup(ht, hash);
	while (node != NULL) {
		if ( qc_offset_eq_adv((OFFSET *)node->key, offset) ) {
			return (OFFSET *)node->key;
		}
		node = node->next;
	}
	return NULL;
}

/**
 * \brief Compares two \link offset Offsets\endlink for equality
 * 
 * \param[in] k1 The first \ref offset for comparison
 * \param[in] k2 The second \ref offset for comparison
 *
 * \return 1 if the two are equal, 0 otherwise.
 */
int qc_offset_eq(void *k1, void *k2) {
	OFFSET *n1, *n2;

	if (k1 == NULL || k2 == NULL) {
		return 0;
	}

	n1 = (OFFSET *)k1;
	n2 = (OFFSET *)k2;

	if (n1->i == n2->i && n1->j == n2->j && n1->big_t == n2->big_t && n1->wt == n2->wt && n1->type == n2->type) {
		return 1;
	}
	return 0;
}

/**
 * \brief Compares two \link offset Offsets\endlink for equality as used by
 * advanced recipes
 *
 * Differs from the regular comparison as it uses t rather than big_t.
 * 
 * \param[in] k1 The first \ref offset for comparison
 * \param[in] k2 The second \ref offset for comparison
 *
 * \return 1 if the two are equal, 0 otherwise.
 */
int qc_offset_eq_adv(void *k1, void *k2) {
	OFFSET *n1, *n2;

	if (k1 == NULL || k2 == NULL) {
		return 0;
	}

	n1 = (OFFSET *)k1;
	n2 = (OFFSET *)k2;

	if (n1->i == n2->i && n1->j == n2->j && n1->t == n2->t && n1->wt == n2->wt && n1->type == n2->type) {
		return 1;
	}
	return 0;
}

/**
 * \brief Compares two \link offset Offsets\endlink to see if the first is less
 * than the second.
 * 
 * \param[in] k1 The first \ref offset for comparison
 * \param[in] k2 The second \ref offset for comparison
 *
 * \return 1 if k1 is less than k2, 0 otherwise.
 */
int qc_offset_lt(void *k1, void *k2) {
	OFFSET *n1, *n2;

	n1 = (OFFSET *)k1;
	n2 = (OFFSET *)k2;

	if (k1 == NULL) {
		return 0;
	}
	if (k2 == NULL) {
		return 1;
	}

	if (n1->i < n2->i) {
		return 1;
	}
	else if (n1->i == n2->i) {
		if (n1->j < n2->j) {
			return 1;
		}
		else if (n1->j == n2->j) {
			if (n1->big_t < n2->big_t) {
				return 1;
			}
		}
	}

	return 0;
}

/**
 * \brief Compares two \link offset Offsets\endlink to see if the first is less
 * than the second, as used by advanced recipes.
 *
 * Differs from the regular comparison as it uses t instead of big_t and also
 * sorts based on the type of the offset.
 * 
 * \param[in] k1 The first \ref offset for comparison
 * \param[in] k2 The second \ref offset for comparison
 *
 * \return 1 if k1 is less than k2, 0 otherwise.
 */
int qc_offset_lt_adv(void *k1, void *k2) {
	OFFSET *n1, *n2;

	n1 = (OFFSET *)k1;
	n2 = (OFFSET *)k2;

	if (k1 == NULL) {
		return 0;
	}
	if (k2 == NULL) {
		return 1;
	}

	if (n1->i < n2->i) {
		return 1;
	}
	else if (n1->i == n2->i) {
		if (n1->j < n2->j) {
			return 1;
		}
		else if (n1->j == n2->j) {
			if (n1->t < n2->t) {
				return 1;
			} else if (n1->t == n2->t) {
				if (n1->type < n2->type) {
					return 1;
				}
			}
		}
	}

	return 0;
}

/**
 * \brief Frees an \ref offset
 * 
 * \param[in] offset The \ref offset to be freed
 */
void qc_free_offset(OFFSET *offset) {
	free(offset);
}

/** @} **/







// BOOTUP FUNCTIONS

/**
 * \brief Boot up a \ref qc according to a \ref recipe.
 * 
 * \param[in] qc The \ref qc to boot up
 * \param[in] recipe The \ref recipe to use to boot up
 * \param[in] n The size of the n dimension of the \ref recipe
 * \param[in] m The size of the m dimension of the \ref recipe
 * \param[in] switch_time The time layer to switch off the boot up process and
 * determine the repeating layer
 *
 * \return 1 if the boot up process completed, 0 otherwise
 */
int qc_boot_up(QC *qc, RECIPE *recipe, int n, int m, long int switch_time) {
	if (recipe->type == RECIPE_INFINITE) {
		return qc_boot_up_infinite(qc, recipe, n, m, switch_time);
	} else {
		printf("Unsupported recipe type.\n");
		exit(1);
	}
}

/**
 * \brief Bood up a \ref qc according to a \ref recipe_adv
 * 
 * \param[in] qc The \ref qc to boot up
 * \param[in] recipe The \ref recipe_adv to use to boot up
 * \param[in] thresh_big_t The minimum big_t that must pass before a cycle is considered. 
 *
 * \return 1 if a cycle has been detected, 0 otherwise.
 */
int qc_boot_up_adv(QC *qc, RECIPE_ADV *recipe, long int thresh_big_t, int thresh_cycle_len) {
	int pr, du;

	// Cannot do anything until there are finalized balls
	if (qc->unfinalized_big_t <= 0) return 0;

	du = qc_create_and_insert_blocks_adv(qc, recipe, DUAL, thresh_big_t, thresh_cycle_len);
	pr = qc_create_and_insert_blocks_adv(qc, recipe, PRIMAL, thresh_big_t, thresh_cycle_len);

	return (pr || du);
}

/**
 * \brief Boot up a \ref qc using an infinite (repeating) computation
 * 
 * \param[in] qc The \ref qc to boot up
 * \param[in] recipe The \ref recipe to use to boot up
 * \param[in] n The size of the n dimension of the \ref recipe
 * \param[in] m The size of the m dimension of the \ref recipe
 * \param[in] switch_time The time layer to switch off the boot up process and
 * determine the repeating layer
 *
 * \return 1 if the boot up process completed, 0 otherwise
 */
int qc_boot_up_infinite(QC *qc, RECIPE *recipe, int n, int m, long int switch_time) {
	LL_NODE *ll;
	LAYER *layer;
	BLOCK *block;
	OFFSET *offset;
	int num_layers, i, j;
	long int layer_big_t;

	num_layers = recipe->num_layers;

	// Cannot create layers until there are finalized balls
	if (qc->unfinalized_big_t <= 0) return 0;

	// Define a new layer for the most recently finalized time step
	layer_big_t = qc->unfinalized_big_t - 1;
	layer = qc_build_layer(qc, recipe, n, m, layer_big_t);
	
	// If errors are off, and we have an empty layer, the boot up is finished
	if (qc->perfect_gates && qc_is_empty_layer(layer, recipe->n, recipe->m)) {

		qc_free_layer(layer, recipe->n);
		
		// Turn perfect gates off as we have finished the bootup process
		qc->perfect_gates = FALSE;
		
		// There is no point in turning off the tracking of errors here as
		// the qc will need to be destroyed at the end of the boot-up phase.
		// qc->track = FALSE;
		
		layer = recipe->layers[recipe->repeat_id];
		for (i = 0; i < recipe->n; i++) {
			for (j = 0; j < recipe->m; j++) {
				block = layer->blocks[i][j];
				if (block == NULL) {
					continue;
				}

				ll = block->offsets;
				while (ll != NULL) {
					offset = (OFFSET *)ll->key;
					if (offset->big_t != 0) {
						ll = ll->next;
						continue;
					}

					if (block->type == PRIMAL && offset->j == 0 && offset->wt < recipe->min_horz_wt_pr) {
						recipe->min_horz_wt_pr = offset->wt;
					}
					if (block->type == DUAL && offset->i == 0 && offset->wt < recipe->min_horz_wt_du) {
						recipe->min_horz_wt_du = offset->wt;
					}

					ll = ll->next;
				}
			}
		}

		return DONE;
	} 
	
	// When qc->big_t is equal to the switch time, we want to turn off errors.
	//
	// Want >= because at very small switch_time, qc->big_t can be larger than
	// switch time before unfinalized_big_t is greater than 0.
	if (qc->perfect_gates == FALSE && qc->big_t >= switch_time) {
		qc->perfect_gates = TRUE;
	}

	// When layer_big_t catches up to switch time, we want to use the previous
	// layer in order to determine the repeating layer. We then want to delete 
	// all replicates of that layer before it, and store the repeat_id.
	if (layer_big_t == switch_time - 1) {
		// Go backward through the layers and delete all that are the same as
		// the repeating one
		for (i = num_layers - 1; i >= 0; --i) {
			if (qc_compare_layers(layer, recipe->layers[i], recipe->n, recipe->m)) {
				qc_free_layer(recipe->layers[i], recipe->n);
				recipe->num_layers--;
				recipe->size--;
			} else {	
				break;
			}
		}

		qc_add_layer(recipe, layer);
		recipe->repeat_id = recipe->num_layers - 1;
		return 0;
	} 
	
	qc_add_layer(recipe, layer);
	return 0;
}

// CONVERT FUNCTIONS

/**
 * \brief Takes a ball and returns a hash of its i, j, t coordinates
 * 
 * \param[in] ball The \ref ball to get the hash for
 *
 * \returns The 32-bit hash corresponding to the ball's position
 */
int qc_hash_ball(BALL *ball) {
	return qc_hash_ijt(ball->i, ball->j, ball->t);
}

/**
 * \brief Takes an i, j, t cordinate and returns a hash
 * 
 * \param[in] i The i coordinate
 * \param[in] j The j coordinate
 * \param[in] t The t coordinate
 *
 * \returns A 32-bit hash of the coordinate.
 */
int qc_hash_ijt(int i, int j, long int t) {
	long int pos[3];

	pos[0] = i;
	pos[1] = j;
	pos[2] = t;

	return abs(fasthash32(pos, 3, 0));
}

/**
 * \brief Gets the block for a given set
 * 
 * \param[in] qc The qc with the recipe 
 * \param[in] set The set the block is required for 
 */
BLOCK *qc_get_block_for_set(QC *qc, SET *set) {
	BLOCK *block, *block2;
	RECIPE_ADV *recipe;
	CDLL_NODE *node, *n;
	int found, cycled;

	recipe = qc->recipe_adv;
	if (!recipe) return NULL;

	node = recipe->block_arr[set->i][set->j];
	if (node->next->key != NULL) {
		node = node->next;
		recipe->block_arr[set->i][set->j] = node;
	} else {
		// We want to start at the beginning of the cycle, then traverse down
		// until the first occurence of a block at this position. 
		found = FALSE;
		n = qc->recipe_adv->cycle_begin;
		while (n != qc->recipe_adv->block_cdll) {
			block = (BLOCK *)n->key;
			if (set->i == block->i && set->j == block->j) {
				found = TRUE;
				break;
			}
			n = n->next;
		}
		assert(found);
		
		// Once we have the block we want to move to, we need to find it in the
		// block_arr cdll. We do this by looping forward from the current
		// position. The justification for this, is that for long repeating
		// cycles, this will be faster than going backwards (which would be
		// faster for very short repeating cycles of 1 'layer' or so). 
		found = FALSE;
		cycled = FALSE;
		while (cycled == FALSE) {
			if (node->key == NULL) {
				node = node->next;	
				continue;
			}
			
			block2 = (BLOCK *)node->key;
			if (block == block2) {
				found = TRUE;
				break;
			}

			if (node == recipe->block_arr[set->i][set->j]) {
				// This avoids an infinite loop in the case of an error. It
				// should be impossible to actually loop the entire block_arr
				// cdll and not find the block we are looking for.
				if (!cycled) {
					cycled = TRUE;
				} else {
					break;
				}
			}

			node = node->next;
		}
		assert(found);

		recipe->block_arr[set->i][set->j] = node;
	}

	return (BLOCK *)node->key;
}

/**
 * \brief Converts a given block to a dot and respective lines using an
 * advanced recipe. 
 * 
 * \param[in] qc The \ref qc to get the recipe from
 * \param[out] matching The \ref matching to insert the dot and lines into
 * \param[in] ball The \ref ball to convert to the dot
 * \param[in] block The \ref block to use the offsets from to construct the
 * lines 
 */
void qc_convert_block_to_dot_and_lines(QC *qc, MATCHING *matching, BALL *ball, BLOCK *block) {
	DOT *dot1, *dot2;
	OFFSET *o;
	BALL *b;
	LINE *line;
	LL_NODE *node, *node2;
	DLL_NODE *dnode;
	CDLL_NODE *n;
	int ii, jj, found;
	long int tt;

	assert(block != NULL);
	assert(ball->dot == NULL);

	// Note this is called before measurements, so we never want to insert a
	// vertex right away.
	dot1 = m_create_dot(ball, NULL);
	ht_insert_key(matching->dot_ht, qc_hash_ball(ball), dot1);
	if (matching->undo_flag == TRUE) {
		n = cdll_create_node(dot1);
		m_create_undo(matching, UNDO_CREATE_DOT, NULL, 0, NULL, NULL, n);	
	}

	assert(dot1);

	// Loop over all the offsets eminating from the current block
	node = block->offsets;
	while (node != NULL) {
		o = (OFFSET *)node->key;

		// If the offset is a boundary, use the user defined get_boundary
		// function inside qc to find the relevant dot.
		if (o->type == PRIMAL_BOUNDARY || o->type == DUAL_BOUNDARY) {
			b = qc->get_boundary(o->i, o->j, o->big_t, o->type, qc->boundaries);
			assert(b != NULL);
			if (b->dot == NULL) {
				dot2 = m_create_and_insert_dot_and_vertex(matching, b);
			} else {
				dot2 = b->dot;
			}
			found = TRUE;
		}

		// Connect to a dot in the current time layer or the previous depending
		// on the offset big_t.
		else {
			ii = ball->i + o->i;
			jj = ball->j + o->j;
			tt = ball->t + o->t;

			found = FALSE;
			dnode = ht_hash_lookup(matching->dot_ht, qc_hash_ijt(ii, jj, tt));
			while (dnode != NULL) {
				dot2 = (DOT *)dnode->key;
				if (dot2->i == ii && dot2->j == jj && dot2->little_t == tt) {
					while (dot2->merge != NULL) {
						dot2 = dot2->merge;
					}
					found = TRUE;
					break;
				}
				dnode = dnode->next;
			}
		}

		// If the destination doesn't exist, then we don't want to use this
		// offset to create a new line
		if (!found) {
			node = node->next;
			continue;
		}

		// Check to see if the line already exists from the destination dot.
		// Note: This code never actually appears to succeed...
		found = FALSE;
		node2 = dot2->lines;
		while (node2 != NULL) {
			line = (LINE *)node2->key;

			// If the line exists, skip to the next offset
			if ((line->a == dot1) || (line->b == dot1)) {
				node2 = node2->next;
				found = TRUE;
				break;
			}
			node2 = node2->next;
		}

		if (!found) {
			// If the line doesn't exist, create it.
			m_create_line_wt(matching, dot1, dot2, o->wt - o->wt%2);
		}
		
		node = node->next;
	}
}

/**
 * \brief Converts balls to dots 
 * 
 * \param[in] qc The \ref qc containing the recipe to use 
 * \param[in] matching The \ref matching to insert the \link dot Dots\endlink into 
 * \param[in] ball_cdlln A ::CDLL_NODE to traverse backwards from to obtain the
 * relevant \link ball Balls\endlink, usually last_converted_ball_cdlln
 * \param[in] num_blocks The number of \ref block that should be traversed
 */
void qc_convert_balls_to_dots(QC *qc, MATCHING *matching, CDLL_NODE *ball_cdll, int num_blocks) {
	CDLL_NODE *cdll_node;
	BALL *ball;
	DOT *dot;
	int i, n, m;

	// If there are no balls left (i.e the list is empty or we reached the
	// end), simply return
	if (ball_cdll->prev->key == NULL) {
		return;
	}

	// ball_cdll is actually a last_converted_ball_cdlln, therefore
	// we want to traverse backwards from the last converted until 
	// we have all the relevant balls to the current layer, the 
	// number of balls is the number of blocks in the layer.
	cdll_node = ball_cdll;
	for (i = 0; i < num_blocks; i++) {
		cdll_node = cdll_node->prev;
		ball = (BALL *)cdll_node->key;

		assert(ball->dot == NULL);
		if (ball->mp < 0) {
			dot = m_create_and_insert_dot_and_vertex(matching, ball);
		} else {
			dot = m_create_and_insert_dot(matching, ball);
		}

		n = ball->i;
		m = ball->j;

		// Store the dot in the recipe, this allows for the easy 
		// creation of lines across time.
		qc->recipe->dotarr2[n][m] = dot;
	}
}

/**
 * \brief Converts the \link offset Offsets\endlink from a \ref qc to \link
 * line Lines\endlink
 * 
 * \param[in] qc The \ref qc containing the recipe to use 
 * \param[in] matching The \ref matching to insert the \link line Lines\endlink into 
 * \param[in] nest The \ref nest that the \link ball Balls\endlink belong to
 * \param[in] ball_cdll The \ref cdll containing the \link ball Balls\endlink
 * to be used to convert to lines
 * \param[in] layer The \ref layer of the \ref recipe to convert the \link
 * offset Offsets\endlink from
 * \param[in] num_blocks The number of \ref block that should be traverse
 */
void qc_convert_offsets_to_lines(QC *qc, MATCHING *matching, NEST *nest, CDLL_NODE *ball_cdll, LAYER *layer, int num_blocks) {
	CDLL_NODE *cdll_node;   
	LL_NODE *node;
	DOT *dot1, *dot2;
	BLOCK *block;
	BALL *a, *b;
	OFFSET *o;
	int i, n, m;

	// If there are no balls left (i.e the list is empty or we reached the
	// end), simply return
	if (ball_cdll->prev->key==NULL) {
		return;
	}

	// ball_cdll is actually a last_converted_ball_cdlln, therefore
	// we want to traverse backwards from the last converted until 
	// we have all the relevant balls to the current layer, the 
	// number of balls is the number of blocks in the layer.
	cdll_node = ball_cdll;
	for (i = 0; i < num_blocks; i++) {
		cdll_node = cdll_node->prev;
		a = (BALL *)cdll_node->key;
		dot1 = a->dot;
		assert(dot1 != NULL);

		n = a->i;
		m = a->j;

		// Get the block that is relevant to the ball currently being considered
		block = layer->blocks[n][m];
		assert(block != NULL);

		// Loop over all the offsets eminating from the current block
		node = block->offsets;
		while (node != NULL) {
			o = (OFFSET *)node->key;
			
			// If the offset is a boundary, use the user defined get_boundary 
			// function inside qc to find the relevant dot.
			if (o->type == PRIMAL_BOUNDARY || o->type == DUAL_BOUNDARY) {
				b = qc->get_boundary(o->i, o->j, o->big_t, o->type, qc->boundaries);
				if (b == NULL) {
					//qc_print_offset((void *)o);
					assert(b != NULL);
				}
				if (b->dot == NULL) {
					m_create_and_insert_dot_and_vertex(matching, b);
				}
				dot2 = b->dot;
			} 
			
			// Connect to a dot in the current time layer or the previous depending
			// on the offset big_t.
			else if (o->big_t == 0) {
				dot2 = qc->recipe->dotarr2[n + o->i][m + o->j];
			} 
			else {
				dot2 = qc->recipe->dotarr1[n + o->i][m + o->j];
			}
		
			assert(dot2 != NULL);
			m_create_and_insert_line_wt(matching, dot1, dot2, o->wt - o->wt%2);
			
			node = node->next;
		}
	}

	// Update the last_converted_ball_cdlln in the nest now that 
	// we have converted more balls into dots and lines
	nest->last_converted_ball_cdlln = cdll_node;
}

// CONVERT INFINITE FUNCTIONS

/**
 * \brief Gets the next \ref layer in an infinite \ref recipe computation
 * 
 * \param[in] qc The \ref qc being used for the computation
 *
 * \return The next ]ref layer to be used for the computation
 */
LAYER *qc_get_next_layer_infinite(QC *qc) {
	RECIPE *recipe;
	BALL *ball;
	int repeat_id, layer_id;
	long int big_t;
	
	recipe = qc->recipe;
	repeat_id = recipe->repeat_id;

	// repeat_id = 2, 5 layers
	// big_t = 0, so layers[0]
	// big_t = 1, so layers[1]
	// big_t = 2, so layers[2], repeated count = 0
	// big_t = 3, so layers[2], repeated count = 1 
	// big_t = 4, so layers[2], repeated_count = 2
	// big_t = 5, so layers[3] == big_t - repeated_count
	// big_t = 6, so layers[4]

	assert(repeat_id < recipe->num_layers);

	if (repeat_id == -1) {
		printf("No repeating layer was discovered.\n");
		printf("You may not have: \n");
		printf("\t- Performed the boot up process correctly;\n");
		printf("\t- Set switch_time to be a large enough value;\n");
		printf("\t- Used the correct form of boot up for your simulation.\n");
		assert(repeat_id >= 0);
	}
 
	// Find the current big_t by looking in the balls
	ball = (BALL *)qc->nest_pr->ball_cdll->next->key;
	if (ball == NULL) {
		ball = (BALL *)qc->nest_du->ball_cdll->next->key;
	}
	assert(ball != NULL);
	big_t = ball->big_t;

	// If big_t is less than the repeat_id, then we have unique layers
	if (big_t <= repeat_id) {
		layer_id = big_t;		
	}
	else 
	{
		// If perfect gates are off, then we want to repeat the repeated layer
		if (qc->perfect_gates == FALSE) {
			layer_id = repeat_id;
			recipe->repeated_count++;
		}

		// Otherwise we are at the end of the computation, use the last of the
		// layers
		else {
			layer_id = big_t - recipe->repeated_count;

			if (layer_id >= recipe->num_layers) {
				printf("Failure: Attempted to use more layers than available.\n");
				printf("Hint: It is likely that measure_stabilizers is being called too many times after turning on perfect_gates.\n");
				exit(1);
			}
		}
	}
	
	return recipe->layers[layer_id];
}

/**
 * \brief Converts the \link nest Nests\endlink from a \ref qc for an infinite
 * computation
 * 
 * \param[in] qc The \ref qc containing the \link nest Nests\endlink to convert
 * \param[in] undo Whether or not undo is enabled for matching 
 */
void qc_convert_nest_infinite(QC *qc, int undo) {
	LAYER *layer;
	DOT ***dotarrt;
	RECIPE *recipe;
	CDLL_NODE *uc_pr, *uc_du;

	// Set the undo flag in the matchings
	qc->m_pr->undo_flag = undo;
	qc->m_du->undo_flag = undo;

	// Rotate the dotarrs such that the previously considered layer
	// of dots is moved into dotarr1 so that dotarr2 can be overwritten.
	recipe = qc->recipe;
	dotarrt = recipe->dotarr1;
	recipe->dotarr1 = recipe->dotarr2;
	recipe->dotarr2 = dotarrt;

	layer = qc_get_next_layer_infinite(qc);

	// We want to convert balls to dots, and offsets to lines from the same 
	// starting location in the ball_cdll, so store this pointer as 
	// qc_convert_balls_to_dots has the side-effect of changing last_converted_ball_cdlln.
	
	uc_pr = qc->nest_pr->last_converted_ball_cdlln;
	uc_du = qc->nest_du->last_converted_ball_cdlln;

	// Convert balls to dots and offsets to lines
	qc_convert_balls_to_dots(qc, qc->m_pr, uc_pr, layer->num_blocks_pr);
	qc_convert_balls_to_dots(qc, qc->m_du, uc_du, layer->num_blocks_du);
	qc_convert_offsets_to_lines(qc, qc->m_pr, qc->nest_pr, uc_pr, layer, layer->num_blocks_pr);
	qc_convert_offsets_to_lines(qc, qc->m_du, qc->nest_du, uc_du, layer, layer->num_blocks_du);
}







// ERROR OP FUNCTIONS

/**
 * \defgroup trans Transformations 
 *
 * @{
 */

/**
 * \brief Transforms an \ref error from an identity gate. Does nothing.
 *
 * \param[in] error The \ref error to transform
 */
void qc_iden_transform_error(__attribute__((unused)) ERROR *error) {
}

/**
 * \brief Transforms an \ref error from a loss gate. Does nothing.
 *
 * \param[in] error The \ref error to transform
 */
void qc_loss_transform_error(__attribute__((unused)) ERROR *error) {
}

/**
 * \brief Transforms an \ref error from an H gate. Converts an X error to a Z error and vica versa. 
 *
 * \param[in] error The \ref error to transform
 */
void qc_H_transform_error(ERROR *error) {
	if (error->op == X) {
		error->op = Z;
	}
	else if (error->op == Z) {
		error->op = X;
	}
}

/**
 * \brief Transforms an \ref error from an X gate. Does nothing.
 *
 * \param[in] error The \ref error to transform
 */
void qc_X_transform_error(__attribute__((unused)) ERROR *error) {
}

/**
 * \brief Transforms an \ref error from an Y gate. Does nothing.
 *
 * \param[in] error The \ref error to transform
 */
void qc_Y_transform_error(__attribute__((unused)) ERROR *error) {
}

/**
 * \brief Transforms an \ref error from an Z gate. Does nothing.
 *
 * \param[in] error The \ref error to transform
 */
void qc_Z_transform_error(__attribute__((unused)) ERROR *error) {
}

/**
 * \brief Transforms an \ref error from a S gate. Converts a X error to a Y error and vica versa
 *
 * \param[in] error The \ref error to transform
 */
void qc_S_transform_error(ERROR *error) {
	if (error->op == X) {
		error->op = Y;
	}
	else if (error->op == Y) {
		error->op = X;
	}
}

/**
 * \brief Transforms an \ref error from a T gate. Converts an X error to a Y error 
 *
 * \param[in] error The \ref error to transform
 */
void qc_T_transform_error(ERROR *error) {
	if (error->op == X) {
		error->op = Y;
	}
}

/**
 * \brief Transforms an \ref error from a TX gate (X gate followed by a T gate). Converts a Z error to a Y error 
 *
 * An alias for performing ::qc_X_transform_error then ::qc_T_transform_error.
 *
 * \param[in] error The \ref error to transform
 */
void qc_TX_transform_error(ERROR *error) {
	if (error->op == Z) {
		error->op = Y;
	}
}

/**
 * \brief Transforms an \ref error from a Controlled NOT gate
 *
 * If qubit 1 has an X error (X or Y) then apply an X error to qubit 2
 * If qubit 2 has a Z error (Z or Y) then apply a Z error to qubit 1
 *
 * \param[in] op1 The error operator from qubit 1
 * \param[in] op2 The error operator from qubit 2
 */
void qc_cnot_transform_errors(int *op1, int *op2) {
	if (*op1 & L || *op2 & L) {
		assert(!(*op1 & L && *op2 & L));
		if (*op1 & L) {
			*op2 ^= Y;
		} else if (*op2 & L) {
			*op1 ^= Y;
		}
		return;
	}

	if (*op1 & X) {
		*op2 ^= X;
	}
	if (*op2 & Z) {
		*op1 ^= Z;
	}
}

/**
 * \brief Transforms an \ref error from a Controlled NOT gate followed by a
 * swap gate. 
 *
 * \param[in] op1 The error operator from qubit 1
 * \param[in] op2 The error operator from qubit 2
 */
void qc_cnots_transform_errors(int *op1, int *op2) {
	int op;
	
	if (*op1 & X) {
		*op2 ^= X;
	}
	if (*op2 & Z) {
		*op1 ^= Z;
	}

	op = *op2;
	*op2 = *op1;
	*op1 = op;
}

/**
 * \brief Transforms an \ref error from a controlled-Z gate. 
 *
 * If qubit 1 has an X error (X or Y) then apply a Z error to qubit 2
 * If qubit 2 has an X error (X or Y) then apply a Z error to qubit 1
 *
 * \param[in] op1 The error operator from qubit 1
 * \param[in] op2 The error operator from qubit 2
 */
void qc_cZ_transform_errors(int *op1, int *op2) {
	if (*op1 & X) {
		*op2 ^= Z;
	}
	if (*op2 & X) {   
		*op1 ^= Z;
	}
}

/**
 * \brief Transforms an \ref error from a swap gate.
 *
 * The operator from qubit 1 is swapped with qubit 2
 *
 * \param[in] op1 The error operator from qubit 1
 * \param[in] op2 The error operator from qubit 2
 */
void qc_swap_transform_errors(int *op1, int *op2) {
	int i;

	i = *op1;
	*op1 = *op2;
	*op2 = i;
}

/** @} */





//   
// PRINT FUNCTIONS
//

/**
 * \ingroup qubit
 * \brief Prints a \ref qubit and any \link error Errors\endlink on it
 *
 * \param[in] q The \ref qubit to be printed 
 */
void qc_print_qubit(QUBIT *q) {
	int i;
	CDLL_NODE *n;

	printf("qubit i: %d, j: %d, k: %d, t: %ld, e: %d, num_errors: %d\n", q->i, q->j, q->k, q->t, q->e, q->num_errors);

	n = q->error_cdll->next;
	for (i = 0; i < q->num_errors; i++) {
		printf("- %3d ", i);
		qc_print_error((ERROR *)n->key);
		n = n->next;
	}
}

/**
 * \ingroup error_model
 * \brief Prints an \ref error_model
 *
 * \param[in] em The \ref error_model to print
 * \param[in] p The probability of an error occuring 
 */
void qc_print_em(ERROR_MODEL *em, double p) {
	int i, j;

	printf("nq: %d, duration: %d, scale: %20.16lf\n", em->num_qubits, em->duration, em->scale);

	for (i = 0; i < em->num_lines; i++) {
		printf("%2d %20.16lf %20.16lf %20.16lf", i, em->raw_rel_p[i], em->sum_raw_rel_p[i], p*em->scale*em->raw_rel_p[i]);
		for (j = 1; j <= em->num_qubits; j++) {
			printf(" %ld", em->raw_em[i][j]);
		}
		printf("\n");
	}
}

/**
 * \ingroup error_model
 * \brief Prints all \link error_model Errors Models\endlink in a \ref qc
 *
 * \param[in] qc The \ref qc to print from
 * \param[in] p The probability of an error occuring 
 */
void qc_print_ems(QC *qc, double p) {
	int i;

	for (i=0; i<qc->num_ems; i++) {
		printf("gate: %d, ", i);
		qc_print_em(qc->ems[i], p);
		printf("\n");
	}
}

/**
 * \ingroup error
 * \brief Prints an \ref error
 *
 * \param[in] error The error to print 
 */
void qc_print_error(ERROR *error) {
	printf("big_t: %ld, q: %p, label: %ld, op: %d, p: %g, qc_cdlln: %p, gate: %d\n",
		error->big_t, error->q, error->label, error->op, error->p, error->qc_cdlln, error->gate);
}

/**
 * \ingroup error
 * \brief Prints an \ref error. Accepts a void pointer to an \ref error.
 *
 * \param[in] key The error to print
 */
void qc_print_void_error(void *key) {
	qc_print_error((ERROR *)key);
}

/**
 * \ingroup error
 * \brief Prints a \ref cdll of \ref error
 *
 * \param[in] sent The sentinel node of a \ref cdll containing errors
 */
void qc_print_errors(CDLL_NODE *sent) {
	CDLL_NODE *n;
	int i;

	n = sent->next;
	i = 1;
	while (n != sent) {
		printf("%d: ", i++);
		qc_print_error((ERROR *)n->key);
		n = n->next;
	}
}

/**
 * \ingroup error
 * \brief Prints an \ref error and any related errors in its error circle.
 *
 * \param[in] error The start of the \ref error circle 
 */
void qc_print_error_circle(ERROR *error) {
	ERROR *start;

	start = error;
	printf("start error ");
	qc_print_error(error);
	error = error->next;
	while (error != start) {
		printf("- ");
		qc_print_error(error);
		error = error->next;
	}
}

/**
 * \ingroup error
 * \brief Prints a \ref cdll of \ref error circles
 *
 * \param[in] sent The sentinel node of a \ref cdll containing \ref error
 * circles
 */
void qc_print_error_circles(CDLL_NODE *sent) {
	CDLL_NODE *n;

	n = sent->next;
	while (n != sent) {
		qc_print_error_circle((ERROR *)n->key);
		n = n->next;
	}
}

/**
 * \ingroup set 
 * \brief Prints a \ref set
 *
 * \param[in] set The \ref set to be printed 
 */
void qc_print_set(SET *set) {
	// If there is no set
	if (set == NULL) {
		printf("set == NULL\n");
	}

	// If the set is a boundary
	else if (set->bdy == NULL) {
		printf("type: %d, %p, i: %d, j: %d, big_t: %ld, t: %ld, num_meas_left: %d, bdy: N\n", set->type, set, set->i, set->j, set->big_t, set->t, set->num_meas_left);
	}
	
	else {
		printf("type: %d, %p, i: %d, j: %d, big_t: %ld, t: %ld, num_meas_left: %d, bdy: %d\n", set->type, set, set->i, set->j, set->big_t, set->t, set->num_meas_left, set->bdy->ball->i);
	}

	if (set->ball) {
		printf("ball: mp: %d\n", set->ball->mp);
	}

	if (set->set_pos_cdll) {
		cdll_print(set->set_pos_cdll, qc_print_void_set_pos);
	}
}

/**
 * \ingroup set
 * \brief Prints a \ref set. Accepts a void pointer to a \ref set.
 *
 * \param[in] key The set to be printed 
 */
void qc_print_void_set(void *key) {
	qc_print_set((SET *)key);
}

void qc_print_void_set_pos(void *key) {
	SET_POS *sp;
	
	sp = ((SET_POS *)key);
	printf("SetPos: %d %d %d\n", sp->i, sp->j, sp->lay);
}

/**
 * \ingroup set
 * \brief Prints the \link set Sets\endlink from the \ref qc
 * 
 * \param[in] qc The \ref qc to print the \link set Sets\endlink from 
 */
void qc_print_set_heap(QC *qc) {
	bh_print(qc->set_heap, qc_print_void_set);
}

/**
 * \ingroup de
 * \brief Prints a \ref de. Prints the \ref set and \ref error from the \ref de
 *
 * \param[in] de The \ref de to be printed
 */
void qc_print_de(DE *de) {
	printf("de:\n- ");
	qc_print_set(de->set);
	printf("- ");
	qc_print_error(de->error);
}

/**
 * \ingroup de
 * \brief Prints a \ref de. Prints the \ref set and \ref error from the \ref
 * de. Accepts a void pointer to a \ref de.
 *
 * \param[in] key The \ref de to be printed
 */
void qc_print_void_de(void *key) {
	qc_print_de((DE *)key);
}

/**
 * \ingroup de
 * \brief Prints the \link de Detection Events\endlink from a \ref qc.
 *
 * \param[in] qc The \ref qc from which to print all the \link de Detection Events\endlink. 
 */
void qc_print_des(QC *qc) {
	cdll_print(qc->de_cdll, qc_print_void_de);
}

/**
 * \ingroup ball
 * \brief Prints the destination \link ball Balls\endlink going down a \ref stick from a \ref ball
 *
 * \param[in] ball The source \ref ball from which to find the destination of 
 * \param[in] stick The \ref stick to print the destination \link ball Balls\endlink from.
 */
void qc_print_dest_coord(BALL *ball, STICK *stick) {
	BALL *a, *b;

	a = stick->a;
	b = stick->b;

	if (a == ball) {
		printf("(%d, %d, %ld, %9.7f) ", b->i, b->j, b->big_t, stick->p_stick);
	}
	else {
		printf("(%d, %d, %ld, %9.7f) ", a->i, a->j, a->big_t, stick->p_stick);
	}
}

/**
 * \ingroup ball
 * \brief Prints a \ref ball
 *
 * \param ball The \ref ball to print.
 */
void qc_print_ball(BALL *ball) {
	HT *ht;
	int i;
	DLL_NODE *n;
	STICK *stick;

	printf("%p [%p] (%d, %d, %ld, %ld, %p) (%ld) mp: %d - ", ball, ball->copy, ball->i, ball->j, ball->t, ball->big_t, ball->dot, ball->big_t, ball->mp);

	ht = ball->stick_ht;
	for (i=0; i<ht->length; i++) {
		n = ht->table[i];
		while (n != NULL) {
			stick = (STICK *)n->key;
			qc_print_dest_coord(ball, stick);
			n = n->next;
		}
	}

	printf("\n");
}

/**
 * \ingroup ball
 * \brief Prints a \ref ball. Accepts a void pointer to a \ref ball.
 *
 * \param key The \ref ball to print.
 */
void qc_print_void_ball(void *key) {
	qc_print_ball((BALL *)key);
}

/**
 * \ingroup stick
 * \brief Prints a \ref nest of \link stick Sticks\endlink
 *
 * \param nest The \ref nest to print the \link stick Sticks\endlink 
 * from
 */
void qc_print_sticks(NEST *nest) {
	CDLL_NODE *n;
	STICK *s;

	n = nest->stick_cdll->next;
	while (n != nest->stick_cdll) {
		s = (STICK *)n->key;
		qc_print_stick(s);
		n = n->next;
	}
}

/**
 * \ingroup stick
 * \brief Prints a \ref stick
 *
 * \param[in] stick The \ref stick to print
 */
void qc_print_stick(STICK *stick) {
	BALL *a, *b;
	LL_NODE *n;
	ERROR *e;

	a = stick->a;
	b = stick->b;

	printf("(%d, %d, %ld) -- (%d, %d, %ld) type: %d %d\n",
		a->i, a->j, a->big_t, b->i, b->j, b->big_t,
		a->type, b->type);
	
	n = stick->error_ll;
	while (n != NULL) {
		e = (ERROR *)n->key;
		printf("- l: %ld, op: %d, p: %20.16lf, g: %d, i: %d, j: %d, k: %d, t: %ld\n", e->label, e->op, e->p, e->gate, e->i, e->j, e->k, e->t);
		n = n->next;
	}
}

/**
 * \ingroup stick
 * \brief Prints a \ref stick. Accepts a void pointer to a \ref stick.
 *
 * \param[in] key The \ref stick to print
 */
void qc_print_void_stick(void *key) {
	qc_print_stick((STICK *)key);
}

/**
 * \ingroup nest
 * \brief Prints a \ref nest
 *
 * \param[in] nest The \ref nest to be printed 
 */
void qc_print_nest(NEST *nest) {
	cdll_reverse_print(nest->ball_cdll, qc_print_void_ball);
}

/**
 * \ingroup qc
 * \brief Prints stats about a \ref qc
 *
 * \param qc The \ref qc to print stats about
 */
void qc_print_qc_stats(QC *qc) {
	printf("stats nq%d nsyn%d nset%d ne%d nne%d nde%d Pr nd%d nl%d nv%d nu%d Du nd%d nl%d nv%d nu%d\n",
		qc->num_qubits,
		qc->syn_heap->num_elem,
		qc->set_heap->num_elem,
		qc->num_errors,
		qc->num_nest_errors,
		qc->num_des,
		qc->m_pr->num_dots,
		qc->m_pr->num_lines,
		qc->m_pr->num_vertices,
		qc->m_pr->num_undos,
		qc->m_du->num_dots,
		qc->m_du->num_lines,
		qc->m_du->num_vertices,
		qc->m_du->num_undos
	);
}

/**
 * \ingroup syndrome
 * \brief Prints a \ref syndrome. Accepts a void pointer to a \ref syndrome.
 *
 * \param[in] key The \ref syndrome to print
 */
void qc_print_void_syndrome(void *key) {
	SYNDROME *syn;
	SET *set;

	syn = (SYNDROME *)key;
	set = syn->set;

	printf("qc_heap_i: %d, t: %ld, set: (%d, %d, %ld) nmr: %d\n", syn->qc_heap_i, syn->t, set->i, set->j, set->big_t, set->num_meas_left);
}

/**
 * \ingroup syndrome
 * \brief Prints a \ref heap of \ref syndrome 
 *
 * \param[in] qc The \ref qc containing the \ref heap of \ref syndrome 
 */
void qc_print_syndrome_heap(QC *qc) {
	bh_print(qc->syn_heap, qc_print_void_syndrome);
}

/**
 * \ingroup offset
 * \brief Prints an \ref offset. Accepts a void pointer to an \ref offset.
 *
 * \param k The \ref offset to be printed.
 */
void qc_print_offset(void *k) {
	OFFSET *o;

	o = (OFFSET *)k;
	printf("id = %d, i = %d, j = %d, t = %ld, big_t = %ld, wt = %d, type=%d\n", o->id, o->i, o->j, o->t, o->big_t, o->wt, o->type);
}

/**
 * \ingroup block
 * \brief Prints a \ref block.
 *
 * \param block The \ref block to be printed.
 */
void qc_print_block(BLOCK *block) {
	printf("Block %d,%d,%ld\n", ((BLOCK *)block)->i, ((BLOCK *)block)->j, ((BLOCK *)block)->t);
	ll_print(((BLOCK *)block)->offsets, qc_print_offset);
}

/**
 * \ingroup block
 * \brief Prints a \ref block. Accepts a void pointer to a \ref block.
 *
 * \param block The \ref block to be printed.
 */
void qc_print_void_block(void *block) {
	qc_print_block((BLOCK *)block);	
}

/**
 * \ingroup block
 * \brief Prints a \ref llist of \link block Blocks\endlink
 *
 * \param blocks_head The head node of a \ref llist of \link block Blocks\endlink
 */
void qc_print_blocks(LL_NODE *blocks_head) {
	ll_print(blocks_head, qc_print_void_block);
}

/**
 * \ingroup layer
 * \brief Prints a \ref layer
 *
 * \param[in] layer The \ref layer to be printed
 * \param[in] n The size of the n dimension of the \ref layer 
 * \param[in] m The size of the m dimension of the \ref layer 
 */
void qc_print_layer(LAYER *layer, int n, int m) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			if (layer->blocks[i][j]) {
				if (layer->blocks[i][j]->type == PRIMAL) {
					printf("%c[0;31m%3d", 0x1B, layer->blocks[i][j]->id);
				} else {
					printf("%c[0;34m%3d", 0x1B, layer->blocks[i][j]->id);
				}
				printf("%c[0m", 0x1B);
			} else {
				printf("   ");
			}
		}
		printf("\n");
	}
}

/**
 * \ingroup recipe
 * \brief Prints a \ref recipe
 *
 * \param[in] recipe The recipe to be printed 
 */
void qc_print_recipe(RECIPE *recipe) {
	int i, p, d, pb, db;
	LL_NODE *node;
	OFFSET *o;

	p = d = pb = db = 0;

	printf("Printing recipe:\n");
	printf("Printing layers, repeat id %d:\n", recipe->repeat_id);
	for (i = 0; i < recipe->size; i++) {
		printf("Layer: %d\n", i);
		qc_print_layer(recipe->layers[i], recipe->n, recipe->m);
	}

	printf("\n");
	printf("Printing blocks:\n");
	node = recipe->blocks;
	while (node != NULL) {
		qc_print_block((BLOCK *)node->key);
		node = node->next;
	}

	printf("\n");
	printf("Printing offsets:\n");
	node = recipe->offsets;
	while (node != NULL) {
		o = (OFFSET *)node->key;
		if (o->type == PRIMAL) ++p;
		if (o->type == DUAL) ++d;
		if (o->type == PRIMAL_BOUNDARY) ++pb;
		if (o->type == DUAL_BOUNDARY) ++db;

		qc_print_offset(o);
		node = node->next;
	}

	printf("p %d d %d pb %d db %d\n", p, d, pb, db);
}

/**
 * \ingroup recipe_adv
 * \brief Prints a \ref recipe_adv
 *
 * \param[in] recipe The advanced recipe to be printed 
 */
void qc_print_recipe_adv(RECIPE_ADV *recipe) {
	CDLL_NODE *node;
	LL_NODE *lln;
	BLOCK *block;

	printf("Type: %d\n", recipe->type);
	printf("Length: %d\n", recipe->cycle_len);
	printf("In Cycle? %d\n", recipe->in_cycle);
	printf("Cycle t0: %ld\n", recipe->cycle_t0);
	printf("Cycle Period: %ld\n", recipe->cycle_period);

	node = recipe->block_cdll->next;
	while (node != recipe->block_cdll) {
		block = (BLOCK *)node->key;
		if (node == recipe->cycle_begin) {
			printf("** CYCLE START **\n");
		}

		printf("BLOCK (type: %d) (%d, %d, %ld, %ld)\n", block->type, block->i, block->j, block->t, block->big_t);
		lln = block->offsets;
		while (lln != NULL) {
			qc_print_offset(lln->key);
			lln = lln->next;
		}
		
		if (node == recipe->cycle_end) {
			printf("** CYCLE END **\n");
			break;
		}
		node = node->next;
	}
}

