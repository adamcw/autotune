#ifndef QC_H
#define QC_H

#include "qc_td.h"
#include "../match/match_td.h"

#include <stdio.h>
#include <stdbool.h>
#include "../llist/llist.h"
#include "../cdllist/cdllist.h"
#include "../hash_table/hash_table.h"
#include "../bheap/bheap.h"
#include "../match/match.h"

// OPERATORS
#define I 0
#define X 1
#define Z 2
#define Y 3
#define L 4

// PRIMAL/DUAL DEFINITIONS
#define PRIMAL_BOUNDARY 2
#define PRIMAL 1
#define UNKNOWN 0
#define DUAL -1
#define DUAL_BOUNDARY -2

// QC_UNDO DEFINITIONS
#define INCREMENT_BIG_T 0
#define SINGLE_QUBIT_GATE 1
#define TWO_QUBIT_GATE 2

// MISCELLANEOUS
#define STICK_HT_SIZE 101
#define TRUE 1
#define FALSE 0
#define FAIL_LABEL -1 // set to value of problematic error label (debugging only)
#define DONE 1

// RECIPE TYPE DEFINITIONS
#define RECIPE_INFINITE 1
#define RECIPE_FINITE 2

// TIMER DEFINITIONS
#define START 1
#define STOP 0
#define TIMERS 50
extern double t_start[TIMERS];
extern double t_total[TIMERS];
extern int t_counts[TIMERS];
extern const char *t_name[TIMERS];
void timer(int timer_num, int action, const char *name);

/** \defgroup qc Quantum Computer */
/** \defgroup qubit Qubit */
/** \defgroup error_model Error Model */
/** \defgroup error Error */
/** \defgroup de Detection Event */
/** \defgroup syndrome Syndrome Qubit */
/** \defgroup set Set */
/** \defgroup measurement Measurement */
/** \defgroup nest Nest */
/** \defgroup ball Ball */
/** \defgroup stick Stick */
/** \defgroup recipe Recipe */
/** \defgroup layer Layer */
/** \defgroup block Block */
/** \defgroup offset Offset */

typedef struct qc QC; 
typedef struct qubit QUBIT;
typedef struct error_model ERROR_MODEL;
typedef struct error ERROR;
typedef struct de DE; 
typedef struct syndrome SYNDROME;
typedef struct set SET;
typedef struct set_pos SET_POS;
typedef struct measurement MEASUREMENT;
typedef struct nest NEST;
typedef struct ball BALL;
typedef struct stick STICK;
typedef struct recipe RECIPE;
typedef struct layer LAYER;
typedef struct block BLOCK;
typedef struct offset OFFSET;
typedef struct recipe_adv RECIPE_ADV;

/**
 * \ingroup qc
 * \brief The primary quantum computer structure
 */
struct qc { 
	/** The first seed to the random number generator */
	int s0;

	/** The second see to the random number generator */
	int s1;

	/** The number of complete rounds of stabilizer measurement */
	long int big_t; 
	
	/** The big_t at which there are or will be unfinalized balls */
	long int unfinalized_big_t; 

	/** The fixed maximum number of qubits used during the computation */
	int num_qubits; 
	
	/** A circular double-linked list with the qubits as keys */
	CDLL_NODE *qubit_cdll;

	/** A heap containing all of the syndrome qubits */
	BHEAP *syn_heap;

	/** A heap containing all of the sets */
	BHEAP *set_heap;

	CDLL_NODE *finalized_balls;
	CDLL_NODE *zombie_dots;

	/** The number of error models in the given architecture */
	int num_ems;

	/** An array of error models */
	ERROR_MODEL **ems;

	/** Cumulative count of the number of random errors that have occured during simulation */
	long int num_es;
	
	/** Next Error ID: The next integer to assign to an error to uniquely identify it */
	long int next_label; 
	
	/** The number of newest errors that have not yet been added to the nest (haven't fully propagated) */
	int num_errors;

	/** A circule double-linked list of the newest errors in the computation */
	CDLL_NODE *error_cdll;

	/** The number of errors that have occured that are in the nest */
	int num_nest_errors;
	
	/** A circular double-linked list of the errors in the nest. Allows for
	 * their deletion from qubits. 
	 */
	CDLL_NODE *nest_error_cdll; 
	
	/** The number of detection events */
	int num_des;
	
	/** A hash table of detection events */
	HT *de_ht;

	/** A circular double-linked list of detection events */
	CDLL_NODE *de_cdll;

	/** The primal nest */
	NEST *nest_pr; 

	/** The dual nest */
	NEST *nest_du; 

	/** The primal matching problem */
	MATCHING *m_pr;

	/** The dual matching problem */
	MATCHING *m_du;

	/** The recipe for the quantum computation */
	RECIPE *recipe;
	RECIPE_ADV *recipe_adv;
	
	/** Whether or not gates can be applied without error. */
	int perfect_gates;

	/** Whether or not to track errors. **/
	int track;

	/** A function that takes a given position and type and calculates the
	 * nearest boundary and returns its BALL. 
	 */
	BALL *(*get_boundary)(int i, int j, long int big_t, int type, void *boundaries);

	/** An array of boundaries */
	void *boundaries;

	int logical_loss_err_pr;
	int last_logical_loss_err_pr;
	int logical_loss_err_du;
	int last_logical_loss_err_du;
};

/**
 * \ingroup qubit
 * \brief A qubit
 */
struct qubit {
	/** \brief The time position of the qubit */
	long int t;

	/** \brief The i coordinate of the qubit */
	int i;
	
	/** \brief The j coordinate of the qubit */
	int j;
	
	/** \brief The k coordinate of the qubit */
	int k; 
	
	/** \brief Used for stochastic errors, value represents the error on the qubit
	 *  - 0 = Identity 
	 *  - 1 = X Error
	 *  - 2 = Y Error
	 *  - 3 = Z Error
	 */
	int e; 

	/** \brief The number of errors in the qubit */
	int num_errors;

	/** \brief A circular double-linked list of all the errors in the qubit */
	CDLL_NODE *error_cdll;

	/** \brief A hash-table of the errors in the qubit */
	HT *error_ht;

	/** \brief \ref qc contains a \ref cdll of all \link qubit Qubits\endlink. This is the
	 * ::CDLL_NODE in that list, that this qubit is in. 
	 */
	CDLL_NODE *qc_cdlln;
	
	/** \brief A stored copy of the qubit, used to restore a qubit to a previous state */
	QUBIT *copy;
};

/**
 * \ingroup error_model
 * \brief A user defined error model for a 1 or 2 qubit gate. Additional
 * information can be found in Appendix B of the paper Topological Code
 * Autotune, 2012.
 */
struct error_model {
	/** The number of qubits the error model effects. */
	int num_qubits;

	/** The value that the relative strengths will be normalised to sum to */
	double scale;

	/** The number of possible errors in the error model */
	int num_lines;

	/** A 2D array of the raw error model. Each row is a possible error, with
	 * the first column being the relative strength of the error, and the
	 * remaining columns the error to be applied to each qubit.
	 */
	long int **raw_em;

	/** The raw probabilties as converted from the relative strengths provided
	 * in the error model file. 
	 *
	 * \f[  
	 *	  raw\_rel\_p_j = \frac{
	 *		  relative\_strength_j
	 *	  }
	 *	  {
	 *		  \sum\limits_{i=0}^{num\_lines-1}{
	 *			  relative\_strength_i
	 *		  }
	 *	  }
	 * \f]
	 */ 
	double *raw_rel_p;

	/** The raw probabilties as converted from the relative strengths provided
	 * in the error model file. 
	 *
	 * \f[  
	 *	  sum\_raw\_rel\_p_j = \sum\limits_{i=0}^{j}{
	 *		  probability_i
	 *	  }
	 * \f]
	 */	 
	double *sum_raw_rel_p;

	/** The duration of the error in arbitrary units */
	int duration;
};

/**
 * \ingroup error
 * \brief An error that has occured due to the use of a gate
 */
struct error {
	/** i coordinates of the error's creation */
	int i;
	
	/** j coordinates of the error's creation */
	int j;
	
	/** k coordinates of the error's creation */
	int k; 

	/** Time of the error's creation */
	long int t; 
	
	/** big_t of the error's creation */
	long int big_t;

	/** The unique ID for this error */
	long int label;
	
	/** The type of the error. PRIMAL, DUAL or UNKNOWN */
	int type; 

	/** The Pauli operator of the error. I, X, Y, Z */
	int op; 

	/** The gate that created the error */
	int gate;

	/** The probability of the error occuring */
	double p;

	/** \brief \ref qc contains a \ref cdll of all \link error Errors\endlink. This is the
	 * ::CDLL_NODE in that list that this \ref error is in. 
	 */
	CDLL_NODE *qc_cdlln;
	
	/** \brief A \ref qubit contains a \ref cdll of all \link error Errors\endlink on it. This is the
	 * ::CDLL_NODE in that list that this \ref error is in. 
	 */
	CDLL_NODE *q_cdlln;
	
	/** The qubit on which the error was created */
	QUBIT *q;

	/** The stick that the error belongs to */
	STICK *stick;

	/** The next error in a combined error event (multiple errors from a single gate */
	ERROR *next;

	/** The previous error in a combined error event (multiple errors from a single gate */
	ERROR *prev; 	
	
	/** A stored copy of the error, used to restore the error to a previous state */
	ERROR *copy;
};

/**
 * \ingroup de
 * \brief A detection event
 */
struct de {
	/** The \ref set the \ref de belongs to */
	SET *set;

	/** The \ref error that was detected */
	ERROR *error;

	/** The label of the \ref error that was detected */
	long int label;

	/** \brief A \ref qc contains a \ref cdll of all \ref de%s in it. This is the
	 * ::CDLL_NODE in that list that this \ref de is in. 
	 */
	CDLL_NODE *qc_cdlln;

	/** A stored copy of the \ref de, used to restore the \ref de to a previous state */
	DE *copy;
};

/**
 * \ingroup syndrome
 * \brief A syndrome qubit
 */
struct syndrome {
	/** The t that the \ref syndrome was last measured */
	long int t;

	/** The index of this \ref syndrome in the \ref syndrome \ref heap in \ref qc */
	int qc_heap_i;

	/** The \ref set that this \ref syndrome belongs to */
	SET *set;

	/** \brief A \ref set contains a \ref cdll of all \ref syndrome%s in it. This is the
	 * ::CDLL_NODE in that list that this \ref syndrome is in. 
	 */
	CDLL_NODE *set_cdlln;

	/** A stored copy of the \ref syndrome, used to restore the \ref syndrome to a previous state */
	SYNDROME *copy;
};

/**
 * \ingroup set
 * \brief A set
 */
struct set {
	/** The type of the \ref set (PRIMAL, DUAL, PRIMAL_BOUNDARY or
	 * DUAL_BOUNDARY) */
	int type;

	/** The i coordinate */
	int i;

	/** The j coordinate */
	int j;

	/** The big_t of the \ref set */
	long int big_t;

	/** The t of the \ref set */
	long int t;

	/** The index of this \ref set in the \ref set \ref heap in \ref qc */
	int qc_heap_i;

	/** The number of \ref measurement%s left for this \ref set to be complete */
	int num_meas_left;

	/** The boundary this \ref set is closest to */
	SET *bdy;

	BLOCK *block;

	/** */
	CDLL_NODE *set_pos_cdll;

	/** A \ref cdll of \ref syndrome%s for this \ref set */
	CDLL_NODE *syn_cdll;

	/** The number of \ref error%s in this \ref set */
	int num_errors;

	/** A \ref cdlln of \ref measurement%s on this \ref set */
	CDLL_NODE *mt_cdll;

	/** The \ref ball associated with this \ref set */
	BALL *ball;

	/** A stored copy of the \ref set, used to restore the \ref set to a previous state */
	SET *copy;
};

/**
 * \ingroup set
 * \brief A set's position
 */
struct set_pos {
    /** The i coordinate */
	int i;

    /** The j coordinate */
    int j;

    /** The layer of the \ref set */
	int lay;

	/** The \ref set of this \ref set_pos */
	SET *set;
};

/**
 * \ingroup measurement
 * \brief A measurement
 */
struct measurement {
	/** The stochastic measurement result (-1, 0, 1) */
	int m; 

	/** The number of \ref set%s that this \ref measurement is in */
	int num_parent_sets;

	/** The number of \ref error%s in the \ref measurement */
	int num_errors;

	/** A \ref cdll of \ref error%s in the \ref measurement */
	CDLL_NODE *error_cdll;

	/** A stored copy of the \ref measurement, used to restore the \ref measurement to a previous state */
	MEASUREMENT *copy;
};

/**
 * \ingroup nest
 * \brief A nest 
 */
struct nest {
	/** A \ref cdll of all the \ref ball%s in the \ref nest */
	CDLL_NODE *ball_cdll;

	/** A ::CDLL_NODE of the last \ref ball to be turned into a \ref block */
	CDLL_NODE *last_blocked_ball_cdlln;

	/** A ::CDLL_NODE of the last \ref ball converted into a \ref dot */
	CDLL_NODE *last_converted_ball_cdlln;

	/** A \ref cdll of all the \ref stick%s in the \ref nest */
	CDLL_NODE *stick_cdll;

	/** A \ref ht of all the \ref stick%s in the \ref nest */ 
	HT *stick_ht;

	/** A stored copy of the \ref nest, used to restore the \ref nest to a previous state */
	NEST *copy;
};

/**
 * \ingroup ball
 * \brief A ball 
 */
struct ball {
	/** The type of ball (PRIMAL/DUAL/PRIMAL_BOUNDARY/DUAL_BOUNDARY) */
	int type;

	/** The i coordinate */
	int i;

	/** The j coordinate */
	int j;

	/** The big_t */
	long int big_t;

	/** The t **/
	long int t;

	int hash_ijt;

	/** A \ref ht of \ref stick%s connected to this \ref ball */
	HT *stick_ht;

	/** The \ref nest contains a \ref cdll of \ref ball%s. This is the
	 * ::CDLL_NODE in that list that corresponds to this \ref ball. 
	 */
	CDLL_NODE *nest_cdlln;

	/** Measurement Product */
	int mp; 

	/** The \ref dot associated with this \ref ball */
	DOT *dot;

	/** A copy of this \ref ball */
	BALL *copy;
};

/**
 * \ingroup stick
 * \brief A stick
 */
struct stick {
	/** The probability of the \ref stick */
	double p_stick;

	/** The \ref ball at the "start" of the \ref stick */
	BALL *a;

	/** The \ref ball at the "end" of the \ref stick */
	BALL *b;

	/** The number of \ref error%s on the \ref stick */
	int num_errors;

	/** A \ref llist of \ref error%s on the \ref stick */
	LL_NODE *error_ll;

	/** The \ref nest contains a \ref cdll of \ref stick%s. This is the
	 * ::CDLL_NODE in that list that corresponds to this \ref stick. 
	 */
	CDLL_NODE *nest_cdlln;

	/** A copy of this \ref stick */
	STICK *copy;
};

/**
 * \ingroup recipe
 * \brief A recipe consisting of the layer, block and offset information
 * require to perform a quantum computation.
 */
struct recipe {
	/** The type of \ref recipe (currently only supports INFINITE) */
	int type;

	/** The size of the \ref recipe (number of possible \ref layer%s) */
	int size;

	/** The size of the n dimension */
	int n;

	/** The size of the m dimension */
	int m;

	/** The minimum weight horizontal edge connecting the primal boundaries */
	int min_horz_wt_pr;
	
	/** The minimum weight horizontal edge connecting the dual boundaries */
	int min_horz_wt_du;

	/** The number of \ref layer%s in the \ref recipe */
	int num_layers;
	
	/** An array of \ref layer%s */
	LAYER **layers;

	/** A \ref llist of \ref block%s used */
	LL_NODE *blocks;

	/** A \ref llist of \ref offset%s used */
	LL_NODE *offsets;
	
	/** First lattice of \ref dot%s */
	DOT ***dotarr1;
	
	/** Second lattice of \ref dot%s */
	DOT ***dotarr2;

	// Infinite: The id of the repeated \ref layer
	int repeat_id;

	/** Infinite: The count of repeated \ref layer%s */
	long int repeated_count;
};

struct recipe_adv {
	int type;
	int n;
	int m;

	/** The minimum weight horizontal edge connecting the primal boundaries */
	int min_horz_wt_pr;
	
	/** The minimum weight horizontal edge connecting the dual boundaries */
	int min_horz_wt_du;

	/** Hash table of unique blocks */
	HT *block_ht;
	/** Hash table of unique offsets */
	HT *offset_ht;

	CDLL_NODE *block_cdll;
	CDLL_NODE ***block_arr;
	long int **t_arr;

	int in_cycle;
	int cycle_len;
	long int cycle_t0;
	long int cycle_period;
	CDLL_NODE *cycle_begin;
	CDLL_NODE *cycle_end;
};

/**
 * \brief A layer of blocks 
 */
struct layer {
	/** The number of primal \ref block%s */
	int num_blocks_pr;

	/** The number of dual \ref block%s */
	int num_blocks_du;

	/** An 2D array of \ref block%s */
	BLOCK ***blocks;
};

/**
 * \brief A block of offsets
 */
struct block {
	/** The id of the \ref block */
	int id;

	/** The type of \ref block (PRIMAL/DUAL) */
	int type;

	/** The position of the block in the i dimension */
	int i;

	/** The position of the block in the j dimension */
	int j;

	long int t;

	long int big_t;

	/** The hash of the block based on its location */
	int hash_ijt;

	/** A \ref llist of \ref offset%s for this \ref block */
	LL_NODE *offsets;
};

/**
 * \brief An offset 
 */
struct offset {
	/** The unique ID of the offset */
	int id;

	/** The offset in the i dimension */
	int i;

	/** The offset in the j dimension */
	int j;

	/** The offset in the big_t dimension */
	long int big_t;

	/** The offset in the t dimension */
	long int t;

	/** The type of \ref offset */
	int type;

	/** The weight of the \ref offset */
	int wt;

	double p;

	/** The hash of this offset */
	int hash_ijt;
};


// QC FUNCTIONS
QC *qc_create_qc(int s0, int s1, int de_ht_size, int stick_ht_size, int t_delete, RECIPE *recipe);
QC *qc_create_qc_adv(int s0, int s1, int de_ht_size, int stick_ht_size, int t_delete, RECIPE_ADV *recipe);
void qc_init_qc(QC *qc, int s0, int s1, int de_ht_size, int stick_ht_size, int t_delete);
void qc_free_uninserted_dots_and_lines(NEST *nest);
void qc_free_unfinalized_dots_and_lines(CDLL_NODE *finalized_balls); 
void qc_free_qc(QC *qc);
void qc_free_qc_copy(QC *qc);
void qc_increment_big_t(QC *qc);
void qc_mwpm(QC *qc, int undo);
void qc_undo_mwpm(QC *qc);
QC *qc_copy_qc(QC *qc);
QUBIT *qc_copy_qubit(QUBIT *q);
void *qc_copy_void_qubit(void *key);
ERROR *qc_copy_error(ERROR *error);
void *qc_copy_void_error(void *key);
STICK *qc_copy_stick(STICK *stick);
void *qc_copy_void_stick(void *key);
BALL *qc_copy_ball(BALL *ball);
void *qc_copy_null_ball(void *key); 
void *qc_copy_null_dot(void *key); 
void *qc_copy_void_ball(void *key);
SYNDROME *qc_copy_syndrome(SYNDROME *syn);
void *qc_copy_void_syndrome(void *key);
SET *qc_copy_set(SET *set);
void *qc_copy_void_set(void *key);
SET_POS *qc_copy_set_pos(SET_POS *sp);
void *qc_copy_void_set_pos(void *key);
CDLL_NODE *qc_copy_syn_cdll(CDLL_NODE *syn_cdll);
MEASUREMENT *qc_copy_measurement(MEASUREMENT *mt);
void *qc_copy_void_measurement(void *key);
DE *qc_copy_de(DE *de);
void *qc_copy_void_de(void *key);
NEST *qc_copy_nest(NEST *nest);

// QUBIT FUNCTIONS
QUBIT *qc_create_qubit(int i, int j, int k, int ht_size);
void qc_free_qubit(QUBIT *q);
void qc_free_void_qubit(void *key);
void qc_insert_qubit(QC *qc, QUBIT *q);
void qc_qubit_set_qc_cdlln(void *key, CDLL_NODE *n);
QUBIT *qc_create_and_insert_qubit(QC *qc, int i, int j, int k, int ht_size);
void qc_delete_qubit(QC *qc, QUBIT *q);
void qc_init_qubit(QC *qc, QUBIT *q, int gate, double p,
	int (*multiply_es)(int e1, int e2));
void qc_init_qubit_invisible(QC *qc, QUBIT *q, int gate, double p,
	int (*multiply_es)(int e1, int e2));
void qc_single_qubit_unitary(QC *qc, QUBIT *q, int gate, double p,
	void (*transform_e)(QUBIT *q),
	int (*multiply_es)(int e1, int e2),
	void (*transform_error)(ERROR *error));
void qc_two_qubit_unitary(QC *qc, QUBIT *q1, QUBIT *q2, int gate, double p,
	void (*transform_es)(QUBIT *q1, QUBIT *q2),
	int (*multiply_es)(int e1, int e2),
	void (*transform_errors)(int *op1, int *op2));
int qc_meas_qubit(QC *qc, QUBIT *q, int gate, int basis, double p,
	SET *set1, SET *set2,
	int (*multiply_es)(int e1, int e2),
	int (*meas_sim)(QUBIT *q));

// ERROR MODEL FUNCTIONS
ERROR_MODEL *qc_create_error_model(const char *fname);
void qc_free_error_model(ERROR_MODEL *em);
int qc_insert_error_model(QC *qc, ERROR_MODEL *em);
void qc_random_raw_e(QC *qc, QUBIT *q, ERROR_MODEL *em, double p,
	int (*multiply_es)(int e1, int e2));
void qc_random_raw_e2(QC *qc, QUBIT *q1, QUBIT *q2, ERROR_MODEL *em, double p,
	int (*multiply_es)(int e1, int e2));
void qc_add_errors(QC *qc, QUBIT *q, int gate, double p);
void qc_add_errors2(QC *qc, QUBIT *q1, QUBIT *q2, int gate, double p);

// ERROR FUNCTIONS
ERROR *qc_create_error(int i, int j, int k, long int t, long int big_t, int label, int op, int gate, double p);
void qc_free_error(ERROR *error);
void qc_free_void_error(void *key);
void qc_null_error_q(void *key);
void qc_error_set_qc_cdlln(void *key, CDLL_NODE *n);
void qc_error_set_q_cdlln(void *key, CDLL_NODE *n);
void qc_insert_error(QC *qc, QUBIT *q, ERROR *error, ERROR *sentinel);
void qc_insert_error_nearby(QC *qc, QUBIT *q, ERROR *error, ERROR *sentinel);
void qc_mod_insert_error_in_ht(HT *ht, ERROR *error);
void qc_delete_error(QC *qc, ERROR *error);
int qc_error_eq_label(void *key1, void *key2);
void qc_move_error_to_nest(QC *qc, ERROR *error);

// DE FUNCTIONS
DE *qc_create_de(SET *set, ERROR *error);
void qc_free_de(DE *de);
void qc_free_void_de(void *key);
void qc_de_set_qc_cdlln(void *key, CDLL_NODE *n);
void qc_insert_de(QC *qc, DE *de);
void qc_delete_de(QC *qc, DE *de);
int qc_de_eq_label(void *key1, void *key2);
DE *qc_find_de(QC *qc, int label, int type);
void qc_finalize_de(QC *qc, DE *de);

// SYNDROME FUNCTIONS
SYNDROME *qc_create_syndrome();
void qc_free_syndrome(SYNDROME *syn);
void qc_insert_syndrome(QC *qc, SYNDROME *syn);
void qc_uninsert_syndrome(QC *qc, SYNDROME *syn);
int qc_syndrome_lt_t(void *key1, void *key2);
void qc_syndrome_swap(BHEAP *h, int i1, int i2);
void qc_syndrome_set_qc_heap_i(void *key, int i);
void qc_associate_syndrome(SET *set, SYNDROME *syn);
void qc_unassociate_syndrome(SYNDROME *syn);
void qc_syndrome_set_set_cdlln(void *key, CDLL_NODE *n);

// SET FUNCTIONS
SET *qc_create_set(int type, int i, int j, int num_meas_left, SET *bdy);
SET *qc_create_set_adv(QC *qc, int type, int i, int j, int num_meas_left, SET *bdy); 
SET_POS *qc_create_set_pos(int i, int j, int lay);
void qc_swap_set_layer(SET *set, int lay);
void qc_free_bdy(SET *set);
void qc_free_set(SET *set);
void qc_free_void_set(void *key);
void qc_free_void_set_and_dot(void *key); 
void qc_insert_set(QC *qc, SET *set);
void qc_uninsert_set(QC *qc, SET *set);
int qc_set_lt_big_t(void *key1, void *key2);
void qc_set_swap(BHEAP *h, int i1, int i2);
void qc_set_set_qc_heap_i(void *key, int i);
int qc_is_boundary_set(SET *set);
void qc_merge_dots_without_boundary(MATCHING *m, DOT *dot1, DOT *dot2, LL_NODE *lln_save); 
int qc_merge_dots_with_boundary(MATCHING *m, DOT *dot1, DOT *dot2, LL_NODE *lln_save); 
SET *qc_merge_sets(QC *qc, SET *set1, SET *set2);
void qc_finalize_set(QC *qc, SET *set, long int t);

// MEASUREMENT FUNCTIONS
MEASUREMENT *qc_create_measurement(int m);
void qc_free_measurement(MEASUREMENT *mt);
void qc_free_void_measurement(void *key);
void qc_insert_measurement(SET *set, MEASUREMENT *mt);

// NEST FUNCTIONS
NEST *qc_create_nest(int stick_ht_size);
void qc_free_nest(NEST *nest);
void qc_finalize_local_nest(BALL *ball);
void qc_finalize_nests(QC *qc, long int big_t);
void qc_trim_nests(QC *qc, long int big_t);
void qc_convert_nests(QC *qc, int undo);
void qc_convert_nest(QC *qc, NEST *nest, MATCHING *m, int undo);

// BALL FUNCTIONS
BALL *qc_create_ball(SET *set, int stick_ht_size);
BALL *qc_create_ball_adv(SET *set, int stick_ht_size);

BALL *qc_create_ball_raw(int type, int i, int j, long int big_t, int stick_ht_size);
BALL *qc_create_ball_raw_adv(int type, int i, int j, long int t, long int big_t, int stick_ht_size);

void qc_free_ball(BALL *ball);
void qc_free_void_ball(void *key);
void qc_insert_ball(NEST *nest, BALL *ball);
void qc_delete_ball(NEST *nest, BALL *ball);
void qc_ball_set_nest_cdlln(void *key, CDLL_NODE *n);

// STICK FUNCTIONS
STICK *qc_create_stick(BALL *a, BALL *b);
void qc_free_stick(STICK *stick);
void qc_free_void_stick(void *key);
STICK *qc_find_stick(QC *qc, BALL *b1, BALL *b2);
void qc_stick_set_nest_cdlln(void *key, CDLL_NODE *n);
void qc_insert_stick(NEST *nest, STICK *stick);
void qc_delete_stick(NEST *nest, STICK *stick);
int qc_stick_eq_balls(void *key1, void *key2);
void qc_finalize_stick(STICK *stick);


// RECIPE FUNCTIONS
RECIPE *qc_create_recipe(int type, int num, int n, int m); 
RECIPE_ADV *qc_create_recipe_adv(int type, int n, int m, int copy);
void qc_free_recipe(RECIPE *recipe);
void qc_free_recipe_adv(RECIPE_ADV *recipe);
void qc_reset_recipe_adv(RECIPE_ADV *recipe); 

// LAYER FUNCTIONS
LAYER *qc_create_layer(int n, int m);
void qc_add_layer(RECIPE *recipe, LAYER *layer); 
int qc_compare_layers(LAYER *l1, LAYER *l2, int n, int m);
int qc_is_empty_layer(LAYER *layer, int n, int m);

LAYER *qc_build_layer(QC *qc, RECIPE *recipe, int n, int m, long int layer_big_t);

void qc_free_layer(LAYER *layer, int n);

// BLOCK FUNCTIONS
BLOCK *qc_create_block(RECIPE *recipe, BALL *ball, int type);
BLOCK *qc_create_block_adv(RECIPE_ADV *recipe, BALL *ball);

void qc_create_and_insert_blocks(QC *qc, RECIPE *recipe, LAYER *layer, long int layer_big_t, int type);
int qc_create_and_insert_blocks_adv(QC *qc, RECIPE_ADV *recipe, int type, long int thresh_big_t, int thresh_cycle_len);

LL_NODE *qc_add_block(LL_NODE *blocks_head, BLOCK *block);
BLOCK *qc_block_lookup(LL_NODE *blocks_head, BLOCK *block);
int qc_block_eq(BLOCK *b1, BLOCK *b2);
int qc_block_eq_adv(BLOCK *b1, BLOCK *b2);
void qc_free_block(BLOCK *block);
void qc_free_void_block(void *block);

// OFFSET FUNCTIONS
OFFSET *qc_create_offset(int i, int j, long int big_t, double wt, int type);
OFFSET *qc_create_offset_adv(int i, int j, long int t, double wt, int type);

LL_NODE *qc_add_offset(LL_NODE *offset_head, OFFSET *offset);
OFFSET *qc_find_offset(BLOCK *block, int i, int j, long int big_t, double wt);
OFFSET *qc_offset_lookup(LL_NODE *offset_head, OFFSET *offset);

OFFSET *qc_offset_lookup_adv(HT *ht, int hash, OFFSET *offset);

int qc_offset_eq(void *k1, void *k2);
int qc_offset_eq_adv(void *k1, void *k2);
int qc_offset_lt(void *k1, void *k2);
int qc_offset_lt_adv(void *k1, void *k2);
void qc_free_offset(OFFSET *offset);

// BOOTUP FUNCTIONS
int qc_boot_up(QC *qc, RECIPE *recipe, int n, int m, long int switch_time);
int qc_boot_up_adv(QC *qc, RECIPE_ADV *recipe, long int thresh_big_t, int thresh_cycle_len);
int qc_boot_up_infinite(QC *qc, RECIPE *recipe, int n, int m, long int switch_time);

// CONVERT FUNCTIONS
int qc_hash_ball(BALL *ball); 
int qc_hash_ijt(int i, int j, long int t);
BLOCK *qc_get_block_for_set(QC *qc, SET *set); 
void qc_convert_block_to_dot_and_lines(QC *qc, MATCHING *matching, BALL *ball, BLOCK *block);
void qc_convert_balls_to_dots(QC *qc, MATCHING *matching, CDLL_NODE *ball_cdll, int num_blocks);
void qc_convert_offsets_to_lines(QC *qc, MATCHING *matching, NEST *nest, CDLL_NODE *ball_cdll, LAYER *layer, int num_blocks);

// CONVERT INFINITE FUNCTIONS

LAYER *qc_get_next_layer_infinite(QC *qc);
void qc_convert_nest_infinite(QC *qc, int undo);

// ERROR OP FUNCTIONS
void qc_iden_transform_error(ERROR *error);
void qc_loss_transform_error(__attribute__((unused)) ERROR *error);
void qc_H_transform_error(ERROR *error);
void qc_X_transform_error(ERROR *error);
void qc_Y_transform_error(ERROR *error);
void qc_Z_transform_error(ERROR *error);
void qc_S_transform_error(ERROR *error);
void qc_T_transform_error(ERROR *error);
void qc_TX_transform_error(ERROR *error);
void qc_cnot_transform_errors(int *op1, int *op2);
void qc_cnots_transform_errors(int *op1, int *op2);
void qc_cZ_transform_errors(int *op1, int *op2);
void qc_swap_transform_errors(int *op1, int *op2);

// PRINT FUNCTIONS
void qc_print_qubit(QUBIT *q);
void qc_print_em(ERROR_MODEL *em, double p);
void qc_print_ems(QC *qc, double p);
void qc_print_error(ERROR *error);
void qc_print_void_error(void *key);
void qc_print_errors(CDLL_NODE *sent);
void qc_print_error_circle(ERROR *error);
void qc_print_error_circles(CDLL_NODE *sent);
void qc_print_set(SET *set);
void qc_print_void_set_pos(void *key); 
void qc_print_void_set(void *key);
void qc_print_set_heap(QC *qc);
void qc_print_de(DE *de);
void qc_print_void_de(void *key);
void qc_print_des(QC *qc);
void qc_print_dest_coord(BALL *ball, STICK *stick);
void qc_print_ball(BALL *ball);
void qc_print_void_ball(void *key);
void qc_print_sticks(NEST *nest);
void qc_print_stick(STICK *stick);
void qc_print_void_stick(void *key);
void qc_print_nest(NEST *nest);
void qc_print_qc_stats(QC *qc);
void qc_print_void_syndrome(void *key);
void qc_print_syndrome_heap(QC *qc);

void qc_print_offset(void *k);
void qc_print_block(BLOCK *block);
void qc_print_void_block(void *block);
void qc_print_blocks(LL_NODE *blocks_head);
void qc_print_layer(LAYER *layer, int n, int m);
void qc_print_recipe(RECIPE *recipe);
void qc_print_recipe_adv(RECIPE_ADV *recipe);

#endif
