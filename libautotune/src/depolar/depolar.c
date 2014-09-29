#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "depolar.h"
#include "../memory/memory.h"

// #define DEBUG

// DP_QC FUNCTIONS

/**
 * \brief Creates a new \ref dp_qc
 *
 * \param[in] s0 The first seed of the random number generator
 * \param[in] s1 The second seed of the random number generator
 * \param[in] de_ht_size The size of the \ref de hash table
 * \param[in] stick_ht_size The size of the \ref stick hash table
 * \param[in] p The probability of a random error
 * \param[in] t_delete How many timesteps in the past to keep before deletion 
 * \param[in] recipe The \ref recipe for the quantum computation 
 *
 * \callgraph
 * \callergraph
 *
 * \return The newly created \ref dp_qc
 */
DP_QC *dp_create_dp_qc(int s0, int s1, int de_ht_size, int stick_ht_size, double p, int t_delete, RECIPE *recipe, char *ems) {
	DP_QC *dp_qc;
	QC *qc;
	
	dp_qc = (DP_QC *)my_malloc(sizeof(DP_QC));
	dp_qc->qc = qc = (QC *)qc_create_qc(s0, s1, de_ht_size, stick_ht_size, t_delete, recipe);
	dp_qc->p = p;

	// Load the default error models
	dp_init_dp_qc(dp_qc, ems);

	return dp_qc;
}

/**
 * \brief 
 * 
 * \param[in] 
 */
DP_QC *dp_create_dp_qc_adv(int s0, int s1, int de_ht_size, int stick_ht_size, double p, int t_delete, RECIPE_ADV *recipe, char *ems) {
	DP_QC *dp_qc;
	QC *qc;

	dp_qc = (DP_QC *)my_malloc(sizeof(DP_QC));
	dp_qc->qc = qc = (QC *)qc_create_qc_adv(s0, s1, de_ht_size, stick_ht_size, t_delete, recipe);
	dp_qc->p = p;

	// Load the default error models
	dp_init_dp_qc(dp_qc, ems);

	return dp_qc;
}

/**
 * \brief 
 * 
 * \param[in] 
 */
void dp_init_dp_qc(DP_QC *dp_qc, char *ems) {
	dp_qc->init_X_id = dp_load_error_model(dp_qc, ems, "init_X_em");
	dp_qc->init_Z_id = dp_load_error_model(dp_qc, ems, "init_Z_em");
	dp_qc->H_id = dp_load_error_model(dp_qc, ems, "H_em");
	dp_qc->cnot_id = dp_load_error_model(dp_qc, ems, "cnot_em");
	dp_qc->cZ_id = dp_load_error_model(dp_qc, ems, "cZ_em");
	dp_qc->swap_id = dp_load_error_model(dp_qc, ems, "swap_em");
	dp_qc->meas_X_id = dp_load_error_model(dp_qc, ems, "meas_X_em");
	dp_qc->meas_Z_id = dp_load_error_model(dp_qc, ems, "meas_Z_em");
	dp_qc->iden_init_X_id = dp_load_error_model(dp_qc, ems, "iden_init_X_em");
	dp_qc->iden_init_Z_id = dp_load_error_model(dp_qc, ems, "iden_init_Z_em");
	dp_qc->iden_H_id = dp_load_error_model(dp_qc, ems, "iden_H_em");
	dp_qc->iden_cnot_id = dp_load_error_model(dp_qc, ems, "iden_cnot_em");
	dp_qc->iden_cZ_id = dp_load_error_model(dp_qc, ems, "iden_cZ_em");
	dp_qc->iden_swap_id = dp_load_error_model(dp_qc, ems, "iden_swap_em");
	dp_qc->iden_meas_X_id = dp_load_error_model(dp_qc, ems, "iden_meas_X_em");
	dp_qc->iden_meas_Z_id = dp_load_error_model(dp_qc, ems, "iden_meas_Z_em");
}

/**
 * \brief Loads an error model
 * 
 * \param[in] dp_qc
 * \param[in] ems The path to the folder containing the error models.
 * \param[in] filename The filename of the error model to load
 */
int dp_load_error_model(DP_QC *dp_qc, char *ems, const char *filename) {
	ERROR_MODEL *em;
	int ems_len = strlen(ems);
	em = qc_create_error_model(strcat(ems, filename));
	ems[ems_len] = '\0';
	return qc_insert_error_model(dp_qc->qc, em);
}

/**
 * \brief Copies a \ref dp_qc
 * 
 * \param[in] dp_qc The \ref dp_qc to be copied
 *
 * \return The newly copied \ref dp_qc
 */
DP_QC *dp_copy_dp_qc(DP_QC *dp_qc) {
	DP_QC *dp_qc2;

	dp_qc2 = (DP_QC *)my_malloc(sizeof(DP_QC));
	*dp_qc2 = *dp_qc;
	dp_qc2->qc = qc_copy_qc(dp_qc->qc);

	return dp_qc2;
}

/**
 * \brief Frees a \ref dp_qc
 * 
 * \param[in] dp_qc The \ref dp_qc to be freed
 */
void dp_free_dp_qc(DP_QC *dp_qc) {
	qc_free_qc(dp_qc->qc);
	free(dp_qc);
}

/**
 * \brief Frees a copy of a \ref dp_qc
 * 
 * \param[in] dp_qc The copied \ref dp_qc to be freed
 */
void dp_free_dp_qc_copy(DP_QC *dp_qc) {
	qc_free_qc_copy(dp_qc->qc);
	free(dp_qc);
}

/**
 * \brief Whether or not an e contains an X component
 * 
 * \param[in] e The e to check for an X component
 *
 * \return true if the e has an X component, false otherwise
 */
bool dp_contains_X(int e) {
	if ((e&X) != 0) return true;
	return false;
}

/**
 * \brief Whether or not an e contains a Z component
 * 
 * \param[in] e The e to check for a Z component
 *
 * \return true if the e has a Z component, false otherwise
 */
bool dp_contains_Z(int e) {
	if ((e & Z) != 0) {
	   return true;
	}
	return false;
}

/**
 * \brief Multiplies two es together
 * 
 * \param[in] e1 The first e
 * \param[in] e2 The second e
 *
 * \return The exclusive-or result of the two es 
 */
int dp_multiply_es(int e1, int e2) {
	return e1^e2;
}



// INIT_X FUNCTIONS

/**
 * \brief Initialize an X qubit
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to initialize
 * \param[in] q The \ref qubit to initialize
 */
void dp_init_X(DP_QC *dp_qc, QUBIT *q) {
	qc_init_qubit(dp_qc->qc, q, dp_qc->init_X_id, dp_qc->p, 
		dp_multiply_es);
}



// INIT_Z FUNCTIONS

/**
 * \brief Initialize a Z qubit
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to initialize
 * \param[in] q The \ref qubit to initialize
 */
void dp_init_Z(DP_QC *dp_qc, QUBIT *q) {
	#ifdef DEBUG
		printf("init_Z (%d, %d, %ld)\n", q->i, q->j, q->t);
	#endif
	qc_init_qubit(dp_qc->qc, q, dp_qc->init_Z_id, dp_qc->p, 
		dp_multiply_es);
}



// IDEN FUNCTIONS

/**
 * \brief Applies an identity gate
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 * \param[in] gate The id of the identity gate corresponding to the gate to
 * perform the equivalent identity of (eg. dp_qc->iden_init_X_id to perform an
 * identity corresponding to an init_X)
 */
void dp_iden(DP_QC *dp_qc, QUBIT *q, int gate) {
	#ifdef DEBUG
		printf("iden (%d) (%d, %d, %ld)\n", gate, q->i, q->j, q->t);
	#endif
	qc_single_qubit_unitary(dp_qc->qc, q, gate, dp_qc->p,
		dp_iden_transform_e, dp_multiply_es, qc_iden_transform_error);
}

/**
 * \brief Applies an identity gate of duration init_X
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dp_iden_init_X(DP_QC *dp_qc, QUBIT *q) {
	dp_iden(dp_qc, q, dp_qc->iden_init_X_id);
}

/**
 * \brief Applies an identity gate of duration init_Z
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dp_iden_init_Z(DP_QC *dp_qc, QUBIT *q) {
	dp_iden(dp_qc, q, dp_qc->iden_init_Z_id);
}

/**
 * \brief Applies an identity gate of duration H
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dp_iden_H(DP_QC *dp_qc, QUBIT *q) {
	dp_iden(dp_qc, q, dp_qc->iden_H_id);
}

/**
 * \brief Applies an identity gate of duration cnot
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dp_iden_cnot(DP_QC *dp_qc, QUBIT *q) {
	dp_iden(dp_qc, q, dp_qc->iden_cnot_id);
}

/**
 * \brief Applies an identity gate of duration cZ
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dp_iden_cZ(DP_QC *dp_qc, QUBIT *q) {
	dp_iden(dp_qc, q, dp_qc->iden_cZ_id);
}

/**
 * \brief Applies an identity gate of duration swap
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dp_iden_swap(DP_QC *dp_qc, QUBIT *q) {
	dp_iden(dp_qc, q, dp_qc->iden_swap_id);
}

/**
 * \brief Applies an identity gate of duration meas_X
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dp_iden_meas_X(DP_QC *dp_qc, QUBIT *q) {
	dp_iden(dp_qc, q, dp_qc->iden_meas_X_id);
}

/**
 * \brief Applies an identity gate of duration meas_Z
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dp_iden_meas_Z(DP_QC *dp_qc, QUBIT *q) {
	dp_iden(dp_qc, q, dp_qc->iden_meas_Z_id);
}

/**
 * \brief Transforms a \ref qubit based on a single qubit identity gate 
 * 
 * \param[in] q The \ref qubit to transform
 */
void dp_iden_transform_e(__attribute__((unused)) QUBIT *q) {
}




// DEAD FUNCTIONS

/**
 * \brief Applies a dead gate (essentially does nothing)
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 * \param[in] gate The id of the gate to perform a dead gate of the same
 * duration (eg. dp_qc->init_X to perform a dead gate or duration init_X) */
void dp_dead(DP_QC *dp_qc, QUBIT *q, int gate) {
	QC *qc;
	ERROR_MODEL *em;
	#ifdef DEBUG
		printf("dead (%d) (%d, %d, %ld)\n", gate, q->i, q->j, q->t);
	#endif
	qc = dp_qc->qc;
	em = qc->ems[gate];
	q->t += em->duration;
}

/**
 * \brief Applies a dead gate (essentially does nothing) of duration init_X
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit
 * \param[in] q The \ref qubit to perform a dead gate on
 */
void dp_dead_init_X(DP_QC *dp_qc, QUBIT *q) {
	dp_dead(dp_qc, q, dp_qc->init_X_id);
}

/**
 * \brief Applies a dead gate (essentially does nothing) of duration init_Z
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit
 * \param[in] q The \ref qubit to perform a dead gate on
 */
void dp_dead_init_Z(DP_QC *dp_qc, QUBIT *q) {
	dp_dead(dp_qc, q, dp_qc->init_Z_id);
}

/**
 * \brief Applies a dead gate (essentially does nothing) of duration H
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit
 * \param[in] q The \ref qubit to perform a dead gate on
 */
void dp_dead_H(DP_QC *dp_qc, QUBIT *q) {
	dp_dead(dp_qc, q, dp_qc->H_id);
}

/**
 * \brief Applies a dead gate (essentially does nothing) of duration cnot
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit
 * \param[in] q The \ref qubit to perform a dead gate on
 */
void dp_dead_cnot(DP_QC *dp_qc, QUBIT *q) {
	dp_dead(dp_qc, q, dp_qc->cnot_id);
}

/**
 * \brief Applies a dead gate (essentially does nothing) of duration cZ
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit
 * \param[in] q The \ref qubit to perform a dead gate on
 */
void dp_dead_cZ(DP_QC *dp_qc, QUBIT *q) {
	dp_dead(dp_qc, q, dp_qc->cZ_id);
}

/**
 * \brief Applies a dead gate (essentially does nothing) of duration swap
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit
 * \param[in] q The \ref qubit to perform a dead gate on
 */
void dp_dead_swap(DP_QC *dp_qc, QUBIT *q) {
	dp_dead(dp_qc, q, dp_qc->swap_id);
}

/**
 * \brief Applies a dead gate (essentially does nothing) of duration meas_X
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit
 * \param[in] q The \ref qubit to perform a dead gate on
 */
void dp_dead_meas_X(DP_QC *dp_qc, QUBIT *q) {
	dp_dead(dp_qc, q, dp_qc->meas_X_id);
}

/**
 * \brief Applies a dead gate (essentially does nothing) of duration meas_Z
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit
 * \param[in] q The \ref qubit to perform a dead gate on
 */
void dp_dead_meas_Z(DP_QC *dp_qc, QUBIT *q) {
	dp_dead(dp_qc, q, dp_qc->meas_Z_id);
}



// H FUNCTIONS

/**
 * \brief Applies a Hadamard gate to a single \ref qubit 
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dp_H(DP_QC *dp_qc, QUBIT *q) {
	#ifdef DEBUG
		printf("H (%d, %d, %ld)\n", q->i, q->j, q->t);
	#endif
	qc_single_qubit_unitary(dp_qc->qc, q, dp_qc->H_id, dp_qc->p,
		dp_H_transform_e, dp_multiply_es, qc_H_transform_error);
}

/**
 * \brief Transforms a \ref qubit based on a single qubit Hadamard gate 
 * 
 * \param[in] q The \ref qubit to transform
 */
void dp_H_transform_e(QUBIT *q) {
	if (q->e == X) {
		q->e = Z;
	}
	else if (q->e == Z) {
		q->e = X;
	}
}



// CNOT FUNCTIONS

/**
 * \brief Applies a controlled not gate to two \ref qubit%s 
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dp_cnot(DP_QC *dp_qc, QUBIT *q1, QUBIT *q2) {
	#ifdef DEBUG
		printf("cnot (%d, %d, %ld) (%d, %d, %ld)\n", q1->i, q1->j, q1->t, q2->i, q2->j, q2->t);
	#endif
	qc_two_qubit_unitary(dp_qc->qc, q1, q2, dp_qc->cnot_id, dp_qc->p,
		dp_cnot_transform_es, dp_multiply_es, qc_cnot_transform_errors);
}

/**
 * \brief Transforms two \ref qubit%s based on a two qubit controlled not gate
 * 
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dp_cnot_transform_es(QUBIT *q1, QUBIT *q2) {
	if (q1->e & X) {
		q2->e ^= X;
	}
	if (q2->e & Z) {
		q1->e ^= Z;
	}
}



// CZ FUNCTIONS

/**
 * \brief Applies a controlled Z gate to two \ref qubit%s 
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dp_cZ(DP_QC *dp_qc, QUBIT *q1, QUBIT *q2) {
	#ifdef DEBUG
		printf("cZ (%d, %d, %ld) (%d, %d, %ld)\n", q1->i, q1->j, q1->t, q2->i, q2->j, q2->t);
	#endif
	qc_two_qubit_unitary(dp_qc->qc, q1, q2, dp_qc->cZ_id, dp_qc->p,
		dp_cZ_transform_es, dp_multiply_es, qc_cZ_transform_errors);
}

/**
 * \brief Transforms two \ref qubit%s based on a two qubit controlled Z gate
 * 
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dp_cZ_transform_es(QUBIT *q1, QUBIT *q2) {
	if (q1->e & X) {
		q2->e ^= Z;
	}
	if (q2->e & X) {
		q1->e ^= Z;
	}
}



// SWAP FUNCTIONS

/**
 * \brief Applies a swap gate to two \ref qubit%s 
 * 
 * \param[in] dp_qc The \ref dp_qc containing the \ref qubit to transform
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dp_swap(DP_QC *dp_qc, QUBIT *q1, QUBIT *q2) {
	qc_two_qubit_unitary(dp_qc->qc, q1, q2, dp_qc->swap_id, dp_qc->p,
		dp_swap_transform_es, dp_multiply_es, qc_swap_transform_errors);
}

/**
 * \brief Transforms two \ref qubit%s based on a two qubit swap gate
 * 
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dp_swap_transform_es(QUBIT *q1, QUBIT *q2) {
	int i;

	i = q1->e;
	q1->e = q2->e;
	q2->e = i;
}



// MEAS_X FUNCTIONS

/**
 * \brief Measures a \ref qubit in the X basis
 * 
 * \param[in] dp_qc The \ref dp_qc that contains the \ref qubit
 * \param[in] q The \ref qubit to be measured
 * \param[in] set1 The first \ref set the \ref qubit belongs to. 
 * \param[in] set2 The second \ref set the \ref qubit belongs to. 
 *
 * \return The measurement result (-1, 0 or 1)
 */
int dp_meas_X(DP_QC *dp_qc, QUBIT *q, SET *set1, SET *set2) {
	return qc_meas_qubit(dp_qc->qc, q, dp_qc->meas_X_id, X,
		dp_qc->p, set1, set2, dp_multiply_es, dp_meas_X_sim);
}

/**
 * \brief Simulates an X measurement
 * 
 * \param[in] q The \ref qubit to measure
 *
 * \return -1 If the \ref qubit contains no X component of error, 1 otherwise
 */
int dp_meas_X_sim(QUBIT *q) {
	if ((q->e & Z) == Z) {
		return -1;
	}
	return 1;
}

// MEAS_Z FUNCTIONS

/**
 * \brief Measures a \ref qubit in the Z basis
 * 
 * \param[in] dp_qc The \ref dp_qc that contains the \ref qubit
 * \param[in] q The \ref qubit to be measured
 * \param[in] set1 The first \ref set the \ref qubit belongs to. 
 * \param[in] set2 The second \ref set the \ref qubit belongs to. 
 *
 * \return The measurement result (-1, 0 or 1)
 */
int dp_meas_Z(DP_QC *dp_qc, QUBIT *q, SET *set1, SET *set2) {
	#ifdef DEBUG
		printf("meas_Z (%d, %d, %ld)\n", q->i, q->j, q->t);
	#endif
	return qc_meas_qubit(dp_qc->qc, q, dp_qc->meas_Z_id, Z,
		dp_qc->p, set1, set2, dp_multiply_es, dp_meas_Z_sim);
}

/**
 * \brief Simulates a Z measurement
 * 
 * \param[in] q The \ref qubit to measure
 *
 * \return -1 If the \ref qubit contains no Z component of error, 1 otherwise
 */
int dp_meas_Z_sim(QUBIT *q) {
	if ((q->e & X) == X) {
		return -1;
	}
	return 1;
}

