#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "depolar_loss.h"
#include "../memory/memory.h"
#include "../random/random.h"

//#define DEBUG 1

static RNG *dpl_qc_rng;

// DPL_QC FUNCTIONS

/**
 * \brief Creates a new \ref dpl_qc
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
 * \return The newly created \ref dpl_qc
 */
DPL_QC *dpl_create_dpl_qc(int s0, int s1, int de_ht_size, int stick_ht_size, double p, double p_loss, double p_depolar, int t_delete, RECIPE_ADV *recipe, char *ems, int has_loss) {
	DPL_QC *dpl_qc;
	QC *qc;
	
	dpl_qc_rng = rng_create();
	rng_randoms(dpl_qc_rng, &s0, &s1);

	dpl_qc = (DPL_QC *)my_malloc(sizeof(DPL_QC));
	dpl_qc->qc = qc = (QC *)qc_create_qc_adv(s0, s1, de_ht_size, stick_ht_size, t_delete, recipe);
	dpl_qc->p = p;

	if (has_loss) {
		dpl_qc->p_loss = p_loss;
		dpl_qc->p_depolar = p_depolar;
	} else {
		dpl_qc->p_loss = 0;
		dpl_qc->p_depolar = 0;
	}
	dpl_qc->has_loss = has_loss;

	// Load the default error models
	dpl_qc->init_X_id = dpl_load_error_model(dpl_qc, ems, "init_X_em");
	dpl_qc->init_Z_id = dpl_load_error_model(dpl_qc, ems, "init_Z_em");
	dpl_qc->H_id = dpl_load_error_model(dpl_qc, ems, "H_em");
	dpl_qc->cnot_id = dpl_load_error_model(dpl_qc, ems, "cnot_em");
	dpl_qc->cZ_id = dpl_load_error_model(dpl_qc, ems, "cZ_em");
	dpl_qc->swap_id = dpl_load_error_model(dpl_qc, ems, "swap_em");
	dpl_qc->meas_X_id = dpl_load_error_model(dpl_qc, ems, "meas_X_em");
	dpl_qc->meas_Z_id = dpl_load_error_model(dpl_qc, ems, "meas_Z_em");
    // TCS only uses one identity gate. More will need to be added if the loss
    // model is used with other correction schemes in the future.
	dpl_qc->iden_cZ_id = dpl_load_error_model(dpl_qc, ems, "iden_cZ_em");

	// Load the loss error models
	dpl_qc->loss_init_X_id = dpl_load_error_model(dpl_qc, ems, "loss_init_X_em");
	dpl_qc->loss_init_Z_id = dpl_load_error_model(dpl_qc, ems, "loss_init_Z_em");
	dpl_qc->loss_H_id = dpl_load_error_model(dpl_qc, ems, "loss_H_em");
	dpl_qc->loss_cnot_id = dpl_load_error_model(dpl_qc, ems, "loss_cnot_em");
	dpl_qc->loss_cZ_id = dpl_load_error_model(dpl_qc, ems, "loss_cZ_em");
	dpl_qc->loss_swap_id = dpl_load_error_model(dpl_qc, ems, "loss_swap_em");
	dpl_qc->loss_meas_X_id = dpl_load_error_model(dpl_qc, ems, "loss_meas_X_em");
	dpl_qc->loss_meas_Z_id = dpl_load_error_model(dpl_qc, ems, "loss_meas_Z_em");
	dpl_qc->loss_iden_cZ_id = dpl_load_error_model(dpl_qc, ems, "loss_iden_cZ_em");

	return dpl_qc;
}

/**
 * \brief Loads an error model
 * 
 * \param[in] dpl_qc
 * \param[in] ems The path to the folder containing the error models.
 * \param[in] filename The filename of the error model to load
 */
int dpl_load_error_model(DPL_QC *dpl_qc, char *ems, const char *filename) {
	ERROR_MODEL *em;
	int ems_len = strlen(ems);
	em = qc_create_error_model(strcat(ems, filename));
	ems[ems_len] = '\0';
	return qc_insert_error_model(dpl_qc->qc, em);
}

/**
 * \brief Copies a \ref dpl_qc
 * 
 * \param[in] dpl_qc The \ref dpl_qc to be copied
 *
 * \return The newly copied \ref dpl_qc
 */
DPL_QC *dpl_copy_dpl_qc(DPL_QC *dpl_qc) {
	DPL_QC *dpl_qc2;

	dpl_qc2 = (DPL_QC *)my_malloc(sizeof(DPL_QC));
	*dpl_qc2 = *dpl_qc;
	dpl_qc2->qc = qc_copy_qc(dpl_qc->qc);

	return dpl_qc2;
}

/**
 * \brief Frees a \ref dpl_qc
 * 
 * \param[in] dpl_qc The \ref dpl_qc to be freed
 */
void dpl_free_dpl_qc(DPL_QC *dpl_qc) {
	rng_free(dpl_qc_rng);
	qc_free_qc(dpl_qc->qc);
	free(dpl_qc);
}

/**
 * \brief Frees a copy of a \ref dpl_qc
 * 
 * \param[in] dpl_qc The copied \ref dpl_qc to be freed
 */
void dpl_free_dpl_qc_copy(DPL_QC *dpl_qc) {
	qc_free_qc_copy(dpl_qc->qc);
	free(dpl_qc);
}

/**
 * \brief Whether or not an e contains an X component
 * 
 * \param[in] e The e to check for an X component
 *
 * \return true if the e has an X component, false otherwise
 */
bool dpl_contains_X(int e) {
	if (e & X) {
        return true;
    }
    return false;
}

/**
 * \brief Whether or not an e contains a Z component
 * 
 * \param[in] e The e to check for a Z component
 *
 * \return true if the e has a Z component, false otherwise
 */

bool dpl_contains_Z(int e) {
	if (e & Z) {
	   return true;
	}
	return false;
}

/**
 * \brief Whether or not an e is lost
 * 
 * \param[in] e The e to check for loss
 *
 * \return true if the e is lost, false otherwise
 */
bool dpl_is_lost(int e) {
	if (e & L) {
	   return true;
	}
	return false;
}

void dpl_save_rng() {
	rng_save(dpl_qc_rng);
}

void dpl_restore_rng() {
	rng_restore(dpl_qc_rng);
}


/**
 * \brief Multiplies two es together
 * 
 * \param[in] e1 The first e
 * \param[in] e2 The second e
 *
 * \return The exclusive-or result of the two es 
 */
int dpl_multiply_es(int e1, int e2) {
	if (dpl_is_lost(e1) || dpl_is_lost(e2)) {
		return L;
	}
	return e1 ^ e2;
}

/**
 * \brief Apply depolarising noise to a qubit
 *
 * When interacting a qubit with a qubit is lost, the state of the lost state
 * is by definition unknown. We therefore have no idea if the state of that
 * qubit will imbue an X, a Z, a Y or even no error at all. We have weighted
 * each of these outcomes as equally likely.
 * 
 * \param[in] q The \ref qubit to be modified
 */
void dpl_depolar(QUBIT *q, double p) {
	double x;

	x = rng_randomd(dpl_qc_rng); 
	x = p * x;

	// No error
	if (x < 0.25) {
		return;
	}

	// X Error
	else if (x < 0.5) {
		q->e ^= X;
		return;
	}

	// Both X and Z error
	else if (x < 0.75) {
		q->e ^= Y;
		return;
	}

	else if (x <= 1.00) {
		// Z Error
		q->e ^= Z;
		return;
	}
}

// INIT_X FUNCTIONS

/**
 * \brief Initialize an X qubit
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to initialize
 * \param[in] q The \ref qubit to initialize
 */
void dpl_init_X(DPL_QC *dpl_qc, QUBIT *q) {
	#ifdef DEBUG
		printf("init_X (%d, %d, %ld)\n", q->i, q->j, q->t);
	#endif
	qc_init_qubit(dpl_qc->qc, q, dpl_qc->init_X_id, dpl_qc->p, dpl_multiply_es);
	dpl_loss_init_X(dpl_qc, q);
}



// INIT_Z FUNCTIONS

/**
 * \brief Initialize a Z qubit
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to initialize
 * \param[in] q The \ref qubit to initialize
 */
void dpl_init_Z(DPL_QC *dpl_qc, QUBIT *q) {
	#ifdef DEBUG
		printf("init_Z (%d, %d, %ld)\n", q->i, q->j, q->t);
	#endif
	qc_init_qubit(dpl_qc->qc, q, dpl_qc->init_Z_id, dpl_qc->p, dpl_multiply_es);
	dpl_loss_init_Z(dpl_qc, q);
}


// LOSS FUNCTIONS

void dpl_set_e(QUBIT *q, int e) {
	q->e = e;
}

void dpl_loss_force(DPL_QC *dpl_qc, QUBIT *q, int gate) {
	if (!dpl_qc->has_loss) {
		return;
	}

	qc_single_qubit_unitary(dpl_qc->qc, q, gate, 1,
		dpl_loss_transform_e, dpl_multiply_es, qc_loss_transform_error);
	printf("Loss: ");
	qc_print_qubit(q);
}

void dpl_loss(DPL_QC *dpl_qc, QUBIT *q, int gate) {
	if (!dpl_qc->has_loss) {
		return;
	}

	qc_single_qubit_unitary(dpl_qc->qc, q, gate, dpl_qc->p_loss,
		dpl_loss_transform_e, dpl_multiply_es, qc_loss_transform_error);
}

/**
 * \brief Applies an loss gate for init X
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dpl_loss_init_X(DPL_QC *dpl_qc, QUBIT *q) {
	dpl_loss(dpl_qc, q, dpl_qc->loss_init_X_id);
}

/**
 * \brief Applies an loss gate for init Z
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dpl_loss_init_Z(DPL_QC *dpl_qc, QUBIT *q) {
	dpl_loss(dpl_qc, q, dpl_qc->loss_init_Z_id);
}

/**
 * \brief Applies an loss gate for H
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dpl_loss_H(DPL_QC *dpl_qc, QUBIT *q) {
	dpl_loss(dpl_qc, q, dpl_qc->loss_H_id);
}

/**
 * \brief Applies an loss gate for cnot
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dpl_loss_cnot(DPL_QC *dpl_qc, QUBIT *q) {
	dpl_loss(dpl_qc, q, dpl_qc->loss_cnot_id);
}

/**
 * \brief Applies an loss gate for cZ
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dpl_loss_cZ(DPL_QC *dpl_qc, QUBIT *q) {
	dpl_loss(dpl_qc, q, dpl_qc->loss_cZ_id);
}

/**
 * \brief Applies an loss gate for swap
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dpl_loss_swap(DPL_QC *dpl_qc, QUBIT *q) {
	dpl_loss(dpl_qc, q, dpl_qc->loss_swap_id);
}

/**
 * \brief Applies an loss gate for meas X
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dpl_loss_meas_X(DPL_QC *dpl_qc, QUBIT *q) {
	dpl_loss(dpl_qc, q, dpl_qc->loss_meas_X_id);
}

/**
 * \brief Applies an loss gate for meas Z
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dpl_loss_meas_Z(DPL_QC *dpl_qc, QUBIT *q) {
	dpl_loss(dpl_qc, q, dpl_qc->loss_meas_Z_id);
}

/**
 * \brief Applies an loss gate for identity of length cZ
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dpl_loss_iden_cZ(DPL_QC *dpl_qc, QUBIT *q) {
	dpl_loss(dpl_qc, q, dpl_qc->loss_iden_cZ_id);
}

/**
 * \brief Transforms a \ref qubit based on a single qubit loss gate 
 * 
 * \param[in] q The \ref qubit to transform
 */
void dpl_loss_transform_e(__attribute__((unused)) QUBIT *q) {
}


// IDEN FUNCTIONS

/**
 * \brief Applies an identity gate
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 * \param[in] gate The id of the identity gate corresponding to the gate to
 * perform the equivalent identity of (eg. dpl_qc->iden_init_X_id to perform an
 * identity corresponding to an init_X)
 */
void dpl_iden(DPL_QC *dpl_qc, QUBIT *q, int gate) {
	#ifdef DEBUG
		printf("iden (%d) (%d, %d, %ld)\n", gate, q->i, q->j, q->t);
	#endif
	qc_single_qubit_unitary(dpl_qc->qc, q, gate, dpl_qc->p,
		dpl_iden_transform_e, dpl_multiply_es, qc_iden_transform_error);
}

/**
 * \brief Applies an identity gate of duration cZ
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dpl_iden_cZ(DPL_QC *dpl_qc, QUBIT *q) {
	dpl_iden(dpl_qc, q, dpl_qc->iden_cZ_id);
    dpl_loss_cZ(dpl_qc, q);
}

/**
 * \brief Transforms a \ref qubit based on a single qubit identity gate 
 * 
 * \param[in] q The \ref qubit to transform
 */
void dpl_iden_transform_e(__attribute__((unused)) QUBIT *q) {
}


// H FUNCTIONS

/**
 * \brief Applies a Hadamard gate to a single \ref qubit 
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q The \ref qubit to transform
 */
void dpl_H(DPL_QC *dpl_qc, QUBIT *q) {
	#ifdef DEBUG
		printf("H (%d, %d, %ld)\n", q->i, q->j, q->t);
	#endif
	qc_single_qubit_unitary(dpl_qc->qc, q, dpl_qc->H_id, dpl_qc->p,
		dpl_H_transform_e, dpl_multiply_es, qc_H_transform_error);
    dpl_loss_H(dpl_qc, q);
}

/**
 * \brief Transforms a \ref qubit based on a single qubit Hadamard gate 
 * 
 * \param[in] q The \ref qubit to transform
 */
void dpl_H_transform_e(QUBIT *q) {
	if (q->e == L) {
		return;
	}

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
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dpl_cnot(DPL_QC *dpl_qc, QUBIT *q1, QUBIT *q2) {
	#ifdef DEBUG
		printf("cnot (%d, %d, %ld) (%d, %d, %ld)\n", q1->i, q1->j, q1->t, q2->i, q2->j, q2->t);
	#endif
	qc_two_qubit_unitary(dpl_qc->qc, q1, q2, dpl_qc->cnot_id, dpl_qc->p,
		dpl_cnot_transform_es, dpl_multiply_es, qc_cnot_transform_errors);
    dpl_loss_cnot(dpl_qc, q1);
    dpl_loss_cnot(dpl_qc, q2);
}

/**
 * \brief Transforms two \ref qubit%s based on a two qubit controlled not gate
 * 
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dpl_cnot_transform_es(QUBIT *q1, QUBIT *q2) {
	if (q1->e & L) {
		// If both qubits are lost, do nothing
		if (q2->e & L) {
			return;
		}
		//dpl_depolar(q2);
	}
	else if (q2->e & L) {
		//dpl_depolar(q1);
	}
	else {
		if (q1->e & X) {
			q2->e ^= X;
		}
		if (q2->e & Z) {
			q1->e ^= Z;
		}
	}
}



// CZ FUNCTIONS

/**
 * \brief Applies a controlled Z gate to two \ref qubit%s 
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dpl_cZ(DPL_QC *dpl_qc, QUBIT *q1, QUBIT *q2) {
	#ifdef DEBUG
		printf("cZ (%d, %d, %ld) (%d, %d, %ld)\n", q1->i, q1->j, q1->t, q2->i, q2->j, q2->t);
	#endif
	ERROR_MODEL *em;

	if (dpl_qc->p_depolar && ((q1->e & L) || (q2->e & L))) {
		em = dpl_qc->qc->ems[dpl_qc->cZ_id];
		q1->t += em->duration;
		q2->t += em->duration;

		// If q1 or q2 isn't lost, then apply depolarising noise
		if (!(q1->e & L)) {
			dpl_depolar(q1, dpl_qc->p_depolar);
		}
		if (!(q2->e & L)) {
			dpl_depolar(q2, dpl_qc->p_depolar);
		} 
	} else {
		qc_two_qubit_unitary(dpl_qc->qc, q1, q2, dpl_qc->cZ_id, dpl_qc->p,
			dpl_cZ_transform_es, dpl_multiply_es, qc_cZ_transform_errors);
	}

    dpl_loss_cZ(dpl_qc, q1);
    dpl_loss_cZ(dpl_qc, q2);
}

/**
 * \brief Transforms two \ref qubit%s based on a two qubit controlled Z gate
 * 
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dpl_cZ_transform_es(QUBIT *q1, QUBIT *q2) {
	if (q1->e & X) {
		if (!(q2->e & L)) {
			q2->e ^= Z;
		}
	}
	if (q2->e & X) {
		if (!(q1->e & L)) {
			q1->e ^= Z;
		}
	}
}



// SWAP FUNCTIONS

/**
 * \brief Applies a swap gate to two \ref qubit%s 
 * 
 * \param[in] dpl_qc The \ref dpl_qc containing the \ref qubit to transform
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dpl_swap(DPL_QC *dpl_qc, QUBIT *q1, QUBIT *q2) {
	qc_two_qubit_unitary(dpl_qc->qc, q1, q2, dpl_qc->swap_id, dpl_qc->p,
		dpl_swap_transform_es, dpl_multiply_es, qc_swap_transform_errors);
    dpl_loss_swap(dpl_qc, q1);
    dpl_loss_swap(dpl_qc, q2);
}

/**
 * \brief Transforms two \ref qubit%s based on a two qubit swap gate
 * 
 * \param[in] q1 The first \ref qubit to transform
 * \param[in] q2 The second \ref qubit to transform
 */
void dpl_swap_transform_es(QUBIT *q1, QUBIT *q2) {
	int i;

	i = q1->e;
	q1->e = q2->e;
	q2->e = i;
}



// MEAS_X FUNCTIONS

/**
 * \brief Measures a \ref qubit in the X basis
 * 
 * \param[in] dpl_qc The \ref dpl_qc that contains the \ref qubit
 * \param[in] q The \ref qubit to be measured
 * \param[in] set1 The first \ref set the \ref qubit belongs to. 
 * \param[in] set2 The second \ref set the \ref qubit belongs to. 
 *
 * \return The measurement result (-1, 0 or 1)
 */
int dpl_meas_X(DPL_QC *dpl_qc, QUBIT *q, SET *set1, SET *set2) {
    dpl_loss_meas_X(dpl_qc, q);
	return qc_meas_qubit(dpl_qc->qc, q, dpl_qc->meas_X_id, X,
		dpl_qc->p, set1, set2, dpl_multiply_es, dpl_meas_X_sim);
}

/**
 * \brief Simulates an X measurement
 * 
 * \param[in] q The \ref qubit to measure
 *
 * \return 0 if the \ref qubit has been lost. -1 If the \ref qubit contains a Z
 * component of error, 1 otherwise
 */
int dpl_meas_X_sim(QUBIT *q) {
    if (q->e & L) {
        return 0;
    }
	if (q->e & Z) {
		return -1;
	}
	return 1;
}

// MEAS_Z FUNCTIONS

/**
 * \brief Measures a \ref qubit in the Z basis
 * 
 * \param[in] dpl_qc The \ref dpl_qc that contains the \ref qubit
 * \param[in] q The \ref qubit to be measured
 * \param[in] set1 The first \ref set the \ref qubit belongs to. 
 * \param[in] set2 The second \ref set the \ref qubit belongs to. 
 *
 * \return The measurement result (-1, 0 or 1)
 */
int dpl_meas_Z(DPL_QC *dpl_qc, QUBIT *q, SET *set1, SET *set2) {
	#ifdef DEBUG
		printf("meas_Z (%d, %d, %ld)\n", q->i, q->j, q->t);
	#endif
    dpl_loss_meas_Z(dpl_qc, q);
	return qc_meas_qubit(dpl_qc->qc, q, dpl_qc->meas_Z_id, Z,
		dpl_qc->p, set1, set2, dpl_multiply_es, dpl_meas_Z_sim);
}

/**
 * \brief Simulates a Z measurement
 * 
 * \param[in] q The \ref qubit to measure
 *
 * \return -1 If the \ref qubit contains and X component of error, 1 otherwise
 */
int dpl_meas_Z_sim(QUBIT *q) {
    if (q->e & L) {
        return 0;
    }
	if (q->e & X) {
		return -1;
	}
	return 1;
}


