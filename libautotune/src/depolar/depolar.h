#ifndef DEPOLAR_H
#define DEPOLAR_H

#include "../qc/qc.h"

typedef struct dp_qc DP_QC;

struct dp_qc {
	/** The \ref qc for the \ref dp_qc */
	QC *qc;
	
	/** The probability of a random error */
	double p;

	/** The id of the init_X error model */
	int init_X_id;
	/** The id of the init_Z error model */
	int init_Z_id;
	/** The id of the H error model */
	int H_id;
	/** The id of the cnot error model */
	int cnot_id;
	/** The id of the cZ error model */
	int cZ_id;
	/** The id of the swap error model */
	int swap_id;
	/** The id of the meas_X error model */
	int meas_X_id;
	/** The id of the meas_Z error model */
	int meas_Z_id;
	/** The id of the iden_init_X error model */
	int iden_init_X_id;
	/** The id of the iden_init_Z error model */
	int iden_init_Z_id;
	/** The id of the iden_H error model */
	int iden_H_id;
	/** The id of the iden_cnot error model */
	int iden_cnot_id;
	/** The id of the iden_cZ error model */
	int iden_cZ_id;
	/** The id of the iden_swap error model */
	int iden_swap_id;
	/** The id of the iden_meas_X error model */
	int iden_meas_X_id;
	/** The id of the iden_meas_Z error model */
	int iden_meas_Z_id;
};

// DP_QC FUNCTIONS
DP_QC *dp_create_dp_qc(int s0, int s1, int de_ht_size, int stick_ht_size, double p, int t_delete, RECIPE *recipe, char *ems);
DP_QC *dp_create_dp_qc_adv(int s0, int s1, int de_ht_size, int stick_ht_size, double p, int t_delete, RECIPE_ADV *recipe, char *ems);
void dp_init_dp_qc(DP_QC *dp_qc, char *ems);
int dp_load_error_model(DP_QC *dp_qc, char *ems, const char *filename);
DP_QC *dp_copy_dp_qc(DP_QC *dp_qc);
void dp_free_dp_qc(DP_QC *dp_qc);
void dp_free_dp_qc_copy(DP_QC *dp_qc);
bool dp_contains_X(int e);
bool dp_contains_Z(int e);
int dp_multiply_es(int e1, int e2);

// INIT_X FUNCTIONS
void dp_init_X(DP_QC *dp_qc, QUBIT *q);

// INIT_Z FUNCTIONS
void dp_init_Z(DP_QC *dp_qc, QUBIT *q);

// IDEN FUNCTIONS
void dp_iden(DP_QC *dp_qc, QUBIT *q, int gate);
void dp_iden_init_X(DP_QC *dp_qc, QUBIT *q);
void dp_iden_init_Z(DP_QC *dp_qc, QUBIT *q);
void dp_iden_H(DP_QC *dp_qc, QUBIT *q);
void dp_iden_cnot(DP_QC *dp_qc, QUBIT *q);
void dp_iden_cZ(DP_QC *dp_qc, QUBIT *q);
void dp_iden_swap(DP_QC *dp_qc, QUBIT *q);
void dp_iden_meas_X(DP_QC *dp_qc, QUBIT *q);
void dp_iden_meas_Z(DP_QC *dp_qc, QUBIT *q);
void dp_iden_transform_e(QUBIT *q);

// DEAD FUNCTIONS
void dp_dead(DP_QC *dp_qc, QUBIT *q, int gate);
void dp_dead_init_X(DP_QC *dp_qc, QUBIT *q);
void dp_dead_init_Z(DP_QC *dp_qc, QUBIT *q);
void dp_dead_H(DP_QC *dp_qc, QUBIT *q);
void dp_dead_cnot(DP_QC *dp_qc, QUBIT *q);
void dp_dead_cZ(DP_QC *dp_qc, QUBIT *q);
void dp_dead_swap(DP_QC *dp_qc, QUBIT *q);
void dp_dead_meas_X(DP_QC *dp_qc, QUBIT *q);
void dp_dead_meas_Z(DP_QC *dp_qc, QUBIT *q);

// H FUNCTIONS
void dp_H(DP_QC *dp_qc, QUBIT *q);
void dp_H_transform_e(QUBIT *q);

// CNOT FUNCTIONS
void dp_cnot(DP_QC *dp_qc, QUBIT *q1, QUBIT *q2);
void dp_cnot_transform_es(QUBIT *q1, QUBIT *q2);

// CZ FUNCTIONS
void dp_cZ(DP_QC *dp_qc, QUBIT *q1, QUBIT *q2);
void dp_cZ_transform_es(QUBIT *q1, QUBIT *q2);

// SWAP FUNCTIONS
void dp_swap(DP_QC *dp_qc, QUBIT *q1, QUBIT *q2);
void dp_swap_transform_es(QUBIT *q1, QUBIT *q2);

// MEAS_X FUNCTIONS
int dp_meas_X(DP_QC *dp_qc, QUBIT *q, SET *set1, SET *set2);
int dp_meas_X_sim(QUBIT *q);

// MEAS_Z FUNCTIONS
int dp_meas_Z(DP_QC *dp_qc, QUBIT *q, SET *set1, SET *set2);
int dp_meas_Z_sim(QUBIT *q);

#endif
