#ifndef DEPOLAR_LOSS_H
#define DEPOLAR_LOSS_H

#include "../qc/qc.h"

typedef struct dpl_qc DPL_QC;

struct dpl_qc {
	/** The \ref qc for the \ref dpl_qc */
	QC *qc;
	
	/** The probability of a random error */
	double p;

	double p_loss;
	double p_depolar;
	int has_loss; 

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
	/** The id of the iden_cZ error model */
	int iden_cZ_id;
	/** The id of the loss_init_X error model */
	int loss_init_X_id;
	/** The id of the loss_init_Z error model */
	int loss_init_Z_id;
	/** The id of the loss_H error model */
	int loss_H_id;
	/** The id of the loss_cnot error model */
	int loss_cnot_id;
	/** The id of the loss_cZ error model */
	int loss_cZ_id;
	/** The id of the loss_swap error model */
	int loss_swap_id;
	/** The id of the loss_meas_X error model */
	int loss_meas_X_id;
	/** The id of the loss_meas_Z error model */
	int loss_meas_Z_id;
	/** The id of the loss_iden_cZ error model */
	int loss_iden_cZ_id;
};

// DPL_QC FUNCTIONS
DPL_QC *dpl_create_dpl_qc(int s0, int s1, int de_ht_size, int stick_ht_size, double p, double p_loss, double p_depolar, int t_delete, RECIPE_ADV *recipe, char *ems, int has_loss);
int dpl_load_error_model(DPL_QC *dpl_qc, char *ems, const char *filename);
DPL_QC *dpl_copy_dpl_qc(DPL_QC *dpl_qc);
void dpl_free_dpl_qc(DPL_QC *dpl_qc);
void dpl_free_dpl_qc_copy(DPL_QC *dpl_qc);
bool dpl_contains_X(int e);
bool dpl_contains_Z(int e);
bool dpl_is_lost(int e);

void dpl_save_rng(); 
void dpl_restore_rng(); 

int dpl_multiply_es(int e1, int e2);
void dpl_depolar(QUBIT *q, double p);

// INIT_X FUNCTIONS
void dpl_init_X(DPL_QC *dpl_qc, QUBIT *q);

// INIT_Z FUNCTIONS
void dpl_init_Z(DPL_QC *dpl_qc, QUBIT *q);

// LOSS FUNCTIONS
void dpl_set_e(QUBIT *q, int e);
void dpl_loss_force(DPL_QC *dpl_qc, QUBIT *q, int gate);
void dpl_loss(DPL_QC *dpl_qc, QUBIT *q, int gate);
void dpl_loss_init_X(DPL_QC *dpl_qc, QUBIT *q);
void dpl_loss_init_Z(DPL_QC *dpl_qc, QUBIT *q);
void dpl_loss_H(DPL_QC *dpl_qc, QUBIT *q);
void dpl_loss_cnot(DPL_QC *dpl_qc, QUBIT *q);
void dpl_loss_cZ(DPL_QC *dpl_qc, QUBIT *q);
void dpl_loss_swap(DPL_QC *dpl_qc, QUBIT *q);
void dpl_loss_meas_X(DPL_QC *dpl_qc, QUBIT *q);
void dpl_loss_meas_Z(DPL_QC *dpl_qc, QUBIT *q);
void dpl_loss_iden_cZ(DPL_QC *dpl_qc, QUBIT *q);
void dpl_loss_transform_e(__attribute__((unused)) QUBIT *q); 

// IDEN FUNCTIONS
void dpl_iden(DPL_QC *dpl_qc, QUBIT *q, int gate);
void dpl_iden_cZ(DPL_QC *dpl_qc, QUBIT *q);
void dpl_iden_transform_e(QUBIT *q);

// H FUNCTIONS
void dpl_H(DPL_QC *dpl_qc, QUBIT *q);
void dpl_H_transform_e(QUBIT *q);

// CNOT FUNCTIONS
void dpl_cnot(DPL_QC *dpl_qc, QUBIT *q1, QUBIT *q2);
void dpl_cnot_transform_es(QUBIT *q1, QUBIT *q2);

// CZ FUNCTIONS
void dpl_cZ(DPL_QC *dpl_qc, QUBIT *q1, QUBIT *q2);
void dpl_cZ_transform_es(QUBIT *q1, QUBIT *q2);

// SWAP FUNCTIONS
void dpl_swap(DPL_QC *dpl_qc, QUBIT *q1, QUBIT *q2);
void dpl_swap_transform_es(QUBIT *q1, QUBIT *q2);

// MEAS_X FUNCTIONS
int dpl_meas_X(DPL_QC *dpl_qc, QUBIT *q, SET *set1, SET *set2);
int dpl_meas_X_sim(QUBIT *q);

// MEAS_Z FUNCTIONS
int dpl_meas_Z(DPL_QC *dpl_qc, QUBIT *q, SET *set1, SET *set2);
int dpl_meas_Z_sim(QUBIT *q);

#endif
