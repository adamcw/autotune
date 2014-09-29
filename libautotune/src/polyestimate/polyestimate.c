#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "polyestimate.h"
#include "../memory/memory.h"

#define AVERAGE 1

DATABASE *pe_load_database(char *fname) {
	FILE *in;
	STAT *stat;
	DATABASE *db;

   in = fopen(fname, "r");
   if (in == NULL) {
      fprintf(stderr, "Database file %s not found\n", fname);
      exit(EXIT_FAILURE);
   }

   db = NULL;
   stat = pe_load_stat(in);
   while (stat != NULL) {
		// pe_print_stat(stat);
		/*
		if (stat->pZ < stat->pX) {
			fprintf(stderr, "WARNING: low statistics in record ");
			pe_print_stat(stat);
		}
		*/
		db = ll_insert(db, stat);
		stat = pe_load_stat(in);
	}
	fclose(in);

	return db;
}

void pe_free_database(DATABASE *db) {
	ll_free(db, free);
}

FORM *pe_calculate_formulae_from_p(double p, DATABASE *db) {
	int i;
	double p8[NUM_EMS], p3x[3], p3z[3];

	for (i=0; i<NUM_EMS; i++) p8[i] = p;

	pe_convert_ps_to_ps3s(p8, p3x, p3z);

	return pe_create_form(p3x, p3z, db);
}

FORM *pe_calculate_formulae_from_3ps(double *p3, DATABASE *db) {
	int i;
	double p3x[3], p3z[3];

	for (i=0; i<3; i++) p3x[i] = p3z[i] = p3[i];

	return pe_create_form(p3x, p3z, db);
}

FORM *pe_calculate_formulae_from_6ps(double *p3x, double *p3z, DATABASE *db) {
	return pe_create_form(p3x, p3z, db);
}

FORM *pe_calculate_formulae_from_8ps(double *p8, DATABASE *db) {
	double p3x[3], p3z[3];

	pe_convert_ps_to_ps3s(p8, p3x, p3z);

	return pe_create_form(p3x, p3z, db);
}

FORM *pe_calculate_formulae_from_p_json(char *dir, DATABASE *db) {
	char fname[MAX_STRING_LENGTH];
	FILE *in;
	cJSON *json;
	double p8[NUM_EMS];

	if (dir[strlen(dir)-1] != '/') strcat(dir, "/");

	sprintf(fname, "%sp.json", dir);
	in = fopen(fname, "r");
	if (in == NULL) {
		fprintf(stderr, "%s not found.\n", fname);
		exit(EXIT_FAILURE);
	}

	json = load_json(in);
	fclose(in);
	if (json == NULL) {
		fprintf(stderr, "%s improperly formatted, load failed.\n", fname);
		exit(EXIT_FAILURE);
	}

	pe_convert_json_to_ps(json, p8);
	cJSON_Delete(json);

	return pe_calculate_formulae_from_8ps(p8, db);
}
	
FORM *pe_calculate_formulae_from_gates_json(char *dir, DATABASE *db) {
	char fname[MAX_STRING_LENGTH];
	FILE *in;
	cJSON *json;
	double p;

	if (dir[strlen(dir)-1] != '/') strcat(dir, "/");

	sprintf(fname, "%sgates.json", dir);
	in = fopen(fname, "r");
	if (in == NULL) {
		fprintf(stderr, "%s not found.\n", fname);
		exit(EXIT_FAILURE);
	}

	json = load_json(in);
	fclose(in);
	if (json == NULL) {
		fprintf(stderr, "%s improperly formatted, load failed.\n", fname);
		exit(EXIT_FAILURE);
	}

	p = pe_convert_json_to_ems_and_p(json, dir);
	cJSON_Delete(json);

	return pe_calculate_formulae_from_ems(p, dir, db);
}

FORM *pe_calculate_formulae_from_ems(double p, char *dir, DATABASE *db) {
	int i;
	char fname[MAX_STRING_LENGTH];
	ERROR_MODEL *em[NUM_EMS];
	double p3x[3], p3z[3];

	if (dir[strlen(dir)-1] != '/') strcat(dir, "/");

   sprintf(fname, "%scnot_em", dir);
   em[0] = qc_create_error_model(fname);
	if (em[0] == NULL) {
		fprintf(stderr, "%s missing.\n", fname);
		exit(EXIT_FAILURE);
	}

	sprintf(fname, "%sH_em", dir);
	em[1] = qc_create_error_model(fname);
	if (em[1] == NULL) {
		fprintf(stderr, "%s missing.\n", fname);
		exit(EXIT_FAILURE);
	}

	sprintf(fname, "%sinit_Z_em", dir);
	em[2] = qc_create_error_model(fname);
	if (em[2] == NULL) {
		fprintf(stderr, "%s missing.\n", fname);
		exit(EXIT_FAILURE);
	}

	sprintf(fname, "%smeas_Z_em", dir);
	em[3] = qc_create_error_model(fname);
	if (em[3] == NULL) {
		fprintf(stderr, "%s missing.\n", fname);
		exit(EXIT_FAILURE);
	}

	sprintf(fname, "%siden_cnot_em", dir);
	em[4] = qc_create_error_model(fname);
	if (em[4] == NULL) {
		fprintf(stderr, "%s missing.\n", fname);
		exit(EXIT_FAILURE);
	}

	sprintf(fname, "%siden_H_em", dir);
	em[5] = qc_create_error_model(fname);
	if (em[5] == NULL) {
		fprintf(stderr, "%s missing.\n", fname);
		exit(EXIT_FAILURE);
	}

	sprintf(fname, "%siden_init_Z_em", dir);
	em[6] = qc_create_error_model(fname);
	if (em[6] == NULL) {
		fprintf(stderr, "%s missing.\n", fname);
		exit(EXIT_FAILURE);
	}

	sprintf(fname, "%siden_meas_Z_em", dir);
	em[7] = qc_create_error_model(fname);
	if (em[7] == NULL) {
		fprintf(stderr, "%s missing.\n", fname);
		exit(EXIT_FAILURE);
	}

	pe_convert_ems_to_ps3s(p, em, p3x, p3z);

	for (i=0; i<NUM_EMS; i++) {
		qc_free_error_model(em[i]);
	}

	return pe_calculate_formulae_from_6ps(p3x, p3z, db);
}

void pe_free_formulae(FORM *f) {
	free(f);
}

double pe_calculate_logical_error(int d, int type, FORM *f) {
	return pe_p_L_estimate(d, f, type);
}

// Convert a json structure containing data describing error models
// for CNOT, Hadamard, initialization, measurement and identity gates
// for each of these non-trivial gates into 8 simple files and a
// global error rate p.
double pe_convert_json_to_ems_and_p(cJSON *json, char *ems_dir) {
	cJSON *ptr;
	double global_error;

	global_error = -1;

	ptr = json->child;
	while (ptr != NULL) {
		if (strcmp(ptr->string, "global_error") == 0) {
			global_error = ptr->valuedouble;
		}
		else if (strcmp(ptr->string, "gates") == 0) {
			pe_convert_json_to_ems(ptr, ems_dir);
		}
		ptr = ptr->next;
	}

	if (global_error <= 0 || global_error >= 1) {
		fprintf(stderr, "A value of global_error between 0 and 1 must be provided.\n");
		exit(EXIT_FAILURE);
	}

	return global_error;
}

void pe_convert_json_to_ems(cJSON *json, char *ems_dir) {
	int i, flag[NUM_EMS];
	cJSON *ptr;

	for (i=0; i<NUM_EMS; i++) {
		flag[i] = 0;
	}

	ptr = json->child;
	while (ptr != NULL) {
		flag[pe_convert_json_to_em(ptr, ems_dir)] = 1;
		ptr = ptr->next;
	}

	for (i=0; i<NUM_EMS; i++) {
		if (flag[i] == 0) {
			switch (i) {
				case 0:
					fprintf(stderr, "Data for CNOT missing.\n");
				case 1:
					fprintf(stderr, "Data for Hadamard missing.\n");
				case 2:
					fprintf(stderr, "Data for initialization missing.\n");
				case 3:
					fprintf(stderr, "Data for measurement missing.\n");
				case 4:
					fprintf(stderr, "Data for identity of duration CNOT missing.\n");
				case 5:
					fprintf(stderr, "Data for identity of duration Hadamard missing.\n");
				case 6:
					fprintf(stderr, "Data for identity of duration initialization missing.\n");
				case 7:
					fprintf(stderr, "Data for identity of duration measurement missing.\n");
			}
			exit(EXIT_FAILURE);
		}
	}
}

int pe_convert_json_to_em(cJSON *json, char *ems_dir) {
	char fname[MAX_STRING_LENGTH];
	int gate, num_qubits;
	FILE *out;
	cJSON *ptr, *error_ptr;
	double gate_time_ns, rel_error;

	gate = -1;
	num_qubits = -1;
	error_ptr = NULL;
	gate_time_ns = -1;
	rel_error = -1;

	if (strcmp(json->string, "CNOT") == 0) {
		sprintf(fname, "%scnot_em", ems_dir);
		gate = 0;
	}
	else if (strcmp(json->string, "Had") == 0) {
		sprintf(fname, "%sH_em", ems_dir);
		gate = 1;
	}
	else if (strcmp(json->string, "Init") == 0) {
		sprintf(fname, "%sinit_Z_em", ems_dir);
		gate = 2;
	}
	else if (strcmp(json->string, "Measure") == 0) {
		sprintf(fname, "%smeas_Z_em", ems_dir);
		gate = 3;
	}
	else if (strcmp(json->string, "IdCNOT") == 0) {
		sprintf(fname, "%siden_cnot_em", ems_dir);
		gate = 4;
	}
	else if (strcmp(json->string, "IdHad") == 0) {
		sprintf(fname, "%siden_H_em", ems_dir);
		gate = 5;
	}
	else if (strcmp(json->string, "IdInit") == 0) {
		sprintf(fname, "%siden_init_Z_em", ems_dir);
		gate = 6;
	}
	else if (strcmp(json->string, "IdMeasure") == 0) {
		sprintf(fname, "%siden_meas_Z_em", ems_dir);
		gate = 7;
	}

	if (gate < 0) {
		fprintf(stderr, "Unknown gate %s.\n", json->string);
		exit(EXIT_FAILURE);
	}

	out = fopen(fname, "w");

	ptr = json->child;
	while (ptr != NULL) {
		if (strcmp(ptr->string, "errors") == 0) {
			error_ptr = ptr;
		}
		else if (strcmp(ptr->string, "num_qubits") == 0) {
			num_qubits = ptr->valueint;
		}
		else if (strcmp(ptr->string, "gate_time_ns") == 0) {
			gate_time_ns = ptr->valuedouble;
		}
		else if (strcmp(ptr->string, "rel_error") == 0) {
			rel_error = ptr->valuedouble;
		}
		ptr = ptr->next;
	}

	if (error_ptr == NULL) {
		fprintf(stderr, "No errors found for gate %s.\n", json->string);
		exit(EXIT_FAILURE);
	}

	if (num_qubits == 0) {
		if (gate == 0) {
			fprintf(stderr, "%s should have num_qubits = 2.\n", json->string);
		}
		else {
			fprintf(stderr, "%s should have num_qubits = 1.\n", json->string);
		}
		exit(EXIT_FAILURE);
	}

	if (gate_time_ns < 0) {
		fprintf(stderr, "%s has negative or missing gate_time_ns.\n", json->string);
		exit(EXIT_FAILURE);
	}

	if (rel_error < 0) {
		fprintf(stderr, "%s has negative or missing rel_error.\n", json->string);
		exit(EXIT_FAILURE);
	}

	fprintf(out, "%d\n", num_qubits);
	fprintf(out, "%g\n", rel_error);
	fprintf(out, "%ld\n", cJSON_GetArraySize(error_ptr));

	ptr = error_ptr->child;
	while (ptr != NULL) {
		pe_convert_json_to_em_line(ptr, out);
		ptr = ptr->next;
	}

	if (gate_time_ns < 1) {
		fprintf(stderr, "WARNING: %s gate_time_ns will be truncated to 0.\n", json->string);
	}

	fprintf(out, "%ld\n", (long int)gate_time_ns);

	fclose(out);

	return gate;
}

void pe_convert_json_to_em_line(cJSON *json, FILE *out) {
	cJSON *ptr;
	long int i;

	ptr = json->child;
	while (ptr != NULL) {
		if (strcmp(ptr->string, "prob") == 0) {
			fprintf(out, "%ld", ptr->valueint);
			break;
		}
		ptr = ptr->next;
	}

	ptr = json->child;
	while (ptr != NULL) {
		if (strcmp(ptr->string, "pauli_indices") == 0) {
			ptr = ptr->child;
			while (ptr != NULL) {
				i = ptr->valueint;
				if (i == 2) i = 3;
				else if (i == 3) i = 2;
				fprintf(out, " %ld", i);
				ptr = ptr->next;
			}
			fprintf(out, "\n");
			break;
		}
		ptr = ptr->next;
	}
}

double pe_max(double a, double b, double c) {
	if (a < b) a = b;
	if (a < c) a = c;

	return a;
}

double pe_av(double a, double b, double c) {
	return (a+b+c)/3;
}

// Convert error model files for CNOT, Hadamard, initialization,
// measurement and corresponding identity gates into two sets of
// three simple error rates.
void pe_convert_ems_to_ps3s(double p, ERROR_MODEL **ems, double *ps3x, double *ps3z) {
	int i, a, b;
	double m;
	double px, py, pz;
	double pix, piy, piz, pxi, pxx, pxy, pxz, pyi, pyx, pyy, pyz, pzi, pzx, pzy, pzz;
	double ppix, ppxi, ppxx;
	double ppiz, ppzi, ppzz;

	pix = piy = piz = pxi = pxx = pxy = pxz = pyi = pyx = pyy = pyz = pzi = pzx = pzy = pzz = 0;

	for (i=0; i<ems[0]->num_lines; i++) {
		a = ems[0]->raw_em[i][1];
		b = ems[0]->raw_em[i][2];
		if (a == I && b == X) pix = ems[0]->raw_rel_p[i];
		else if (a == I && b == Y) piy = ems[0]->raw_rel_p[i];
		else if (a == I && b == Z) piz = ems[0]->raw_rel_p[i];
		else if (a == X && b == I) pxi = ems[0]->raw_rel_p[i];
		else if (a == X && b == X) pxx = ems[0]->raw_rel_p[i];
		else if (a == X && b == Y) pxy = ems[0]->raw_rel_p[i];
		else if (a == X && b == Z) pxz = ems[0]->raw_rel_p[i];
		else if (a == Y && b == I) pyi = ems[0]->raw_rel_p[i];
		else if (a == Y && b == X) pyx = ems[0]->raw_rel_p[i];
		else if (a == Y && b == Y) pyy = ems[0]->raw_rel_p[i];
		else if (a == Y && b == Z) pyz = ems[0]->raw_rel_p[i];
		else if (a == Z && b == I) pzi = ems[0]->raw_rel_p[i];
		else if (a == Z && b == X) pzx = ems[0]->raw_rel_p[i];
		else if (a == Z && b == Y) pzy = ems[0]->raw_rel_p[i];
		else if (a == Z && b == Z) pzz = ems[0]->raw_rel_p[i];
	}

	ppix = pix + piy + pzx + pzy;
	ppxi = pxi + pxz + pyi + pyz;
	ppxx = pxx + pxy + pyx + pyy;

	if (AVERAGE) {
		m = pe_av(ppix, ppxi, ppxx);
	}
	else {
		m = pe_max(ppix, ppxi, ppxx);
	}
	ps3x[2] = 15*ems[0]->scale*m/4;

	ppiz = piz + piy + pxz + pxy;
	ppzi = pzi + pzx + pyi + pyx;
	ppzz = pzz + pzy + pyz + pyy;

	if (AVERAGE) {
		m = pe_av(ppiz, ppzi, ppzz);
	}
	else {
		m = pe_max(ppiz, ppzi, ppzz);
	}
	ps3z[2] = 15*ems[0]->scale*m/4;

	ps3x[1] = ps3z[1] = 0;

	px = py = pz = 0;
	for (i=0; i<ems[5]->num_lines; i++) {
		a = ems[5]->raw_em[i][1];
		if (a == X) px = ems[5]->raw_rel_p[i];
		if (a == Y) py = ems[5]->raw_rel_p[i];
		if (a == Z) pz = ems[5]->raw_rel_p[i];
	}

	ps3x[1] += 2*ems[5]->scale*(px + py);
	ps3z[1] += 2*ems[5]->scale*(py + pz);

	px = py = pz = 0;
	for (i=0; i<ems[6]->num_lines; i++) {
		a = ems[6]->raw_em[i][1];
		if (a == X) px = ems[6]->raw_rel_p[i];
		if (a == Y) py = ems[6]->raw_rel_p[i];
		if (a == Z) pz = ems[6]->raw_rel_p[i];
	}

	ps3x[1] += ems[6]->scale*(px + py);
	ps3z[1] += ems[6]->scale*(py + pz);

	px = py = pz = 0;
	for (i=0; i<ems[7]->num_lines; i++) {
		a = ems[7]->raw_em[i][1];
		if (a == X) px = ems[7]->raw_rel_p[i];
		if (a == Y) py = ems[7]->raw_rel_p[i];
		if (a == Z) pz = ems[7]->raw_rel_p[i];
	}

	ps3x[1] += ems[7]->scale*(px + py);
	ps3z[1] += ems[7]->scale*(py + pz);

	// 1.5 scale up to error rate including all types of error
	// 1/4 scale down to error rate per identity gate
	ps3x[1] *= 1.5/4;
	ps3z[1] *= 1.5/4;

	ps3x[0] = ps3z[0] = 0;

	px = py = pz = 0;
	for (i=0; i<ems[2]->num_lines; i++) {
		a = ems[2]->raw_em[i][1];
		if (a == X) px = ems[2]->raw_rel_p[i];
		if (a == Y) py = ems[2]->raw_rel_p[i];
		if (a == Z) pz = ems[2]->raw_rel_p[i];
	}

	ps3x[0] += ems[2]->scale*px;
	ps3z[0] += ems[2]->scale*px;

	px = py = pz = 0;
	for (i=0; i<ems[3]->num_lines; i++) {
		a = ems[3]->raw_em[i][1];
		if (a == X) px = ems[3]->raw_rel_p[i];
		if (a == Y) py = ems[3]->raw_rel_p[i];
		if (a == Z) pz = ems[3]->raw_rel_p[i];
	}

	ps3x[0] += ems[3]->scale*px;
	ps3z[0] += ems[3]->scale*px;

	px = py = pz = 0;
	for (i=0; i<ems[1]->num_lines; i++) {
		a = ems[1]->raw_em[i][1];
		if (a == X) px = ems[1]->raw_rel_p[i];
		if (a == Y) py = ems[1]->raw_rel_p[i];
		if (a == Z) pz = ems[1]->raw_rel_p[i];
	}

	ps3z[0] += ems[1]->scale*(px + 2*py + pz);

	for (i=0; i<3; i++) {
		ps3x[i] *= p;
		ps3z[i] *= p;
	}
}

void pe_convert_json_to_ps(cJSON *json, double *ps) {
	int i;
	cJSON *ptr;

	for (i=0; i<NUM_EMS; i++) ps[i] = -1;

	ptr = json->child;

	while (ptr != NULL) {
		if (strcmp(ptr->string, "p_cnot") == 0) ps[0] = ptr->valuedouble;
		else if (strcmp(ptr->string, "p_H") == 0) ps[1] = ptr->valuedouble;
		else if (strcmp(ptr->string, "p_init_Z") == 0) ps[2] = ptr->valuedouble;
		else if (strcmp(ptr->string, "p_meas_Z") == 0) ps[3] = ptr->valuedouble;
		else if (strcmp(ptr->string, "p_iden_cnot") == 0) ps[4] = ptr->valuedouble;
		else if (strcmp(ptr->string, "p_iden_H") == 0) ps[5] = ptr->valuedouble;
		else if (strcmp(ptr->string, "p_iden_init_Z") == 0) ps[6] = ptr->valuedouble;
		else if (strcmp(ptr->string, "p_iden_meas_Z") == 0) ps[7] = ptr->valuedouble;
		else {
			fprintf(stderr, "Unknown error rate %s.\n", ptr->string);
			exit(EXIT_FAILURE);
		}
		ptr = ptr->next;
	}

	for (i=0; i<NUM_EMS; i++) {
		if (ps[i] == -1) {
			switch (i) {
				case 0:
					fprintf(stderr, "p_cnot missing.\n");
				case 1:
					fprintf(stderr, "p_H missing.\n");
				case 2:
					fprintf(stderr, "p_init_Z missing.\n");
				case 3:
					fprintf(stderr, "p_meas_Z missing.\n");
				case 4:
					fprintf(stderr, "p_iden_cnot missing.\n");
				case 5:
					fprintf(stderr, "p_iden_H missing.\n");
				case 6:
					fprintf(stderr, "p_iden_init_Z missing.\n");
				case 7:
					fprintf(stderr, "p_iden_meas_Z missing.\n");
			}
			exit(EXIT_FAILURE);
		}
	}
}

void pe_convert_ps_to_ps3s(double *ps, double *ps3x, double *ps3z) {
	// int i;

	ps3x[2] = ps3z[2] = ps[0];
	ps3x[1] = ps3z[1] = (2*ps[5] + ps[6] + ps[7])/4;
	ps3x[0] = ps[2] + ps[3];
	ps3z[0] = ps[2] + ps[3] + 4*ps[1]/3;

	// for (i=0; i<3; i++) printf("px[%d] = %g, pz[%d] = %g\n", i, ps3x[i], i, ps3z[i]);
}

STAT *pe_load_stat(FILE *in) {
	STAT *stat;
	int d;
	double p0, p1, p2, pX, pXlow, pXhigh, pZ, pZlow, pZhigh;

	if (fscanf(in, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", &d, &p0, &p1, &p2, &pX, &pXlow, &pXhigh, &pZ, &pZlow, &pZhigh) == EOF) return NULL;

	stat = (STAT *)my_malloc(sizeof(STAT));

	stat->d = d;
	stat->p0 = p0;
	stat->p1 = p1;
	stat->p2 = p2;
	stat->pX = pX;
	stat->pZ = pZ;

	return stat;
}

void pe_print_stat(STAT *stat) {
	printf("d: %d, p0: %g, p1: %g, p2: %g, pX: %g, pZ: %g\n", stat->d, stat->p0, stat->p1, stat->p2, stat->pX, stat->pZ);
	// printf("%d %g %g %g %g %g\n", stat->d, stat->p0, stat->p1, stat->p2, stat->pX, stat->pZ);
}

// Convert a code distance d and three error rates p[] of type X or Z
// into a logical error rate inferred from the database stat_ll
double pe_p_L_interpolate(int d, double *p, LL_NODE *stat_ll, int type) {
	int i;
	double p0_, p1_, p2_;
	double p0u, p0d, p1u, p1d, p2u, p2d, p0, p1, p2, del_p0, del_p1, del_p2;
	double pddd, pddu, pdud, pduu, pudd, pudu, puud, puuu;
	double pdd, pdu, pud, puu;
	double pd, pu, pc;
	LL_NODE *n;
	STAT *s, *sddd, *sddu, *sdud, *sduu, *sudd, *sudu, *suud, *suuu;

	if (d < 3 || d > 6) {
		fprintf(stderr, "Database only contains data for 3 <= d <= 6. Argument d = %d out of range.\n", d);
		exit(EXIT_FAILURE);
	}

	for (i=0; i<3; i++) {
		if (p[i] <= 0 || p[i] >= 0.5) {
			fprintf(stderr, "p[%d] not in range [0, 0.5). Error model invalid.\n", i);
			exit(EXIT_FAILURE);
		}
	}

	// for (i=0; i<3; i++) printf("p[%d] = %g\n", i, p[i]);

	// convert absolute input to relative
	p0_ = p[0]/p[2];
	p1_ = p[1]/p[2];
	p2_ = p[2];

	// printf("p0_ = %g, p1_ = %g, p2_ = %g\n", p0_, p1_, p2_);

	p0u = p1u = p2u = DBL_MAX;
	p0d = p1d = p2d = 0;

	n = stat_ll;
	while (n != NULL) {
		s = (STAT *)n->key;
		// pe_print_stat(s);
		p0 = s->p0;
		p1 = s->p1;
		p2 = s->p2;
		if (p0 >= p0_ && p0 < p0u) p0u = p0;
		if (p1 >= p1_ && p1 < p1u) p1u = p1;
		if (p2 >= p2_ && p2 < p2u) p2u = p2;
		if (p0 <= p0_ && p0 > p0d) p0d = p0;
		if (p1 <= p1_ && p1 > p1d) p1d = p1;
		if (p2 <= p2_ && p2 > p2d) p2d = p2;
		n = n->next;
	}

	// printf("p0d = %g, p0u = %g\np1d = %g, p1u = %g\np2d = %g, p2u = %g\n", p0d, p0u, p1d, p1u, p2d, p2u);

	if (p0u == DBL_MAX) {
		fprintf(stderr, "p[0]/p[2] = %g lies outside the database.\n", p0_);
		exit(EXIT_FAILURE);
	}

	if (p1u == DBL_MAX) {
		fprintf(stderr, "p[1]/p[2] = %g lies outside the database.\n", p1_);
		exit(EXIT_FAILURE);
	}

	if (p2u == DBL_MAX) {
		fprintf(stderr, "p[2] = %g lies outside the database.\n", p2_);
		exit(EXIT_FAILURE);
	}

	if (p0d == 0) p0d = p0u;
	if (p1d == 0) p1d = p1u;

	if (p2d == 0) {
		fprintf(stderr, "p[2] = %g lies outside the database.\n", p2_);
		exit(EXIT_FAILURE);
	}

	// printf("p0d = %g, p0u = %g\np1d = %g, p1u = %g\np2d = %g, p2u = %g\n", p0d, p0u, p1d, p1u, p2d, p2u);

	sddd = sddu = sdud = sduu = sudd = sudu = suud = suuu = NULL;

	n = stat_ll;
	while (n != NULL) {
		s = (STAT *)n->key;
		p0 = s->p0;
		p1 = s->p1;
		p2 = s->p2;
		if (d == s->d && p0d == p0 && p1d == p1 && p2d == p2) sddd = s;
		if (d == s->d && p0d == p0 && p1d == p1 && p2u == p2) sddu = s;
		if (d == s->d && p0d == p0 && p1u == p1 && p2d == p2) sdud = s;
		if (d == s->d && p0d == p0 && p1u == p1 && p2u == p2) sduu = s;
		if (d == s->d && p0u == p0 && p1d == p1 && p2d == p2) sudd = s;
		if (d == s->d && p0u == p0 && p1d == p1 && p2u == p2) sudu = s;
		if (d == s->d && p0u == p0 && p1u == p1 && p2d == p2) suud = s;
		if (d == s->d && p0u == p0 && p1u == p1 && p2u == p2) suuu = s;
		n = n->next;
	}

	if (sddd == NULL) {
		fprintf(stderr, "No data in database for d = %d, r0 = %g, r1 = %g, p2 = %g\n", d, p0d, p1d, p2d);
		exit(EXIT_FAILURE);
	}

	if (sddu == NULL) {
		fprintf(stderr, "No data in database for d = %d, r0 = %g, r1 = %g, p2 = %g\n", d, p0d, p1d, p2u);
		exit(EXIT_FAILURE);
	}

	if (sdud == NULL) {
		fprintf(stderr, "No data in database for d = %d, r0 = %g, r1 = %g, p2 = %g\n", d, p0d, p1u, p2d);
		exit(EXIT_FAILURE);
	}

	if (sduu == NULL) {
		fprintf(stderr, "No data in database for d = %d, r0 = %g, r1 = %g, p2 = %g\n", d, p0d, p1u, p2u);
		exit(EXIT_FAILURE);
	}

	if (sudd == NULL) {
		fprintf(stderr, "No data in database for d = %d, r0 = %g, r1 = %g, p2 = %g\n", d, p0u, p1d, p2d);
		exit(EXIT_FAILURE);
	}

	if (sudu == NULL) {
		fprintf(stderr, "No data in database for d = %d, r0 = %g, r1 = %g, p2 = %g\n", d, p0u, p1d, p2u);
		exit(EXIT_FAILURE);
	}

	if (suud == NULL) {
		fprintf(stderr, "No data in database for d = %d, r0 = %g, r1 = %g, p2 = %g\n", d, p0u, p1u, p2d);
		exit(EXIT_FAILURE);
	}

	if (suuu == NULL) {
		fprintf(stderr, "No data in database for d = %d, r0 = %g, r1 = %g, p2 = %g\n", d, p0u, p1u, p2u);
		exit(EXIT_FAILURE);
	}

	/*
	pe_print_stat(sddd);
	pe_print_stat(sddu);
	pe_print_stat(sdud);
	pe_print_stat(sduu);
	pe_print_stat(sudd);
	pe_print_stat(sudu);
	pe_print_stat(suud);
	pe_print_stat(suuu);
	*/

	if (type == X) {
		pddd = sddd->pX;
		pddu = sddu->pX;
		pdud = sdud->pX;
		pduu = sduu->pX;
		pudd = sudd->pX;
		pudu = sudu->pX;
		puud = suud->pX;
		puuu = suuu->pX;
	}
	else {
		pddd = sddd->pZ;
		pddu = sddu->pZ;
		pdud = sdud->pZ;
		pduu = sduu->pZ;
		pudd = sudd->pZ;
		pudu = sudu->pZ;
		puud = suud->pZ;
		puuu = suuu->pZ;
	}

	p0 = sddd->p0;
	p1 = sddd->p1;
	p2 = sddd->p2;

	del_p0 = sudd->p0 - p0;
	del_p1 = sdud->p1 - p1;
	del_p2 = sddu->p2 - p2;

	/*
	printf("del_p0 = %g\n", del_p0);
	printf("del_p1 = %g\n", del_p1);
	printf("del_p2 = %g\n", del_p2);
	*/

	if (del_p0 != 0) {
		pdd = (pudd - pddd)*(p0_ - p0)/del_p0 + pddd;
		pdu = (pudu - pddu)*(p0_ - p0)/del_p0 + pddu;
		pud = (puud - pdud)*(p0_ - p0)/del_p0 + pdud;
		puu = (puuu - pduu)*(p0_ - p0)/del_p0 + pduu;
	}
	else {
		pdd = pudd;
		pdu = pudu;
		pud = puud;
		puu = puuu;
	}

	// printf("pdd = %g, pdu = %g, pud = %g, puu = %g\n", pdd, pdu, pud, puu);

	if (del_p1 != 0) {
		pd = (pud - pdd)*(p1_ - p1)/del_p1 + pdd;
		pu = (puu - pdu)*(p1_ - p1)/del_p1 + pdu;
	}
	else {
		pd = pud;
		pu = puu;
	}

	if (del_p2 != 0) {
		pc = (pu - pd)*(p2_ - p2)/del_p2 + pd;
	}
	else {
		pc = pu;
	}

	return pc;
}

// Given two sets of three error rates and a database, calculate
// functional forms of the logical X and Z error rates for both even
// and odd code distances d.
FORM *pe_create_form(double *ps3x, double *ps3z, LL_NODE *stat_ll) {
	FORM *f;
	double p3x, p3z, p4x, p4z, p5x, p5z, p6x, p6z;

	f = (FORM *)my_malloc(sizeof(FORM));

	p3x = pe_p_L_interpolate(3, ps3x, stat_ll, X);
	p3z = pe_p_L_interpolate(3, ps3z, stat_ll, Z);
	p4x = pe_p_L_interpolate(4, ps3x, stat_ll, X);
	p4z = pe_p_L_interpolate(4, ps3z, stat_ll, Z);
	p5x = pe_p_L_interpolate(5, ps3x, stat_ll, X);
	p5z = pe_p_L_interpolate(5, ps3z, stat_ll, Z);
	p6x = pe_p_L_interpolate(6, ps3x, stat_ll, X);
	p6z = pe_p_L_interpolate(6, ps3z, stat_ll, Z);

	f->qx = p5x/p3x;
	f->qx2 = p6x/p4x;
	f->qz = p5z/p3z;
	f->qz2 = p6z/p4z;

	f->Ax = p3x/(f->qx*f->qx);
	f->Ax2 = p4x/(f->qx2*f->qx2);
	f->Az = p3z/(f->qz*f->qz);
	f->Az2 = p4z/(f->qz2*f->qz2);

	return f;
}

double pe_p_L_estimate(int d, FORM *f, int type) {
	int de;

	de = (d+1)/2;

	if (type == X) {
		if (d%2 == 0) return f->Ax2*pow(f->qx2, de);
		return f->Ax*pow(f->qx, de);
	}

	if (d%2 == 0) return f->Az2*pow(f->qz2, de);

	return f->Az*pow(f->qz, de);
}
