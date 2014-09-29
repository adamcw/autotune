#ifndef POLYESTIMATE_H
#define POLYESTIMATE_H

#include <stdlib.h>
#include "../cJSON/cJSON.h"
#include "../qc/qc.h"

#define MAX_STRING_LENGTH 1024
#define NUM_EMS 8

typedef struct stat STAT;
typedef struct form FORM;
typedef LL_NODE DATABASE;

struct stat {
   int d;
   double p0, p1, p2, pX, pZ;
};

struct form {
   double Ax, Az, Ax2, Az2;
   double qx, qz, qx2, qz2;
};

// User functions
DATABASE *pe_load_database(char *fname);
void pe_free_database(DATABASE *db);
FORM *pe_calculate_formulae_from_p(double p, DATABASE *db);
FORM *pe_calculate_formulae_from_3ps(double *p3, DATABASE *db);
FORM *pe_calculate_formulae_from_6ps(double *p3x, double *p3z, DATABASE *db);
FORM *pe_calculate_formulae_from_8ps(double *p8, DATABASE *db);
FORM *pe_calculate_formulae_from_p_json(char *dir, DATABASE *db);
FORM *pe_calculate_formulae_from_gates_json(char *dir, DATABASE *db);
FORM *pe_calculate_formulae_from_ems(double p, char *dir, DATABASE *db);
void pe_free_formulae(FORM *f);
double pe_calculate_logical_error(int d, int type, FORM *f);

// Internal functions
double pe_convert_json_to_ems_and_p(cJSON *json, char *ems_dir);
void pe_convert_json_to_ems(cJSON *json, char *ems_dir);
int pe_convert_json_to_em(cJSON *json, char *ems_dir);
void pe_convert_json_to_em_line(cJSON *json, FILE *out);

double pe_max(double a, double b, double c);
double pe_av(double a, double b, double c);
void pe_convert_ems_to_ps3s(double p, ERROR_MODEL **ems, double *ps3x, double *ps3z);

void pe_convert_json_to_ps(cJSON *json, double *ps);
void pe_convert_ps_to_ps3s(double *ps, double *ps3x, double *ps3z);

STAT *pe_load_stat(FILE *in);
void pe_print_stat(STAT *stat);

double pe_p_L_interpolate(int d, double *p, LL_NODE *stat_ll, int type);
FORM *pe_create_form(double *ps3x, double *ps3z, LL_NODE *stat_ll);
double pe_p_L_estimate(int d, FORM *f, int type);

#endif
