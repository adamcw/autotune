#define QUBIT_HT_SIZE 1000
#define DE_HT_FACTOR 100
#define STICK_HT_FACTOR 200
#define MAX_STRING_LENGTH 256

/** Surface Code Depolarizing Quantum Computer */
typedef struct sc_dp_qc SC_DP_QC; 
struct sc_dp_qc {
	DP_QC *dp_qc;
	int d, n;
	int **frame;
	QUBIT ***q_arr;
	SET *bdy_s1_pr, *bdy_s2_pr, *bdy_s1_du, *bdy_s2_du; // spatial boundaries
	SET *bdy_t1_pr, *bdy_t2_pr, *bdy_t1_du, *bdy_t2_du; // temporal boundaries
	SYNDROME ***syn_arr_pr, ***syn_arr_du;
	SET ****set_arr_pr, ****set_arr_du;
	double t1;
	int num_checks, last_X_check, last_Z_check, num_X_changes, num_Z_changes;
	BALL **boundaries;
	char *ems;
};

/** Arguments */
typedef struct args ARGS;
struct args {
	float p;
	int d;
	int s0;
	int s1;
	int t_check;
	int t_delete;
	int max_num_X;
	int max_num_Z;
	int verbose;
	int switch_time;
	int screen;
	int t_out;
	long int big_t;
	long int big_t_max;
	int t_check_scale;
	int boot;
	int boot_num_X;
	int boot_num_Z;
	int cap_time;
	char ems[MAX_STRING_LENGTH];
	FILE *out_raw;
	int primal;
};

// Arguments Functions
ARGS *init_args();
void load_args(ARGS *args, int argc, char **argv);

// Boundary Functions
BALL *get_boundary(int i, int j, long int big_t, int type, void *boundaries); 

// SC_DP_QC Functions
SC_DP_QC *create_sc_dp_qc(ARGS *args, RECIPE *recipe);
SET *create_boundary_set(BALL **boundaries, int id, int type, int i, int j, int t);
void create_initial_syndrome_and_set_arrays(SC_DP_QC *sc_dp_qc, int type, int imax, int jmax, SET *bdy_s1, SET *bdy_s2);
SC_DP_QC *copy_sc_dp_qc(SC_DP_QC *sc_dp_qc); // hopefully remove
void free_sc_dp_qc(SC_DP_QC *sc_dp_qc);
void free_sc_dp_qc_copy(SC_DP_QC *sc_dp_qc); // hopefully remove

// User Functions
void measure_stabilizers(SC_DP_QC *sc_dp_qc, long int big_t);
