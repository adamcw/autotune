#define QUBIT_HT_SIZE 1001
#define DE_HT_FACTOR 201
#define STICK_HT_FACTOR 201
#define MAX_STRING_LENGTH 256

/** Surface Code Depolarizing Quantum Computer */
typedef struct sc_dpl_qc SC_DPL_QC;
struct sc_dpl_qc {
	DPL_QC *dpl_qc;
	int d, n;
	int **frame;
	QUBIT ***q_arr, ***q_arr2;
	SET *bdy_s1_pr, *bdy_s2_pr, *bdy_s1_du, *bdy_s2_du; // spatial boundaries
	SET *bdy_t1_pr, *bdy_t2_pr, *bdy_t1_du, *bdy_t2_du; // temporal boundaries
	SYNDROME ***syn_arr_pr, ***syn_arr_du;
	SET ****set_arr_pr, ****set_arr_du;
	double t1;
	int num_checks, last_X_check, last_Z_check, num_X_changes, num_Z_changes;
	long int bdy_meas_pr, bdy_edges_pr, bdy_meas_du, bdy_edges_du;
	BALL **boundaries;
	char *ems;

	int flush_file;
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
	int t_delay;
    int max_num_X;
    int max_num_Z;
    int verbose;
    long int switch_time;
    int screen;
    int t_out;
	int track;
	int loss;
	double p_loss;
	double p_depolar;
	long int big_t;
    long int big_t_max;
	int t_check_scale;
    int boot;
    int boot_num_X;
    int boot_num_Z;
    int cap_time;
	char ems[MAX_STRING_LENGTH];
    FILE *out_raw;
	int flush_file;
};

// Argument Functions
ARGS *init_args();
void load_args(ARGS *args, int argc, char **argv);

// Recipe Functions
void generate_recipe(ARGS *args, RECIPE_ADV *recipe);
void simulate_recipe(ARGS *args, RECIPE_ADV *recipe);

// Boot Up Functions
void calculate_t_check(ARGS *args, RECIPE_ADV *recipe); 

// Boundary Functions
BALL *get_boundary(int i, int j, long int big_t, int type, void *boundaries); 

// Core Functions
SC_DPL_QC *create_sc_dpl_qc(ARGS *args, RECIPE_ADV *recipe, int has_loss);
SC_DPL_QC *copy_sc_dpl_qc(SC_DPL_QC *sc_dpl_qc);
void free_sc_dpl_qc(SC_DPL_QC *sc_dpl_qc);
void free_sc_dpl_qc_copy(SC_DPL_QC *sc_dpl_qc);

void measurement(SC_DPL_QC *sc_dpl_qc, int i, int j, QUBIT ***q_arr, SET ****set_arr, SET *set1, SET *set2);
void final_measurement(SC_DPL_QC *sc_dpl_qc, SET ****set_arr, SYNDROME ***syn_arr, QUBIT ***q_arr, int i, int j, int lay1, int lay2, int type, SET *bdy, int num_meas_left);

int adjust_meas(int mt, int i, int j, int n, int **frame);
void validate_set_arrays(int d, SET ****set_arr_pr, SET ****set_arr_du); 
void measure_stabilizers(SC_DPL_QC *sc_dpl_qc);
void advance_boundary(SC_DPL_QC *sc_dpl_qc, long int t);
void process_aug_edge(AUG_EDGE *ae, SC_DPL_QC *sc_dpl_qc, int type);
void process_aug_edges(MATCHING *m, SC_DPL_QC *sc_dpl_qc, int type);
void correct_mts(SC_DPL_QC *sc_dpl_qc);
void test_correct(SC_DPL_QC *sc_dpl_qc, FILE **out);

// Output Functions
void print_stats(SC_DPL_QC *sc_dpl_qc, FILE **out);
void print_all_es(int n, QUBIT ***q_arr);
void print_all_errors(int n, QUBIT ***q_arr);
