// surface code example
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>

#include "depolar/depolar.h"
#include "memory/memory.h"
#include "my_time/my_time.h"
#include "ex3.h"

typedef struct counter COUNTER;
typedef struct path PATH;
typedef struct selector SELECTOR;
typedef struct tuple TUPLE;
typedef struct pattern PATTERN;

struct counter {
	int digits;
	int base;
	int *x;
};

struct path {
	int length;
	STICK **s_arr;
};

struct selector {
	int n, k;
	int *indices;
};

struct tuple {
	long int big_t;
	int i;
};

struct pattern {
	int n;
	TUPLE *t_arr;
};

// Local functions
int compare_balls(BALL *a, BALL *b);
int compare_sticks(STICK **s1_ptr, STICK **s2_ptr);
COUNTER *create_counter(int digits, int base);
void free_counter(COUNTER *c);
int inc_counter(COUNTER *c);
void print_counter(COUNTER *c);
PATH *create_path(int length);
void free_path(PATH *path);
PATH *create_chosen_path(STICK *s, COUNTER *c);
DLL_NODE *get_ith_ht_elem(HT *ht, int i);
STICK *get_bdy_stick(BALL *b);
SELECTOR *create_selector(int n, int k);
void free_selector(SELECTOR *sel);
void reset_selector(SELECTOR *sel);
int increment_selector(SELECTOR *sel);
void print_selector(SELECTOR *sel);
PATTERN *create_pattern(PATH *path, SELECTOR *sel, int n, STICK **tc_array);
void free_pattern(PATTERN *pat);
void print_pattern(PATTERN *pat);

/**
 * \brief Example 1: Implements an abstract surface code using an 8-step cycle.   
 * 
 * \param[in] argc The number of command line arguments
 * \param[in] argv An array of the command line arguments 
 *
 * \return 0
 */
int main(int argc, char **argv) {
	ARGS *args;
	SC_DP_QC *sc_dp_qc;
	DP_QC *dp_qc;
	QC *qc;
	long int big_t, big_t_max;
	NEST *nest;
	int i, j; // loop variables
	int d; // code distance
	int n; // number of sticks in a repeating layer of the nest
	int nb; // number of balls connected to one boundary of a repeating layer of the nest
	int primal; // TRUE or FALSE depending on whether the primal or dual nest is being considered
	int tc; // central value of big_t
	int ms; // maximum number of sticks emanating from any ball
	int len, len_min, len_max; // minimum and maximum length opposing boundary connecting paths considered
	int np; // number of non-NULL paths
	CDLL_NODE *cdlln;
	STICK *s, ***all_sticks, **bdy_sticks;
	BALL *ball;
	COUNTER *c;
	PATH *path;
	SELECTOR *sel;
	PATTERN *pat;

	// Initialise, then load the args
	args = init_args();
	load_args(args, argc, argv);

	d = args->d;
	big_t_max = 4*((d+1)/2) + 1;
	tc = 4*((d+1)/2);
	primal = args->primal;

	sc_dp_qc = create_sc_dp_qc(args, NULL);
	dp_qc = sc_dp_qc->dp_qc;
	qc = dp_qc->qc;

	// Enable error tracking
	qc->track = TRUE;
	
	// Let QC know about the boundaries
	qc->get_boundary = get_boundary;
	qc->boundaries = (void *)sc_dp_qc->boundaries;

	fprintf(args->out_raw, 
		"p = %g\n"
		"d = %d\n"
		"big_t_max = %ld\n"
		"new t_check = %d\n"
		"t_delete = %d\n"
		"max_num_X = %d\n"
		"max_num_Z = %d\n"
		"verbose = %d\n"
		"s0 = %d;\n"
		"s1 = %d;\n"
		"-s0 %d -s1 %d\n"
		"%s\n",
		args->p, d, big_t_max, args->t_check, args->t_delete, 
		args->max_num_X, args->max_num_Z,
		args->verbose, qc->s0, qc->s1, qc->s0, qc->s1,
		primal ? "-primal" : "-dual");

	// Create nests of sufficient size
	big_t = 0;
	while (big_t <= big_t_max + 3) {
		measure_stabilizers(sc_dp_qc, big_t);
		qc_convert_nests(qc, FALSE);
		assert(qc->big_t == big_t+1);
		big_t++;
	}

	if (primal) nest = qc->nest_pr;
	else nest = qc->nest_du;
	
	// Count sticks in one repeating layer
	assert(d > 2);
	n = 0;
	cdlln = nest->stick_cdll->next;
	while (cdlln != nest->stick_cdll) {
		s = (STICK *)cdlln->key;
		if ((s->a->big_t == tc && s->b->big_t >= tc) || (s->a->big_t == tc+1 && s->b->big_t == tc)) n++;
		cdlln = cdlln->next;
	}
	printf("n: %d\n", n);

	// Store each layer of sticks in an array
	all_sticks = (STICK ***)my_2d_calloc(big_t_max, n, sizeof(STICK *));
	for (i=0; i<big_t_max; i++) {
		j = 0;
		cdlln = nest->stick_cdll->next;
		while (cdlln != nest->stick_cdll) {
			s = (STICK *)cdlln->key;
			if ((s->a->big_t == i && s->b->big_t >= i) || (s->a->big_t == i+1 && s->b->big_t == i)) {
				all_sticks[i][j] = s;
				j++;
			}
			cdlln = cdlln->next;
		}
	}

	// Sort stick pointers in each layer
	for (i=0; i<big_t_max; i++) {
		j = 0;
		s = all_sticks[i][j];
		while (s != NULL) {
			j++;
			if (j == n) break;
			s = all_sticks[i][j];
		}
		qsort(all_sticks[i], j, sizeof(STICK *), (int (*)(const void *, const void *))compare_sticks);
	}

	// Count boundary sticks in one boundary in one repeating layer
	nb = 0;
	cdlln = nest->stick_cdll->next;
	while (cdlln != nest->stick_cdll) {
		s = (STICK *)cdlln->key;
		if (s->a->big_t == tc && s->b->i == -1) nb++;
		cdlln = cdlln->next;
	}
	printf("nb: %d\n", nb);

	// Store tc layer of -1 boundary sticks in an array
	bdy_sticks = (STICK **)my_calloc(nb, sizeof(STICK *));
	i = 0;
	cdlln = nest->stick_cdll->next;
	while (cdlln != nest->stick_cdll) {
		s = (STICK *)cdlln->key;
		if (s->a->big_t == tc && s->b->i == -1) {
			bdy_sticks[i] = s;
			// qc_print_stick(s);
			i++;
		}
		cdlln = cdlln->next;
	}

	// Find maximum number of sticks emanating from any ball
	ms = 0;
	cdlln = nest->ball_cdll->next;
	while (cdlln != nest->ball_cdll) {
		ball = (BALL *)cdlln->key;
		if (ball->stick_ht->num_elem > ms) ms = ball->stick_ht->num_elem;
		cdlln = cdlln->next;
	}
	printf("ms: %d\n", ms);

	len_min = d;
	if (d%2 == 1) len_max = d+1;
	else len_max = d;

	np = 0;
	for (len = len_min; len <= len_max; len++) {
		for (i=0; i<nb; i++) {
			s = bdy_sticks[i];
			c = create_counter(len, ms);
			while (inc_counter(c)) {
				// print_counter(c);
				path = create_chosen_path(s, c);
				if (path != NULL) {
					np++;
					sel = create_selector(len, (d+1)/2);
					while (increment_selector(sel)) {
						pat = create_pattern(path, sel, n, all_sticks[tc]);
						// print_pattern(pat);
						free_pattern(pat);
					}
					free_selector(sel);
					free_path(path);
				}
			}
			free_counter(c);
		}
	}
	printf("np: %d\n", np);

	// Clean up memory
	free(bdy_sticks);
	my_2d_free(big_t_max, (void **)all_sticks);
	free_sc_dp_qc(sc_dp_qc);
	if (args->out_raw != stdout) { 
		fclose(args->out_raw);
	}
	free(args);

	return 0;
}

/**
 * \brief Creates and initialises an \ref args with default values
 * 
 * \return The initialised \ref args
 */
ARGS *init_args() {
	ARGS *args = (ARGS *)my_malloc(sizeof(ARGS)); 

	args->p = 0.01;
	args->d = 3;
	args->big_t = 0;
	args->big_t_max = LONG_MAX;
	args->s0 = 0;
	args->s1 = 0;
	args->t_check = 1;
	args->t_delete = 100;
	args->max_num_X = 10000;
	args->max_num_Z = 10000;
	args->verbose = 0;
	args->screen = 0;
	args->t_out = 0;
	args->t_check_scale = 10;
	args->boot = 1;
	args->boot_num_X = 50;
	args->boot_num_Z = 50;
	args->cap_time = 0;
	args->switch_time = 2;
	args->primal = TRUE;

	strcpy(args->ems, "../ems/");
	
	return args;
}

/**
 * \brief Loads the command lines args into the args struct
 * 
 * \param[out] args The args struct to be loaded into
 * \param[in] argc The number of arguments to read
 * \param[in] argv The array of arguments 
 */
void load_args(ARGS *args, int argc, char **argv) {
	for (int j = 1; j < argc; j++) {
		if (!strcmp(argv[j], "-h")) {
			printf("Usage:\n"
					"	ex1\n\n"
					"Options:\n"
					"-p P                          The probability of a random error. [default: %f]\n"
					"-d D                          The distance of the surface code. [default: %d]\n"
					"-big_t_max BIG_T_MAX          The maximum value of big_t. [default: %ld]\n"
					"-s0 S0                        The first random seed  (Used to redo past simulations). [default: %d]\n"
					"-s1 S1                        The second random seed (Used to redo past simulations). [default: %d]\n"

					"-t_check T_CHECK              The number of time steps between each round of perfect stabilizer measurements. [default: %d]\n"
					"-t_delete T_DELETE            The number of time steps in the past to keep during the computation. [default: %d]\n"
					"-max_num_X MAX_NUM_X          The maximum number of X changes before ending the simulation. [default: %d]\n"
					"-max_num_Z MAX_NUM_Z          The maximum number of Z changes before ending the simulation. [default: %d]\n"
					"-verbose VERBOSE              Verbose output. [default: %d]\n"
					"-switch_time SWITCH_TIME      When to switch to the repeating layer in the boot up phase [default: %d]\n"
					"-t_check_scale T_CHECK_SCALE  A scaling factor applied to t_check to ensure it is reasonable. [default: %d]\n"
					"-boot BOOT	                   Whether or not to have a boot up phase to calculate t_check [default: %d]\n"
					"-boot_num_X BOOT_NUM_X        The number of X changes to wait for during the boot up phase. [default: %d]\n"
					"-boot_num_Z BOOT_NUM_Z        The number of Z changes to wait for during the boot up phase. [default: %d]\n"
					"-cap_time CAP_TIME            [default: %d]\n"
					"-ems EMS                      The folder containing the error models to use. [default: %s]\n"
					"-screen SCREEN                Whether or not to output to the screen. [default: %d]\n"
					"-t_out T_OUT                  How often to output the current status. Useful at low error rates. [default: %d]\n" 
					"-primal                       Process the primal lattice. [default: TRUE]\n"
					"-dual                         Process the dual lattice. [default: FALSE]\n",

			args->p, args->d, args->big_t_max, args->s0, args->s1, 
			args->t_check, args->t_delete, args->max_num_X, args->max_num_Z, args->verbose, 
			args->switch_time, args->t_check_scale, args->boot, args->boot_num_X, args->boot_num_Z, 
			args->cap_time, args->ems, args->screen, args->t_out);

			exit(0);
		}
		else if (!strcmp(argv[j], "-p")) {
			args->p = atof(argv[++j]);
		}
		else if (!strcmp(argv[j], "-d")) {
			args->d = atoi(argv[++j]);
		}
		else if (!strcmp(argv[j], "-big_t_max")) {
			args->big_t_max = atol(argv[++j]);
		}
		else if (!strcmp(argv[j], "-s0")) {
			args->s0 = atoi(argv[++j]);
		}
		else if (!strcmp(argv[j], "-s1")) {
			args->s1 = atoi(argv[++j]);
		}
		else if (!strcmp(argv[j], "-t_check")) {
			args->t_check = atoi(argv[++j]);
		}
		else if (!strcmp(argv[j], "-t_delete")) {
			args->t_delete = atoi(argv[++j]);
		}
		else if (!strcmp(argv[j], "-max_num_X")) {
			args->max_num_X = atoi(argv[++j]);
		}
		else if (!strcmp(argv[j], "-max_num_Z")) {
			args->max_num_Z = atoi(argv[++j]);
		}
		else if (!strcmp(argv[j], "-verbose")) {
			args->verbose = 1;
		}	
		else if (!strcmp(argv[j], "-switch_time")) {
			args->switch_time = atoi(argv[++j]);
		}	
		else if (!strcmp(argv[j], "-t_check_scale")) {
			args->t_check_scale = atoi(argv[++j]);
		}	
		else if (!strcmp(argv[j], "-boot")) {
			args->boot = atoi(argv[++j]);
		}	
		else if (!strcmp(argv[j], "-boot_num_X")) {
			args->boot_num_X = atoi(argv[++j]);
		}	
		else if (!strcmp(argv[j], "-boot_num_Z")) {
			args->boot_num_Z = atoi(argv[++j]);
		}	
		else if (!strcmp(argv[j], "-cap_time")) {
			args->cap_time = atoi(argv[++j]);
		}	
		else if (!strcmp(argv[j], "-ems")) {
			++j;

			// Prevent buffer overflow on long directories.
			if (strlen(argv[j]) < MAX_STRING_LENGTH) {
				strcpy(args->ems, argv[j]);
			} else {
				printf("-ems path can be no longer than %d characters\n", MAX_STRING_LENGTH);
				exit(1);
			}
		}
		else if (!strcmp(argv[j], "-screen")) {
			args->screen = 1;
		}		
		else if (!strcmp(argv[j], "-t_out")) { // a flag to preiodically output to screen, helpful at low error rates 
			args->t_out = atoi(argv[++j]);
		}
		else if (!strcmp(argv[j], "-primal")) {
			args->primal = TRUE;
		}		
		else if (!strcmp(argv[j], "-dual")) {
			args->primal = FALSE;
		}		
		else {
			printf("Unknown switch: %s\n", argv[j]);
			exit(0);
		}
	}

	// Set the output file
	args->out_raw = (args->screen) ? stdout : (FILE *)fopen("out_raw", "w");

	// Ensure that the folder for the error models ends in a trailing slash
	if (args->ems[(int)strlen(args->ems)-1] != '/') {
		strcat(args->ems, "/");
	}
}

/**
 * \brief Gets the boundary corresponding to given i, j, big_t coordinates and
 * type.
 *
 * \param[in] i The i coordinate of the boundary
 * \param[in] j The j coordinate of the boundary
 * \param[in] big_t The big_t coordinate of the boundary (this should be
 * LONG_MAX)
 * \param[in] type The type of the boundary (PRIMAL_BOUNDARY, DUAL_BOUNDARY)
 * \param[in] boundaries An array of boundaries to search for the boundary in
 *
 * \return The boundary if found, NULL otherwise.
 */
BALL *get_boundary(int i, int j, long int big_t, int type, void *boundaries) {
	BALL **bdys;

	bdys = (BALL **)boundaries;

	// Boundaries do not have a big_t set, they are initialised to LONG_MAX.
	if (big_t != LONG_MAX) {
		return NULL;
	}

	// This function essentially takes the coordinates of a boundary and then
	// returns the boundary. This could be implemented as a map of some sort, in
	// this instance, given the limited number of boundaries, it is easier,
	// cleaner and faster to just nest if statements in this way.
	if (i == -1) {
		if (j == 0) {
			if (type == PRIMAL_BOUNDARY) {
				return bdys[0]; // -1, 0, 0 PRIMAL
			} else {
				return bdys[4]; // -1, 0, 0 DUAL
			}
		} else {
			if (type == PRIMAL_BOUNDARY) {
				return bdys[2]; // -1, 1, 0 PRIMAL
			} else {
				return bdys[6]; // -1, 1, 0 DUAL
			}
		}
	} else {
		if (j == 0) {
			if (type == PRIMAL_BOUNDARY) {
				return bdys[1]; // -2, 0, 0 PRIMAL
			} else {
				return bdys[5]; // -2, 0, 0 DUAL
			}
		} else {
			if (type == PRIMAL_BOUNDARY) {
				return bdys[3]; // -2, 1, 0, PRIMAL
			} else {
				return bdys[7]; // -2, 1, 0, DUAL
			}
		}
	}

	// There is no boundary at this coordinate
	return NULL;
}

/**
 * \brief Creates a new \ref sc_dp_qc
 * 
 * \param[in] args The args structure as filled from the command line
 * \param[in] recipe The recipe the \ref sc_dp_qc will use for computation 
 */
SC_DP_QC *create_sc_dp_qc(ARGS *args, RECIPE *recipe) {
	SC_DP_QC *sc_dp_qc;
	int i, j, d, n;
	QUBIT ***q_arr;

	sc_dp_qc = (SC_DP_QC *)my_malloc(sizeof(SC_DP_QC));

	sc_dp_qc->d = d = args->d;
	sc_dp_qc->n = n = 2*d - 1;
	sc_dp_qc->ems = args->ems;

	sc_dp_qc->t1 = double_time();

	sc_dp_qc->num_checks = 0;
	sc_dp_qc->last_X_check = 0;
	sc_dp_qc->last_Z_check = 0;
	sc_dp_qc->num_X_changes = 0;
	sc_dp_qc->num_Z_changes = 0;

	// Create a new dp_qc for the sc_dp_qc
	sc_dp_qc->dp_qc = dp_create_dp_qc(args->s0, args->s1, DE_HT_FACTOR * n * n, 
		STICK_HT_FACTOR * d * d, args->p, args->t_delete, recipe, args->ems);
	
	// Create an n by n Pauli frame
	sc_dp_qc->frame = (int **)my_2d_calloc(n, n, sizeof(int));

	// Create an n by n array of qubits
	sc_dp_qc->q_arr = q_arr = (QUBIT ***)my_2d_calloc(n, n, sizeof(QUBIT *));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			q_arr[i][j] = qc_create_and_insert_qubit(sc_dp_qc->dp_qc->qc, i, j, 0, QUBIT_HT_SIZE);
		}
	}
	
	// Allocate space for boundaries, and then create and insert them
	// The location of boundaries is irrelevant, but a sensible choice is good for debugging
	sc_dp_qc->boundaries = (BALL **)my_malloc(8 * sizeof(BALL *));
	sc_dp_qc->bdy_s1_pr = create_boundary_set(sc_dp_qc->boundaries, 0, PRIMAL_BOUNDARY, -1, 0, 0);
	sc_dp_qc->bdy_s2_pr = create_boundary_set(sc_dp_qc->boundaries, 1, PRIMAL_BOUNDARY, -2, 0, 0);
	sc_dp_qc->bdy_t1_pr = create_boundary_set(sc_dp_qc->boundaries, 2, PRIMAL_BOUNDARY, -1, 1, 0);
	sc_dp_qc->bdy_t2_pr = create_boundary_set(sc_dp_qc->boundaries, 3, PRIMAL_BOUNDARY, -2, 1, 0);
	sc_dp_qc->bdy_s1_du = create_boundary_set(sc_dp_qc->boundaries, 4, DUAL_BOUNDARY, -1, 0, 0);
	sc_dp_qc->bdy_s2_du = create_boundary_set(sc_dp_qc->boundaries, 5, DUAL_BOUNDARY, -2, 0, 0);
	sc_dp_qc->bdy_t1_du = create_boundary_set(sc_dp_qc->boundaries, 6, DUAL_BOUNDARY, -1, 1, 0);
	sc_dp_qc->bdy_t2_du = create_boundary_set(sc_dp_qc->boundaries, 7, DUAL_BOUNDARY, -2, 1, 0);

	// Create the Primal and Dual syndrome and set arrays 
	create_initial_syndrome_and_set_arrays(sc_dp_qc, PRIMAL, d, d-1, sc_dp_qc->bdy_s1_pr, sc_dp_qc->bdy_s2_pr);
	create_initial_syndrome_and_set_arrays(sc_dp_qc, DUAL, d-1, d, sc_dp_qc->bdy_s1_du, sc_dp_qc->bdy_s2_du);
	
	return sc_dp_qc;
}

/**
 * \brief Creates a boundary \ref set and inserts its \ref ball it into an array of boundaries
 * 
 * \param[out] boundaries The array of \ref ball%s to insert the new \ref set%'s
 * \ref ball into
 * \param[in] id The index the the boundaries array to insert into
 * \param[in] i The i coordinate of the boundary
 * \param[in] j The j coordinate of the boundary
 * \param[in] t The t coordinate of the boundary
 *
 * \return The newly created \ref set for the boundary
 */
SET *create_boundary_set(BALL **boundaries, int id, int type, int i, int j, int t) {
	SET *set;

	set = qc_create_set(type, i, j, t, NULL);
	boundaries[id] = set->ball;

	return set;
}

/**
 * \brief Creates the initial; syndrome and set arrays
 * 
 * Creates an array of syndromes with a new syndrome created at each point. It
 * then creates initial sets associated with the first round of measurements.
 *
 * \param[in] scp_dp_qc The \ref sc_dp_qc to build the \ref syndrome and \ref set arrays in
 * \param[in] type The type of arrays to make (PRIMAL, DUAL)
 * \param[in] imax The size of the i dimension 
 * \param[in] jmax The size of the j dimension 
 * \param[in] bdy_s1 The first space boundary to connect to
 * \param[in] bdy_s2 The second space boundary to connect to 
 */
void create_initial_syndrome_and_set_arrays(SC_DP_QC *sc_dp_qc, int type, int imax, int jmax, SET *bdy_s1, SET *bdy_s2) {
	int seti, setj;
	int is, js, bs, bmax;
	SYNDROME ***syn_arr;
	SET ****set_arr;

	syn_arr = (SYNDROME ***)my_2d_calloc(imax, jmax, sizeof(SYNDROME *));
	set_arr = (SET ****)my_3d_calloc(imax, jmax, 2, sizeof(SET *));

	for (is = 0; is < imax; is++) {
		for (js = 0; js < jmax; js++) {
			// Create and insert a new syndrome into the syndrome array
			syn_arr[is][js] = qc_create_syndrome();
			qc_insert_syndrome(sc_dp_qc->dp_qc->qc, syn_arr[is][js]);

			// Calculate the coordinate relative to the full frame
			if (type == PRIMAL) {
				seti = (is << 1);
				setj = (js << 1) + 1;
				bs = js;
				bmax = jmax;
			}
			else {
				seti = (is << 1) + 1;
				setj = (js << 1);
				bs = is;
				bmax = imax;
			}

			if (bs == 0) {
				// Create a set with the first boundary
				set_arr[is][js][1] = qc_create_set(type, seti, setj, 1, bdy_s1);
				set_arr[is][js][0] = qc_create_set(type, seti, setj, 2, bdy_s1);
			}
			else if (bs == bmax-1) {
				// Create a set with the second boundary
				set_arr[is][js][1] = qc_create_set(type, seti, setj, 1, bdy_s2);
				set_arr[is][js][0] = qc_create_set(type, seti, setj, 2, bdy_s2);
			}
			else {
				// Create a set with no boundary
				set_arr[is][js][1] = qc_create_set(type, seti, setj, 1, NULL);
				set_arr[is][js][0] = qc_create_set(type, seti, setj, 2, NULL);
			}

			// Insert the new set into the set array, and associate the set with a syndrome
			qc_insert_set(sc_dp_qc->dp_qc->qc, set_arr[is][js][1]);
			qc_associate_syndrome(set_arr[is][js][1], syn_arr[is][js]);
			qc_insert_set(sc_dp_qc->dp_qc->qc, set_arr[is][js][0]);
		}
	}
	
	// Store the set and syndrome arrays in the sc_dp_qc
	if (type == PRIMAL) {
		sc_dp_qc->syn_arr_pr = syn_arr;
		sc_dp_qc->set_arr_pr = set_arr;
	} 
	else {
		sc_dp_qc->syn_arr_du = syn_arr;
		sc_dp_qc->set_arr_du = set_arr;
	}
}


/**
 * \brief Copies a \ref sc_dp_qc
 * 
 * \param[in] sc_dp_qc The \ref sc_dp_qc to be copied
 *
 * \return The newly copied \ref sc_dp_qc
 */
SC_DP_QC *copy_sc_dp_qc(SC_DP_QC *sc_dp_qc) {
	SC_DP_QC *sc_dp_qc2;
	int i, j;
	int is, js;
	int d, n;

	d = sc_dp_qc->d;
	n = sc_dp_qc->n;
	
	sc_dp_qc2 = (SC_DP_QC *)my_malloc(sizeof(SC_DP_QC));
	
	// Copy values from sc_dp_qc to sc_dp_qc2
	*sc_dp_qc2 = *sc_dp_qc;

	// Copy the dp_qc
	sc_dp_qc2->dp_qc = dp_copy_dp_qc(sc_dp_qc->dp_qc);

	// Copy the frame
	sc_dp_qc2->frame = (int **)my_2d_calloc(n, n, sizeof(int));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			sc_dp_qc2->frame[i][j] = sc_dp_qc->frame[i][j];
		}
	}

	// Copy the qubits
	sc_dp_qc2->q_arr = (QUBIT ***)my_2d_calloc(n, n, sizeof(QUBIT *));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			sc_dp_qc2->q_arr[i][j] = sc_dp_qc->q_arr[i][j]->copy;
		}
	}

	// Copy the PRIMAL syndrome and set arrays 
	sc_dp_qc2->syn_arr_pr = (SYNDROME ***)my_2d_calloc(d, d-1, sizeof(SYNDROME *));
	sc_dp_qc2->set_arr_pr = (SET ****)my_3d_calloc(d, d-1, 2, sizeof(SET *));
	for (is = 0; is < d; is++) {
		for (js = 0; js < d-1; js++) {
			sc_dp_qc2->syn_arr_pr[is][js] = sc_dp_qc->syn_arr_pr[is][js]->copy;
			sc_dp_qc2->set_arr_pr[is][js][0] = sc_dp_qc->set_arr_pr[is][js][0]->copy;
			sc_dp_qc2->set_arr_pr[is][js][1] = sc_dp_qc->set_arr_pr[is][js][1]->copy;
		}
	}
	
	// Copy the DUAL syndrom and set arrays
	sc_dp_qc2->syn_arr_du = (SYNDROME ***)my_2d_calloc(d-1, d, sizeof(SYNDROME *));
	sc_dp_qc2->set_arr_du = (SET ****)my_3d_calloc(d-1, d, 2, sizeof(SET *));
	for (is = 0; is < d-1; is++) {
		for (js = 0; js < d; js++) {
			sc_dp_qc2->syn_arr_du[is][js] = sc_dp_qc->syn_arr_du[is][js]->copy;
			sc_dp_qc2->set_arr_du[is][js][0] = sc_dp_qc->set_arr_du[is][js][0]->copy;
			sc_dp_qc2->set_arr_du[is][js][1] = sc_dp_qc->set_arr_du[is][js][1]->copy;
		}
	}

	return sc_dp_qc2;
}

/**
 * \brief Frees a \ref sc_dp_qc
 * 
 * \param[in] sc_dp_qc The \ref sc_dp_qc to be freed
 */
void free_sc_dp_qc(SC_DP_QC *sc_dp_qc) {
	int is, js;
	int d, n;

	d = sc_dp_qc->d;
	n = sc_dp_qc->n;

	// Free the \ref dp_qc
	dp_free_dp_qc(sc_dp_qc->dp_qc);
	
	// Free the Pauli frame
	my_2d_free(n, (void **)sc_dp_qc->frame);

	// Free the array of qubits (the qubits themselves are owned by the qc)
	my_2d_free(n, (void **)sc_dp_qc->q_arr);
	
	// Free the primal syndrome array and the syndromes
	for (is = 0; is < d; is++) {
		for (js = 0; js < d-1; js++) {
			qc_free_syndrome(sc_dp_qc->syn_arr_pr[is][js]);
		}
	}
	my_2d_free(d, (void **)sc_dp_qc->syn_arr_pr);
	
	// Free the dual syndrome array and the syndromes
	for (is = 0; is < d-1; is++) {
		for (js = 0; js < d; js++) {
			qc_free_syndrome(sc_dp_qc->syn_arr_du[is][js]);
		}
	}
	my_2d_free(d-1, (void **)sc_dp_qc->syn_arr_du);

	// Free the primal and dual set arrays
	my_3d_free(d, d-1, (void ***)sc_dp_qc->set_arr_pr);
	my_3d_free(d-1, d, (void ***)sc_dp_qc->set_arr_du);

	// Free the sets associated with the boundaries
	qc_free_set(sc_dp_qc->bdy_s1_pr);
	qc_free_set(sc_dp_qc->bdy_s2_pr);
	qc_free_set(sc_dp_qc->bdy_t1_pr);
	qc_free_set(sc_dp_qc->bdy_t2_pr);
	qc_free_set(sc_dp_qc->bdy_s1_du);
	qc_free_set(sc_dp_qc->bdy_s2_du);
	qc_free_set(sc_dp_qc->bdy_t1_du);
	qc_free_set(sc_dp_qc->bdy_t2_du);

	// Free the boundaries array
	free(sc_dp_qc->boundaries);

	free(sc_dp_qc);
}

/**
 * \brief Frees a copied \ref sc_dp_qc
 * 
 * \param[in] sc_dp_qc The \ref sc_dp_qc copy to be freed
 */
void free_sc_dp_qc_copy(SC_DP_QC *sc_dp_qc) {
	/*
	 * \remark No need to free the boundaries, these are not copied
	 */
	int is, js, d, n;

	d = sc_dp_qc->d;
	n = sc_dp_qc->n;

	// Free the copied dp_qc
	dp_free_dp_qc_copy(sc_dp_qc->dp_qc);
	
	// Free the Pauli frame and the qubit array (The qubits are stored in qc)
	my_2d_free(n, (void **)sc_dp_qc->frame);
	my_2d_free(n, (void **)sc_dp_qc->q_arr);
	
	// Free the primal syndrome array and the syndromes
	for (is = 0; is < d; is++) {
		for (js = 0; js < d-1; js++) {
			qc_free_syndrome(sc_dp_qc->syn_arr_pr[is][js]);
		}
	}
	my_2d_free(d, (void **)sc_dp_qc->syn_arr_pr);
	
	// Free the dual syndrome array and the syndromes
	for (is = 0; is < d-1; is++) {
		for (js = 0; js < d; js++) {
			qc_free_syndrome(sc_dp_qc->syn_arr_du[is][js]);
		}
	}
	my_2d_free(d-1, (void **)sc_dp_qc->syn_arr_du);

	// Free the primal and dual set arrays
	my_3d_free(d, d-1, (void ***)sc_dp_qc->set_arr_pr);
	my_3d_free(d-1, d, (void ***)sc_dp_qc->set_arr_du);

	free(sc_dp_qc);
}

void measure_stabilizers(SC_DP_QC *sc_dp_qc, long int big_t) {
	int i, j, is, js, lay1, lay2, n;
	DP_QC *dp_qc;
	QC *qc;
	QUBIT ***q_arr;
	SET *bdy_s1_pr, *bdy_s2_pr; // primal boundaries
	SET *bdy_s1_du, *bdy_s2_du; // dual boundaries
	SYNDROME ***syn_arr_pr, ***syn_arr_du;
	SET ****set_arr_pr, ****set_arr_du;

	n = sc_dp_qc->n;

	dp_qc = sc_dp_qc->dp_qc;
	qc = dp_qc->qc;

	q_arr = sc_dp_qc->q_arr;

	bdy_s1_pr = sc_dp_qc->bdy_s1_pr;
	bdy_s2_pr = sc_dp_qc->bdy_s2_pr;

	bdy_s1_du = sc_dp_qc->bdy_s1_du;
	bdy_s2_du = sc_dp_qc->bdy_s2_du;

	syn_arr_pr = sc_dp_qc->syn_arr_pr;
	syn_arr_du = sc_dp_qc->syn_arr_du;

	set_arr_pr = sc_dp_qc->set_arr_pr;
	set_arr_du = sc_dp_qc->set_arr_du;

	// Initialise all syndrome qubits	
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			if (i%2 == 0 && j%2 == 1) {
				dp_init_Z(dp_qc, q_arr[i][j]);
			}
			else if (i%2 == 1 && j%2 == 0) {
				dp_dead_H(dp_qc, q_arr[i][j]);
			}
			else {
				dp_iden_init_Z(dp_qc, q_arr[i][j]);
			}
		}
	}

	// change basis of initialization of X-syndrome qubits
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			if (i%2 == 0 && j%2 == 1) {
				dp_H(dp_qc, q_arr[i][j]);
			}
			else if (i%2 == 1 && j%2 == 0) {
				dp_init_Z(dp_qc, q_arr[i][j]);
			}
			else {
				dp_iden_H(dp_qc, q_arr[i][j]);
			}
		}
	}
	
	// North
	for (j=1; j<n; j+=2) {
		dp_iden_cnot(dp_qc, q_arr[0][j]);
	}
	for (i=0; i<n-1; i++) {
		for (j=0; j<n; j+=2) {
			dp_cnot(dp_qc, q_arr[i][j], q_arr[i+1][j]);
		}
		i++;
		for (j=1; j<n; j+=2) {
			dp_cnot(dp_qc, q_arr[i+1][j], q_arr[i][j]);
		}
	}
	for (j=0; j<n; j+=2) {
		dp_iden_cnot(dp_qc, q_arr[i][j]);
	}
	
	// West
	for (i=0; i<n-1; i++) {
		for (j=0; j<n-1; j+=2) {
			dp_cnot(dp_qc, q_arr[i][j+1], q_arr[i][j]);
		}
		dp_iden_cnot(dp_qc, q_arr[i][j]);
		i++;
		dp_iden_cnot(dp_qc, q_arr[i][0]);
		for (j=1; j<n; j+=2) {
			dp_cnot(dp_qc, q_arr[i][j], q_arr[i][j+1]);
		}
	}
	for (j=0; j<n-1; j+=2) {
		dp_cnot(dp_qc, q_arr[i][j+1], q_arr[i][j]);
	}
	dp_iden_cnot(dp_qc, q_arr[i][j]);
	
	// East
	for (i=0; i<n-1; i++) {
		dp_iden_cnot(dp_qc, q_arr[i][0]);
		for (j=1; j<n; j+=2) {
			dp_cnot(dp_qc, q_arr[i][j], q_arr[i][j+1]);
		}
		i++;
		for (j=0; j<n-1; j+=2) {
			dp_cnot(dp_qc, q_arr[i][j+1], q_arr[i][j]);
		}
		dp_iden_cnot(dp_qc, q_arr[i][j]);
	}
	dp_iden_cnot(dp_qc, q_arr[i][0]);
	for (j=1; j<n; j+=2) {
		dp_cnot(dp_qc, q_arr[i][j], q_arr[i][j+1]);
	}

	// South
	for (j=0; j<n; j+=2) {
		dp_iden_cnot(dp_qc, q_arr[0][j]);
	}
	for (i=0; i<n-1; i++) {
		for (j=1; j<n; j+=2) {
			dp_cnot(dp_qc, q_arr[i][j], q_arr[i+1][j]);
		}
		i++;
		for (j=0; j<n; j+=2) {
			dp_cnot(dp_qc, q_arr[i+1][j], q_arr[i][j]);
		}
	}
	for (j=1; j<n; j+=2) {
		dp_iden_cnot(dp_qc, q_arr[i][j]);
	}
	
	lay1 = (big_t+1)%2;
	lay2 = big_t%2;

	// change basis of measurement of X-syndrome qubits
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			if (i%2 == 0 && j%2 == 1) {
				dp_H(dp_qc, q_arr[i][j]);
			}
			else if (i%2 == 1 && j%2 == 0) {
				is = i>>1;
				js = j>>1;

				// Z-syndrome
				dp_meas_Z(dp_qc, q_arr[i][j], set_arr_du[is][js][lay1], set_arr_du[is][js][lay2]);
				qc_unassociate_syndrome(syn_arr_du[is][js]);
				qc_associate_syndrome(set_arr_du[is][js][lay2], syn_arr_du[is][js]);
				if (i == 1) {
					set_arr_du[is][js][lay1] = qc_create_set(DUAL, i, j, 2, bdy_s1_du);
				}
				else if (i == n-2) {
					set_arr_du[is][js][lay1] = qc_create_set(DUAL, i, j, 2, bdy_s2_du);
				}
				else {
					set_arr_du[is][js][lay1] = qc_create_set(DUAL, i, j, 2, NULL);
				}						
				qc_insert_set(qc, set_arr_du[is][js][lay1]);
			}
			else {
				dp_iden_H(dp_qc, q_arr[i][j]);
			}
		}
	}
		
	// Identity and Measure Z
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			if (i%2 == 0 && j%2 == 1) {
				is = i>>1;
				js = j>>1;
				
				// X-syndrome
				dp_meas_Z(dp_qc, q_arr[i][j], set_arr_pr[is][js][lay1], set_arr_pr[is][js][lay2]);
				qc_unassociate_syndrome(syn_arr_pr[is][js]);
				qc_associate_syndrome(set_arr_pr[is][js][lay2], syn_arr_pr[is][js]);
				if (j == 1) {
					set_arr_pr[is][js][lay1] = qc_create_set(PRIMAL, i, j, 2, bdy_s1_pr);
				}
				else if (j == n-2) {
					set_arr_pr[is][js][lay1] = qc_create_set(PRIMAL, i, j, 2, bdy_s2_pr);
				}
				else {
					set_arr_pr[is][js][lay1] = qc_create_set(PRIMAL, i, j, 2, NULL);
				}						
				qc_insert_set(qc, set_arr_pr[is][js][lay1]);
			}
			else if (i%2 == 1 && j%2 == 0) {
				dp_dead_H(dp_qc, q_arr[i][j]);
			}
			else {
				dp_iden_meas_Z(dp_qc, q_arr[i][j]);
			}
		}
	}

	// printf("%ld %ld\n\n", qc->big_t, big_t);
	assert(qc->big_t == big_t+1);
}

int compare_balls(BALL *a, BALL *b) {
	assert(a != NULL);
	assert(b != NULL);

	if (a->i < b->i) return -1;
	if (a->i > b->i) return 1;
	if (a->j < b->j) return -1;
	if (a->j > b->j) return 1;
	if (a->big_t < b->big_t) return -1;
	if (a->big_t > b->big_t) return 1;

	return 0;
}

int compare_sticks(STICK **s1_ptr, STICK **s2_ptr) {
	int c;
	STICK *s1, *s2;
	BALL *a1, *b1, *a2, *b2;

	s1 = *s1_ptr;
	s2 = *s2_ptr;

	assert(s1 != NULL);
	assert(s2 != NULL);

	assert(compare_balls(s1->a, s1->b) != 0);
	assert(compare_balls(s2->a, s2->b) != 0);

	if (compare_balls(s1->a, s1->b) < 0) {
		a1 = s1->a;
		b1 = s1->b;
	}
	else {
		a1 = s1->b;
		b1 = s1->a;
	}

	if (compare_balls(s2->a, s2->b) < 0) {
		a2 = s2->a;
		b2 = s2->b;
	}
	else {
		a2 = s2->b;
		b2 = s2->a;
	}

	c = compare_balls(a1, a2);
	if (c != 0) return c;

	return compare_balls(b1, b2);
}

COUNTER *create_counter(int digits, int base) {
	COUNTER *c;

	c = (COUNTER *)my_malloc(sizeof(COUNTER));

	c->digits = digits;
	c->base = base;
	c->x = (int *)my_calloc(digits, sizeof(int));

	return c;
}

void free_counter(COUNTER *c) {
	free(c->x);
	free(c);
}

int inc_counter(COUNTER *c) {
	int i;

	i = c->digits - 1;

	c->x[i]++;
	while (c->x[i] == c->base && i > 0) {
		c->x[i] = 0;
		i--;
		c->x[i]++;
	}

	if (c->x[0] == c->base) return FALSE;

	return TRUE;
}

void print_counter(COUNTER *c) {
	int i;

	for (i=0; i<c->digits; i++) printf("%d ", c->x[i]);
	printf("\n");
}

PATH *create_path(int length) {
	PATH *path;

	path = (PATH *)my_malloc(sizeof(PATH));

	path->length = length;
	path->s_arr = (STICK **)my_calloc(length, sizeof(STICK *));

	return path;
}

void free_path(PATH *path) {
	free(path->s_arr);
	free(path);
}

PATH *create_chosen_path(STICK *s, COUNTER *c) {
	int i;
	PATH *path;
	BALL *ball;
	DLL_NODE *n;

	path = create_path(c->digits + 2);

	assert(s != NULL);
	path->s_arr[0] = s;

	assert(s->b->i == -1);
	ball = s->a;

	i = 0;
	while (i < c->digits) {
		n = get_ith_ht_elem(ball->stick_ht, c->x[i]);
		if (n == NULL) {
			free_path(path);
			return NULL;
		}
		i++;

		s = (STICK *)n->key;
		path->s_arr[i] = s;

		ball = ball == s->a ? s->b : s->a;
		if (ball->i < 0) {
			free_path(path);
			return NULL;
		}
	}

	s = get_bdy_stick(ball);
	if (s == NULL || s->b->i != -2) {
		free_path(path);
		return NULL;
	}
	i++;
	path->s_arr[i] = s;

	return path;
}

DLL_NODE *get_ith_ht_elem(HT *ht, int i) {
	int j, k;
	DLL_NODE *n;

	k = 0;
	for (j=0; j<ht->length; j++) {
		n = ht->table[j];
		while (n != NULL) {
			if (k == i) return n;
			k++;
			n = n->next;
		}
	}

	return NULL;
}

STICK *get_bdy_stick(BALL *b) {
	int i;
	DLL_NODE *n;
	STICK *s;
	HT *ht;

	ht = b->stick_ht;
	for (i=0; i<ht->length; i++) {
		n = ht->table[i];
		while (n != NULL) {
			s = (STICK *)n->key;
			if (s->b->i < 0) return s;
			n = n->next;
		}
	}

	return NULL;
}

SELECTOR *create_selector(int n, int k) {
	int i;
	SELECTOR *sel;

	sel = (SELECTOR *)my_malloc(sizeof(SELECTOR));

	sel->n = n;
	sel->k = k;
	sel->indices = (int *)my_malloc(k*sizeof(int));

	for (i=0; i<k; i++) sel->indices[i] = i;

	return sel;
}

void free_selector(SELECTOR *sel) {
	free(sel->indices);
	free(sel);
}

void reset_selector(SELECTOR *sel) {
	int i;

	for (i=0; i<sel->k; i++) sel->indices[i] = i;
}

int increment_selector(SELECTOR *sel) {
	int i;

	i = sel->k-1;

	while (i >= 0) {
		sel->indices[i]++;
		if (sel->indices[i] < sel->n - (sel->k - 1 - i)) break;
		i--;
	}

	if (i < 0) return FALSE;

	i++;
	while (i < sel->k) {
		sel->indices[i] = sel->indices[i-1] + 1;
		i++;
	}

	return TRUE;
}

void print_selector(SELECTOR *sel) {
	int i;

	for (i=0; i<sel->k; i++) {
		printf("%d ", sel->indices[i]);
	}
	printf("\n");
}

int compare_tuples(TUPLE *tp1, TUPLE *tp2) {
	if (tp1->big_t < tp2->big_t) return -1;
	if (tp1->big_t > tp2->big_t) return 1;
	if (tp1->i < tp2->i) return -1;
	if (tp1->i > tp2->i) return 1;

	return 0;
}

PATTERN *create_pattern(PATH *path, SELECTOR *sel, int n, STICK **tc_array) {
	long int min_t, big_t;
	int i, j, k;
	STICK *s, *st;
	PATTERN *pat;

	min_t = LONG_MAX;
	for (i=0; i<sel->k; i++) {
		s = path->s_arr[sel->indices[i]];
		if (s->a->big_t < min_t) min_t = s->a->big_t;
		if (s->b->big_t < min_t) min_t = s->b->big_t;
	}

	pat = (PATTERN *)my_malloc(sizeof(PATTERN));

	pat->n = sel->k;
	pat->t_arr = (TUPLE *)my_malloc(pat->n*sizeof(TUPLE));

	for (j=0; j<sel->k; j++) {
		s = path->s_arr[sel->indices[j]];
		for (k=0; k<n; k++) {
			st = tc_array[k];
			if (s->a->i != st->a->i) continue;
			if (s->a->j != st->a->j) continue;
			if (s->b->i != st->b->i) continue;
			if (s->b->j != st->b->j) continue;
			if ((s->b->big_t - s->a->big_t) != (st->b->big_t - st->a->big_t)) continue;
			break;
		}
		big_t = s->a->big_t < s->b->big_t ? s->a->big_t - min_t : s->b->big_t - min_t;
		pat->t_arr[j].big_t = big_t;
		pat->t_arr[j].i = k;
	}

	qsort(pat->t_arr, sel->k, sizeof(TUPLE), (int (*)(const void *, const void *))compare_tuples);

	return pat;
}

void free_pattern(PATTERN *pat) {
	free(pat->t_arr);
	free(pat);
}

void print_pattern(PATTERN *pat) {
	int i;

	for (i=0; i<pat->n; i++) printf("(%ld, %d) ", pat->t_arr[i].big_t, pat->t_arr[i].i);
	printf("\n");
}

