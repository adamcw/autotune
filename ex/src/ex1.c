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
#include "ex1.h"

/**
 * \brief Example 1: Implements an abstract surface code using an 8-step cycle.   
 * 
 * \param[in] argc The number of command line arguments
 * \param[in] argv An array of the command line arguments 
 *
 * \return 0
 */
int main(int argc, char **argv) {
	RECIPE_ADV *recipe;

	// Initialise, then load the args
	ARGS *args = init_args();
	load_args(args, argc, argv);

	// Create the recipe for an infinite simulation
	recipe = qc_create_recipe_adv(RECIPE_INFINITE, 2*args->d-1, 2*args->d-1, FALSE);
	
	// Generate the recipe 
	generate_recipe(args, recipe);

	// Boot up process, calculate optimal t_check
	if (args->boot == 1) {
		calculate_t_check(args, recipe);
		qc_reset_recipe_adv(recipe);
	}

	// Run the simulation
	simulate_recipe(args, recipe);
	
	// Free the recipe
	qc_free_recipe_adv(recipe);
	
	// If we opened a file for the output, it needs to be closed
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
					"-t_out T_OUT                  How often to output the current status. Useful at low error rates. [default: %d]\n", 

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
			args->big_t_max = atoi(argv[++j]);
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
 * \brief Calculates the optimal value of t_check by running a preliminary boot
 * up simulation
 * 
 * \param[in] args The arguments to the simulation
 * \param[in,out] recipe The \ref recipe to be used to calculate t_check 
 */
void calculate_t_check(ARGS *args, RECIPE_ADV *recipe) {
	long int big_t;
	int changes;
	SC_DP_QC *sc_dp_qc;
	QC *qc; 

	big_t = 0;

	// Initialise the Surface Code Depolarizing Quantum Computer
	sc_dp_qc = create_sc_dp_qc(args, recipe);
	qc = sc_dp_qc->dp_qc->qc;

	// Let QC know about the boundaries
	qc->track = FALSE;
	qc->get_boundary = get_boundary; 
	qc->boundaries = (void *)sc_dp_qc->boundaries; 
	
	// Output the initial variable state
	fprintf(args->out_raw, 
		"p = %g\n"
		"d = %d\n"
		"big_t_max = %ld\n"
		"t_check = %d\n"
		"t_delete = %d\n"
		"max_num_X = %d\n"
		"max_num_Z = %d\n"
		"boot = %d\n"
		"boot_num_X = %d\n"
		"boot_num_Z = %d\n"
		"t_check_scale = %d\n"
		"verbose = %d\n"
		"s0 = %d;\n"
		"s1 = %d;\n"
		"-s0 %d -s1 %d\n",
	args->p, args->d, args->big_t_max, args->t_check, args->t_delete, 
	args->max_num_X, args->max_num_Z, args->boot, args->boot_num_X, args->boot_num_Z, 
	args->t_check_scale, args->verbose, qc->s0, qc->s1, qc->s0, qc->s1);

	while (TRUE) {
		measure_stabilizers(sc_dp_qc, big_t);
		qc_convert_nests(qc, FALSE);
		qc_mwpm(qc, FALSE);
		correct_mts(sc_dp_qc);
		qc_trim_nests(qc, qc->unfinalized_big_t - 20);
		m_time_delete(qc->m_pr);	
		m_time_delete(qc->m_du);

		if (big_t > 0 && (big_t % args->t_check) == 0) {
			test_correct(sc_dp_qc, args->out_raw);	

			// If we have gotten to 50 time checks and have yet to find a change, 
			// double t_check or we'll be here all year.
			if (sc_dp_qc->num_X_changes == 0 && sc_dp_qc->num_Z_changes == 0 && sc_dp_qc->num_checks >= 50) {
				args->t_check <<= 1;
				fprintf(args->out_raw, "Doubling tcheck: %d\n", args->t_check);
				sc_dp_qc->num_checks = 0;
			}

			if (sc_dp_qc->num_X_changes >= args->boot_num_X || sc_dp_qc->num_Z_changes >= args->boot_num_Z) {
				changes = (sc_dp_qc->num_X_changes >= sc_dp_qc->num_Z_changes) ? sc_dp_qc->num_X_changes : sc_dp_qc->num_Z_changes; 
				fprintf(args->out_raw, "Calculating tcheck: %d %d %d %d\n", args->t_check, sc_dp_qc->num_checks, changes, args->t_check_scale);

				args->t_check = args->t_check * sc_dp_qc->num_checks / changes / args->t_check_scale;

				if (args->t_check < 1) { 
					args->t_check = 1;
				}
				
				break;
			}
		}

		// If we are verbose, then print out the current SC_DP_QC statistics
		if (args->verbose) {
			print_stats(sc_dp_qc, args->out_raw);
		}

		// Ensure that the qc has been properly advanced in big_t
		assert(qc->big_t == big_t+1);
		
		big_t++;
	}

	// Free the sc_dp_qc now that the t_check has been calculated
	free_sc_dp_qc(sc_dp_qc);
}

/**
 * \brief Generates a \ref recipe based on provided \ref args
 * 
 * \param[in] args The \ref args to generate the \ref recipe
 * \param[out] recipe The \ref recipe to be generated 
 */
void generate_recipe(ARGS *args, RECIPE_ADV *recipe) {
	SC_DP_QC *sc_dp_qc;
	long int big_t;

	sc_dp_qc = create_sc_dp_qc(args, NULL);

	big_t = 0;
	do {
		measure_stabilizers(sc_dp_qc, big_t++);
	} while (qc_boot_up_adv(sc_dp_qc->dp_qc->qc, recipe, args->switch_time, 12) != DONE);

	free_sc_dp_qc(sc_dp_qc);
}

/**
 * \brief Simulates the given \ref recipe based on provided \ref args
 * 
 * \param[in] args The \ref args to the simulation
 * \param[in] recipe  A generated \ref recipe to simulate
 */
void simulate_recipe(ARGS *args, RECIPE_ADV *recipe) {
	SC_DP_QC *sc_dp_qc;
	DP_QC *dp_qc;
	QC *qc;

	long int big_t;

	sc_dp_qc = create_sc_dp_qc(args, recipe);
	dp_qc = sc_dp_qc->dp_qc;
	qc = dp_qc->qc;

	// Disable error tracking
	qc->track = FALSE;
	
	// Let QC know about the boundaries
	qc->get_boundary = get_boundary;
	qc->boundaries = (void *)sc_dp_qc->boundaries;

	// Run the simulation
	// Now that the bootup process is complete, simulate as usual
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
		"-s0 %d -s1 %d\n",
	args->p, args->d, args->big_t_max, args->t_check, args->t_delete, 
	args->max_num_X, args->max_num_Z,
	args->verbose, qc->s0, qc->s1, qc->s0, qc->s1);

	big_t = 0;
	while (big_t <= args->big_t_max && (sc_dp_qc->num_X_changes < args->max_num_X || sc_dp_qc->num_Z_changes < args->max_num_Z)) {
		measure_stabilizers(sc_dp_qc, big_t);
		qc_convert_nests(qc, FALSE);
		qc_mwpm(qc, FALSE);
		correct_mts(sc_dp_qc);

		// not quite sure why following needs -2 rather than -1
		qc_trim_nests(qc, qc->unfinalized_big_t - 20);
		m_time_delete(qc->m_pr);
		m_time_delete(qc->m_du);

		if (args->cap_time > 0 && big_t == args->cap_time) {
			break;
		}

		if (big_t > 0 && big_t%args->t_check == 0) {
			test_correct(sc_dp_qc, args->out_raw);
		}

		if (args->verbose || (args->t_out > 0 && big_t%args->t_out == 0)) print_stats(sc_dp_qc, args->out_raw);

		assert(qc->big_t == big_t+1);
		big_t++;
	}

	/*
	printf("not in test_correct:\n");
	printf("primal:\n");
	// m_print_lattice(qc->m_pr);
	m_print_graph(qc->m_pr);
	printf("dual:\n");
	// m_print_lattice(qc->m_du);
	m_print_graph(qc->m_du);
	*/

	free_sc_dp_qc(sc_dp_qc);
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
SC_DP_QC *create_sc_dp_qc(ARGS *args, RECIPE_ADV *recipe) {
	SC_DP_QC *sc_dp_qc;
	QC *qc;
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
	sc_dp_qc->dp_qc = dp_create_dp_qc_adv(args->s0, args->s1, DE_HT_FACTOR * n * n,
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
	qc = sc_dp_qc->dp_qc->qc;

	sc_dp_qc->boundaries = (BALL **)my_malloc(8 * sizeof(BALL *));
	sc_dp_qc->bdy_s1_pr = create_boundary_set(qc, sc_dp_qc->boundaries, 0, PRIMAL_BOUNDARY, -1, 0, 0);
	sc_dp_qc->bdy_s2_pr = create_boundary_set(qc, sc_dp_qc->boundaries, 1, PRIMAL_BOUNDARY, -2, 0, 0);
	sc_dp_qc->bdy_t1_pr = create_boundary_set(qc, sc_dp_qc->boundaries, 2, PRIMAL_BOUNDARY, -1, 1, 0);
	sc_dp_qc->bdy_t2_pr = create_boundary_set(qc, sc_dp_qc->boundaries, 3, PRIMAL_BOUNDARY, -2, 1, 0);
	sc_dp_qc->bdy_s1_du = create_boundary_set(qc, sc_dp_qc->boundaries, 4, DUAL_BOUNDARY, -1, 0, 0);
	sc_dp_qc->bdy_s2_du = create_boundary_set(qc, sc_dp_qc->boundaries, 5, DUAL_BOUNDARY, -2, 0, 0);
	sc_dp_qc->bdy_t1_du = create_boundary_set(qc, sc_dp_qc->boundaries, 6, DUAL_BOUNDARY, -1, 1, 0);
	sc_dp_qc->bdy_t2_du = create_boundary_set(qc, sc_dp_qc->boundaries, 7, DUAL_BOUNDARY, -2, 1, 0);

	qc->get_boundary = get_boundary;
	qc->boundaries = (void *)sc_dp_qc->boundaries;

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
SET *create_boundary_set(QC *qc, BALL **boundaries, int id, int type, int i, int j, int t) {
	SET *set;

	set = qc_create_set_adv(qc, type, i, j, t, NULL);
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
				set_arr[is][js][1] = qc_create_set_adv(sc_dp_qc->dp_qc->qc, type, seti, setj, 1, bdy_s1);
				set_arr[is][js][0] = qc_create_set_adv(sc_dp_qc->dp_qc->qc, type, seti, setj, 2, bdy_s1);
			}
			else if (bs == bmax-1) {
				// Create a set with the second boundary
				set_arr[is][js][1] = qc_create_set_adv(sc_dp_qc->dp_qc->qc, type, seti, setj, 1, bdy_s2);
				set_arr[is][js][0] = qc_create_set_adv(sc_dp_qc->dp_qc->qc, type, seti, setj, 2, bdy_s2);
			}
			else {
				// Create a set with no boundary
				set_arr[is][js][1] = qc_create_set_adv(sc_dp_qc->dp_qc->qc, type, seti, setj, 1, NULL);
				set_arr[is][js][0] = qc_create_set_adv(sc_dp_qc->dp_qc->qc, type, seti, setj, 2, NULL);
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
	qc_free_bdy(sc_dp_qc->bdy_s1_pr);
	qc_free_bdy(sc_dp_qc->bdy_s2_pr);
	qc_free_bdy(sc_dp_qc->bdy_t1_pr);
	qc_free_bdy(sc_dp_qc->bdy_t2_pr);
	qc_free_bdy(sc_dp_qc->bdy_s1_du);
	qc_free_bdy(sc_dp_qc->bdy_s2_du);
	qc_free_bdy(sc_dp_qc->bdy_t1_du);
	qc_free_bdy(sc_dp_qc->bdy_t2_du);

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

/**
 * \brief Measure all the stabilizers
 *
 * This process is best outlined in Figure 1 of "Surface Codes: Towards
 * Practical Large-Scale Quantum Computation", arXiv:1208.0928v2
 *
 * \param[in] sc_dp_qc The \ref sc_dp_qc to measure the stabilizers from
 * \param[in] big_t The big_t of the measurements 
 */
void measure_stabilizers(SC_DP_QC *sc_dp_qc, long int big_t) {
	int n;
	DP_QC *dp_qc;
	QUBIT ***q_arr;

	n = sc_dp_qc->n;
	dp_qc = sc_dp_qc->dp_qc;
	q_arr = sc_dp_qc->q_arr;

	// The following can be summarised as so:
	//
	// Initialise X-syndromes in the Z basis, perform a Hadamard, perform CNOT
	// from the X-syndromes to the surrounding X stablizers in the order North,
	// West, East, South; then perform another Hadamard, then Measure in the Z
	// basis.
	//
	// Perform Identity gate on Z-syndromes, initialize in the Z basis, perform
	// CNOT from surrounding Z-stabilizers to each Z-syndrome in the order
	// Nother, West, East, South; then Measure in the Z basis, then perform the
	// identity gate.

	// Initialize X-syndromes in Z-basis, Dead on Z-syndromes
	init_X_syn_dead_Z_syn(dp_qc, n, q_arr);

	// Change of basis (Hadamard) on X-syndromes, Initialize Z-syndromes in Z-basis
	H_X_syn_init_Z_syn(dp_qc, n, q_arr);

	// Apply CNOT gates in order North, West, East, South 
	syn_north_cnots(dp_qc, n, q_arr);
	syn_west_cnots(dp_qc, n, q_arr);
	syn_east_cnots(dp_qc, n, q_arr);
	syn_south_cnots(dp_qc, n, q_arr);

	// Change of basis (Hadamard) on X-syndromes, Measure Z-syndromes in Z-basis
	H_X_syn_meas_Z_syn(sc_dp_qc, big_t, dp_qc, n, q_arr);
	
	// Measure X-syndromes in Z-basis, Dead on Z-syndromes
	meas_X_syn_dead_Z_syn(sc_dp_qc, big_t, dp_qc, n, q_arr);

	assert(dp_qc->qc->big_t == big_t+1);
}

/**
 * \brief Initialise X-syndromes, Dead on Z-syndromes
 * 
 * Initialize X-syndromes in Z-basis, Dead on Z-syndromes
 *
 * \param[in] dp_qc The \ref dp_qc to perform the gates with
 * \param[in] n The size of the \ref qubit array 
 * \param[in] q_arr The array of \ref qubit%s to perform the gates on 
 */
void init_X_syn_dead_Z_syn(DP_QC *dp_qc, int n, QUBIT ***q_arr) {
   	int i, j;
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			// X-syndrome (intersection) qubits
			if (i % 2 == 0 && j % 2 == 1) {
				dp_init_Z(dp_qc, q_arr[i][j]);
			}

			// Z-syndrome (face) qubits
			else if (i % 2 == 1 && j % 2 == 0) {
				dp_dead_H(dp_qc, q_arr[i][j]); //
			}

			// Data qubits (iden for as long as the X-syndromes)
			else {
				dp_iden_init_Z(dp_qc, q_arr[i][j]);
			}
		}
	} 
}

/**
 * \brief Change of basis on X-syndromes, Initialize Z-syndromes
 * 
 * Change of basis (Hadamard) on X-syndromes, Initialize Z-syndromes in Z-basis
 *
 * \param[in] dp_qc The \ref dp_qc to perform the gates with
 * \param[in] n The size of the \ref qubit array 
 * \param[in] q_arr The array of \ref qubit%s to perform the gates on 
 */
void H_X_syn_init_Z_syn(DP_QC *dp_qc, int n, QUBIT ***q_arr) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			// X-syndromes, perform a Hadamard gate
			if (i % 2 == 0 && j % 2 == 1) {
				dp_H(dp_qc, q_arr[i][j]);
			}

			// Z-syndrome, initialize (note that init Z may take longer, or
			// shorter than hadamard, this is why Z-syndromes in the previous
			// loop were dead for the length of H rather than init Z) 
			else if (i % 2 == 1 && j % 2 == 0) {
				dp_init_Z(dp_qc, q_arr[i][j]);
			}

			// Data qubits (iden for as long as the X-syndromes)
			else {
				dp_iden_H(dp_qc, q_arr[i][j]);
			}
		}
	}
}

/**
 * \brief Applies the CNOT gates North of each syndrome qubit.
 * 
 * \param[in] dp_qc The \ref dp_qc to perform the gates with
 * \param[in] n The size of the \ref qubit array 
 * \param[in] q_arr The array of \ref qubit%s to perform the gates on 
 */
void syn_north_cnots(DP_QC *dp_qc, int n, QUBIT ***q_arr) {
	int i, j;

	// X-syndromes along the top have no data qubit to their north, so perform an
	// identity instead
	for (j = 1; j < n; j += 2) {
		dp_iden_cnot(dp_qc, q_arr[0][j]);
	}
	for (i = 0; i < n-1; i++) {
		// Horizontal edge data qubits to Z-syndromes
		for (j = 0; j < n; j += 2) {
			dp_cnot(dp_qc, q_arr[i][j], q_arr[i+1][j]);
		}

		i++;

		// X-syndromes to vertical edge data qubits
		for (j = 1; j < n; j += 2) {
			dp_cnot(dp_qc, q_arr[i+1][j], q_arr[i][j]);
		}
	}
	// Data qubits along the bottom have no syndrome to their south, so perform
	// an identity instead
	for (j = 0; j < n; j += 2) {
		dp_iden_cnot(dp_qc, q_arr[i][j]);
	}
}

/**
 * \brief Applies the CNOT gates West of each syndrome qubit.
 * 
 * \param[in] dp_qc The \ref dp_qc to perform the gates with
 * \param[in] n The size of the \ref qubit array 
 * \param[in] q_arr The array of \ref qubit%s to perform the gates on 
 */
void syn_west_cnots(DP_QC *dp_qc, int n, QUBIT ***q_arr) {
	int i, j;

	for (i = 0; i < n-1; i++) {
		// From X-syndromes to horizontal edge data qubits
		for (j = 0; j < n-1; j += 2) {
			dp_cnot(dp_qc, q_arr[i][j+1], q_arr[i][j]);
		}
		// The data qubit at right end of the grid do not have a syndrome to
		// their east, so perform an identity instead
		dp_iden_cnot(dp_qc, q_arr[i][j]);
		
		i++;

		// The Z stabilizers at the left of the grid do not have a data qubit
		// to their west, so perform an identity instead.
		dp_iden_cnot(dp_qc, q_arr[i][0]);
		// From vertical edge data qubits to Z-syndromes
		for (j = 1; j < n; j += 2) {
			dp_cnot(dp_qc, q_arr[i][j], q_arr[i][j+1]);
		}
	}
	// The perform the first part of the above loop, as there is an odd number of
	// rows from top to bottom and the above loop stops one row short.
	for (j = 0; j < n-1; j += 2) {
		dp_cnot(dp_qc, q_arr[i][j+1], q_arr[i][j]);
	}
	dp_iden_cnot(dp_qc, q_arr[i][j]);
}

/**
 * \brief Applies the CNOT gates East of each syndrome qubit.
 * 
 * \param[in] dp_qc The \ref dp_qc to perform the gates with
 * \param[in] n The size of the \ref qubit array 
 * \param[in] q_arr The array of \ref qubit%s to perform the gates on 
 */
void syn_east_cnots(DP_QC *dp_qc, int n, QUBIT ***q_arr) {
	int i, j;
	
	for (i = 0; i < n-1; i++) {
		// The data qubits to the left of the grid do not have a syndrome to
		// their west, so perform an identity instead
		dp_iden_cnot(dp_qc, q_arr[i][0]);
		// From X-syndromes to horizontal edge data qubits
		for (j = 1; j < n; j += 2) {
			dp_cnot(dp_qc, q_arr[i][j], q_arr[i][j+1]);
		}

		i++;
		
		// From vertical edge qubits to Z-syndromes
		for (j = 0; j < n-1; j += 2) {
			dp_cnot(dp_qc, q_arr[i][j+1], q_arr[i][j]);
		}
		// The Z syndromes at the right of the grid do not have a data qubit
		// to their east, so perform an identity instead.
		dp_iden_cnot(dp_qc, q_arr[i][j]);
	}

	// The perform the first part of the above loop, as there is an odd number of
	// rows from top to bottom and the above loop stops one row short.
	dp_iden_cnot(dp_qc, q_arr[i][0]);
	for (j = 1; j < n; j += 2) {
		dp_cnot(dp_qc, q_arr[i][j], q_arr[i][j+1]);
	}
}

/**
 * \brief Applies the CNOT gates South of each syndrome qubit.
 * 
 * \param[in] dp_qc The \ref dp_qc to perform the gates with
 * \param[in] n The size of the \ref qubit array 
 * \param[in] q_arr The array of \ref qubit%s to perform the gates on 
 */
void syn_south_cnots(DP_QC *dp_qc, int n, QUBIT ***q_arr) {
	int i, j;
	
	// Data qubits at the top of the grid have no syndromes to their north, so
	// perform an identity instead.
	for (j = 0; j < n; j += 2) {
		dp_iden_cnot(dp_qc, q_arr[0][j]);
	}

	for (i = 0; i < n-1; i++) {
		// From X-syndromes to vertical edge data qubits
		for (j = 1; j < n; j += 2) {
			dp_cnot(dp_qc, q_arr[i][j], q_arr[i+1][j]);
		}

		i++;

		// From horizontal data qubits to Z-syndromes
		for (j = 0; j < n; j += 2) {
			dp_cnot(dp_qc, q_arr[i+1][j], q_arr[i][j]);
		}
	}

	// X-syndromes along the bottom of the grid have no data qubits to their
	// south, so perform an identity instead
	for (j = 1; j < n; j += 2) {
		dp_iden_cnot(dp_qc, q_arr[i][j]);
	}
}

/**
 * \brief Change of basis on X-syndromes, Measure Z-syndromes
 * 
 * \param[in] sc_dp_qc The \ref dp_qc to perform the gates with
 * \param[in] big_t The big_t 
 * \param[in] dp_qc The \ref dp_qc to perform the gates with
 * \param[in] n The size of the \ref qubit array 
 * \param[in] q_arr The array of \ref qubit%s to perform the gates on 
 */
void H_X_syn_meas_Z_syn(SC_DP_QC *sc_dp_qc, long int big_t, DP_QC *dp_qc, int n, QUBIT ***q_arr) {
	int i, j, is, js, lay1, lay2;
	SET *bdy_s1_du, *bdy_s2_du; 
	SYNDROME ***syn_arr_du;
	SET ****set_arr_du;
	
	bdy_s1_du = sc_dp_qc->bdy_s1_du;
	bdy_s2_du = sc_dp_qc->bdy_s2_du;
	syn_arr_du = sc_dp_qc->syn_arr_du;
	set_arr_du = sc_dp_qc->set_arr_du;

	lay1 = (big_t + 1) % 2;
	lay2 = big_t % 2;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {

			// X-syndromes, perform Hadamard
			if (i % 2 == 0 && j % 2 == 1) {
				dp_H(dp_qc, q_arr[i][j]);
			}

			// Z-syndromes, measure.
			else if (i % 2 == 1 && j % 2 == 0) {
				is = i >> 1;
				js = j >> 1;

				// Measure the Z-syndrome
				dp_meas_Z(dp_qc, q_arr[i][j], set_arr_du[is][js][lay1], set_arr_du[is][js][lay2]);

				
				qc_unassociate_syndrome(syn_arr_du[is][js]);
				qc_associate_syndrome(set_arr_du[is][js][lay2], syn_arr_du[is][js]);
				if (i == 1) {
					set_arr_du[is][js][lay1] = qc_create_set_adv(sc_dp_qc->dp_qc->qc, DUAL, i, j, 2, bdy_s1_du);
				}
				else if (i == n-2) {
					set_arr_du[is][js][lay1] = qc_create_set_adv(sc_dp_qc->dp_qc->qc, DUAL, i, j, 2, bdy_s2_du);
				}
				else {
					set_arr_du[is][js][lay1] = qc_create_set_adv(sc_dp_qc->dp_qc->qc, DUAL, i, j, 2, NULL);
				}						
				qc_insert_set(dp_qc->qc, set_arr_du[is][js][lay1]);
			}

			// Data qubits, identity
			else {
				dp_iden_H(dp_qc, q_arr[i][j]);
			}
		}
	}
}

/**
 * \brief Measure X-syndromes, Dead Z-syndromes
 * 
 * \param[in] sc_dp_qc The \ref dp_qc to perform the gates with
 * \param[in] big_t The big_t 
 * \param[in] dp_qc The \ref dp_qc to perform the gates with
 * \param[in] n The size of the \ref qubit array 
 * \param[in] q_arr The array of \ref qubit%s to perform the gates on 
 */
void meas_X_syn_dead_Z_syn(SC_DP_QC *sc_dp_qc, long int big_t, DP_QC *dp_qc, int n, QUBIT ***q_arr) {
	int i, j, is, js, lay1, lay2;
	SET *bdy_s1_pr, *bdy_s2_pr;
	SYNDROME ***syn_arr_pr;
	SET ****set_arr_pr;

	bdy_s1_pr = sc_dp_qc->bdy_s1_pr;
	bdy_s2_pr = sc_dp_qc->bdy_s2_pr;
	syn_arr_pr = sc_dp_qc->syn_arr_pr;
	set_arr_pr = sc_dp_qc->set_arr_pr;

	lay1 = (big_t + 1) % 2;
	lay2 = big_t % 2;

	// Identity and Measure Z
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i % 2 == 0 && j % 2 == 1) {
				is = i >> 1;
				js = j >> 1;
				
				// X-syndrome
				dp_meas_Z(dp_qc, q_arr[i][j], set_arr_pr[is][js][lay1], set_arr_pr[is][js][lay2]);
				qc_unassociate_syndrome(syn_arr_pr[is][js]);
				qc_associate_syndrome(set_arr_pr[is][js][lay2], syn_arr_pr[is][js]);
				if (j == 1) {
					set_arr_pr[is][js][lay1] = qc_create_set_adv(dp_qc->qc, PRIMAL, i, j, 2, bdy_s1_pr);
				}
				else if (j == n-2) {
					set_arr_pr[is][js][lay1] = qc_create_set_adv(dp_qc->qc, PRIMAL, i, j, 2, bdy_s2_pr);
				}
				else {
					set_arr_pr[is][js][lay1] = qc_create_set_adv(dp_qc->qc, PRIMAL, i, j, 2, NULL);
				}						
				qc_insert_set(dp_qc->qc, set_arr_pr[is][js][lay1]);
			}
			else if (i % 2 == 1 && j % 2 == 0) {
				dp_dead_H(dp_qc, q_arr[i][j]); //
			}
			else {
				dp_iden_meas_Z(dp_qc, q_arr[i][j]);
			}
		}
	}
}

/**
 * \brief Processes an augmented edge by applying corrections along it 
 * 
 * \param[in] ae The \ref ae to process
 * \param[in] d The distance code
 * \param[in] frame The Pauli frame to process the edge across 
 */
void process_aug_edge(AUG_EDGE *ae, int d, int **frame) {
	int i1, j1, i2, j2;

	// Get the coordinates at the start and end of the edge
	i1 = ae->va->i;
	j1 = ae->va->j;
	i2 = ae->vb->i;
	j2 = ae->vb->j;

	//printf("%d - %d - %d - %d\n", i1, j1, i2, j2);
	assert(i1 >= 0 && j1 >= 0);

	// Handle connections to Top/Left/Inital Time boundaries
	if (i2 == -1) {
		// Spatial
		if (j2 == 0) {
			// X-syndrome, therefore Left boundary
			if (i1%2 == 0) {
				i2 = i1;
				j2 = -1;
			}
			
			// Z-syndrome, therefore Top boundary
			else {
				i2 = -1;
				j2 = j1;
			}
		} 
		
		// Temporal
		else {
			// Could flip measurements to temporal boundary
			return;
		}
	}

	// Handle connections to the Bottom/Right/Final Time boundaries
	else if (i2 == -2) {
		// Spatial
		if (j2 == 0) {
			// X-syndrome, therefore Right boundary
			if (i1%2 == 0) {
				i2 = i1;
				j2 = 2*d - 1;
			}
			
			// Z-syndrome, therefore Bottom boundary
			else {
				i2 = 2*d - 1;
				j2 = j1;
			}
		}

		// Temporal
		else {
			// Could flip measurements to temporal boundary
			return;
		}
	}

	// X-syndromes exist on even qubits, apply Z corrections
	if (i1 % 2 == 0) {
		// Step from (i1, j1) to (i2, j2) applying Z corrections to each qubit
		// along the path. 
		while (i1 < i2) {
			i1++;
			frame[i1][j1] ^= Z;
			i1++;
		}
		while (i1 > i2) {
			i1--;
			frame[i1][j1] ^= Z;
			i1--;
		}
		while (j1 < j2) {
			j1++;
			frame[i1][j1] ^= Z;
			j1++;
		}
		while (j1 > j2) {
			j1--;
			frame[i1][j1] ^= Z;
			j1--;
		}
	}
	
	// Z-syndromes exist on odd qubits, apply X corrections
	else {
		// Step from (i1, j1) to (i2, j2) applying X corrections to each qubit
		// along the path. 
		while (i1 < i2) {
			i1++;
			frame[i1][j1] ^= X;
			i1++;
		}
		while (i1 > i2) {
			i1--;
			frame[i1][j1] ^= X;
			i1--;
		}
		while (j1 < j2) {
			j1++;
			frame[i1][j1] ^= X;
			j1++;
		}
		while (j1 > j2) {
			j1--;
			frame[i1][j1] ^= X;
			j1--;
		}
	}
}

/**
 * \brief Process augmented edges from a \ref matching
 * 
 * \param[in] m The \ref matching containing the augmented edges
 * \param[in] d The distance of the matching
 * \param[in] frame The Pauli frame to process the edges on
 */
void process_aug_edges(MATCHING *m, int d, int **frame) {
	AUG_EDGE *ae;

	ae = m_get_aug_edge(m);
	while (ae != NULL) {
		process_aug_edge(ae, d, frame);
		m_delete_aug_edge(ae);
		ae = m_get_aug_edge(m);
	}
}

/**
 * \brief Corrects the \ref measurements in a ref sc_dp_qc
 * 
 * Processes the augmented edges from the primal and dual matchings
 *
 * \param[in] sc_dp_qc The \ref sc_dp_qc to correct 
 */
void correct_mts(SC_DP_QC *sc_dp_qc) {
	QC *qc;
	qc = sc_dp_qc->dp_qc->qc;

	process_aug_edges(qc->m_pr, sc_dp_qc->d, sc_dp_qc->frame);
	process_aug_edges(qc->m_du, sc_dp_qc->d, sc_dp_qc->frame);
}

/**
 * \brief Performs a round of perfect stabilizer measurement, then determines if there
 * has been a logical X and/or Z error.
 * 
 * \param[in] sc_dp_qc The \ref sc_dp_qc to test
 * \param[out] out A file pointer to where the output stats of the test should
 * be sent.
 */
void test_correct(SC_DP_QC *sc_dp_qc, FILE *out) {
	SC_DP_QC *sc_dp_qc2;
	DP_QC *dp_qc2;
	QC *qc2;
	int n, i, j, count;

	//printf("Test correct\n");
	sc_dp_qc2 = copy_sc_dp_qc(sc_dp_qc);
	dp_qc2 = sc_dp_qc2->dp_qc;
	qc2 = dp_qc2->qc;
	n = sc_dp_qc2->n;

	// Enable perfect gates and perform a perfect round of stabilizer
	// measurements, followed by matching, then correction.
	qc2->perfect_gates = true;

	qc2->m_pr->undo_flag = TRUE;
	qc2->m_du->undo_flag = TRUE;

	measure_stabilizers(sc_dp_qc2, qc2->big_t);

	qc_finalize_nests(qc2, qc2->big_t);
	qc_convert_nests(qc2, TRUE);

	qc_mwpm(qc2, TRUE);
	correct_mts(sc_dp_qc2);

	// Track how many test corrects have been performed for this simulation
	sc_dp_qc->num_checks++;

	// Count the number of Z errors that have occured along the vertical
	// boundary. An Z-error is a disagreement between the state of the qubit in
	// the qc and the state of the correction in the frame. 
	//
	// If there is a Z error on the qubit, but no Z correction, or the reverse
	// with a Z correction but no Z error on the qubit, then an error along
	// the boundary has occured. 
	count = 0;
	for (i = 0; i < n; i += 2) {
		if ((dp_contains_Z(sc_dp_qc2->q_arr[i][0]->e) && ((sc_dp_qc2->frame[i][0]&Z) != Z)) ||
			(!dp_contains_Z(sc_dp_qc2->q_arr[i][0]->e) && ((sc_dp_qc2->frame[i][0]&Z) == Z))) {
			count++;
		}
	}

	// If the parity of the number of errors along the vertical boundary has
	// changed since the last check, then a logical error has occured.
	if (count % 2 != sc_dp_qc->last_Z_check) {
		sc_dp_qc->last_Z_check = count%2;
		sc_dp_qc->num_Z_changes++;
		print_stats(sc_dp_qc, out);
	}

	// Count the number of X errors that have occured along the horizontal
	// boundary. An X-error is a disagreement between the state of the qubit in
	// the qc and the state of the correction in the frame. 
	//
	// If there is an X error on the qubit, but no X correction, or the reverse
	// with an X correction but no X error on the qubit, then an error along
	// the boundary has occured. 
	count = 0;
	for (j = 0; j < n; j += 2) {
		if ((dp_contains_X(sc_dp_qc2->q_arr[0][j]->e) && ((sc_dp_qc2->frame[0][j]&X) != X)) ||
			(!dp_contains_X(sc_dp_qc2->q_arr[0][j]->e) && ((sc_dp_qc2->frame[0][j]&X) == X))) {
			count++;
		}
	}

	// If the parity of the number of errors along the horizontal boundary has
	// changed since the last check, then a logical error has occured.	
	if (count % 2 != sc_dp_qc->last_X_check) {
		sc_dp_qc->last_X_check = count%2;
		sc_dp_qc->num_X_changes++;
		print_stats(sc_dp_qc, out);
	}

	qc_undo_mwpm(qc2);

	free_sc_dp_qc_copy(sc_dp_qc2);

	sc_dp_qc->dp_qc->qc->m_pr->undo_flag = FALSE;
	sc_dp_qc->dp_qc->qc->m_du->undo_flag = FALSE;

	//printf("End test correct\n");
}

/**
 * \brief Outputs the current statistics of the simulation including: Current
 * time, current big_t, the parity of X and Z in the last check, the number of
 * checks, and the number of X and Z logical changes.
 * 
 * \param[in] sc_dp_qc The \ref sc_dp_qc to print the statistics of
 * \param[out] out A file pointer to where the output stats of the test should
 * be sent.
 */
void print_stats(SC_DP_QC *sc_dp_qc, FILE *out) {
	double t2;

	t2 = double_time();
	fprintf(out, ">>> %.2f | %ld | Last: %d %d | Checks: %d %d | Changes: %d %d\n", 
		t2 - sc_dp_qc->t1, 
		sc_dp_qc->dp_qc->qc->big_t, 
		sc_dp_qc->last_X_check, 
		sc_dp_qc->last_Z_check, 
		sc_dp_qc->num_checks, 
		sc_dp_qc->num_checks, 
		sc_dp_qc->num_X_changes, 
		sc_dp_qc->num_Z_changes
	);
	fflush(out);
}

/**
 * \brief Prints all es from the \ref qubit array
 * 
 * \param[in] n The size of the qubit array
 * \param[in] q_arr The \ref qubit array to print the es of 
 */
void print_all_es(int n, QUBIT ***q_arr) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%d ", q_arr[i][j]->e);
		}
		printf("\n");
	}
	printf("\n");
}

/**
 * \brief Prints all errors from a \ref qubit array
 * 
 * \param[in] n The size of the qubit array
 * \param[in] q_arr The \ref qubit array to print the errors from 
 */
void print_all_errors(int n, QUBIT ***q_arr) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("q(%d, %d)\n", i, j);
			qc_print_errors(q_arr[i][j]->error_cdll);
		}
	}
}

/**
 * \brief Infers the direction of the boundary from a given coordinate.
 * 
 * \param[in] dots A \ref ht of \ref dot\s 
 * \param[in] i1 The i coordinate
 * \param[in] j1 The j coordinate
 * \param[in] t1 The t coordinate
 * \param[in] t_offset The t_offset of the overall simulation 
 */
int *infer_boundary(HT *dots, const int i1, const int j1, const int t1, const int t_offset) {
	int i2, j2, t2;
   
	int found_bdy, found_dot;
	int hash;
	int k;

	int *coord;

	const int directions[78] = {
		 0,  0, -1,
		 0,  0,  1,
		 2,  0,  0,
		-2,  0,  0,
		 0,  2,  0,
		 0, -2,  0,
		 2,  2,  0,
		-2, -2,  0,
		 2, -2,  0,
		-2,  2,  0,
		 2,  0,  1,
		-2,  0,  1,
		 0,  2,  1,
		 0, -2,  1,
		 2,  0, -1,
		-2,  0, -1,
		 0,  2, -1,
		 0, -2, -1,
		 2,  2,  1,
		-2, -2,  1,
		 2, -2,  1,
		-2,  2,  1,
		 2,  2, -1,
		-2, -2, -1,
		 2, -2, -1,
		-2,  2, -1
	};

	coord = (int *)my_malloc(sizeof(int) * 3);

	// Look in all 26 directions for a dot-less position
	// This is inferred to be where the boundary is connected to.
	found_bdy = FALSE;
	for (k = 0; k < 78; ) {
		i2 = i1 + directions[k++];
		j2 = j1 + directions[k++];
		t2 = t1 + directions[k++];

		// Don't look for past temporal boundaries if we have used t_delete
		if (t_offset > 0 && t2 < 0) {
			continue;
		}

		found_dot = FALSE;
		hash = (i2 + j2 + t2) % INT_MAX;
		DLL_NODE *node = dots->table[hash];
		while (node != NULL) {
			DOT *dot = (DOT *)node->key;

			if (dot->i == i2 && dot->j == j2 && dot->t == t2) {
				found_dot = TRUE;
				break;
			}

			node = node->next;
		}

		// If we didn't find a dot, then we know where the boundary is.
		if (found_dot == FALSE) {
			found_bdy = TRUE;

			coord[0] = i2;
			coord[1] = j2;
			coord[2] = t2;
			return coord;
		}
	}

	// If no boundary location is found, we have an error.
	if (!found_bdy) {
		printf("Boundary not found for: (%d, %d, %d)\n", i1, j1, t1);
		assert(found_bdy);
	}

	coord[0] = -1;
	coord[1] = -1;
	coord[2] = -1;
	return coord;
}

/**
 * \brief Verifies the structure and links of balls and dots
 *
 * Ensures that all dots point to the correct balls and vica-versa.
 * 
 * \param[in] qc The \ref qc to be verified
 */
void verify_balls_and_dots(QC *qc) {
	CDLL_NODE *n;
	BALL *b;
	DOT *d;

	// Verify that for every ball in the primal nest, that if the ball has an
	// associated dot, that the dot's ball is either set to this ball, or the
	// dot's ball's copy is set to the ball. 
	n = qc->nest_pr->ball_cdll->next;
	while (n != qc->nest_pr->ball_cdll) {
		b = (BALL *)n->key;
		if (b->dot != NULL) assert(b->dot->ball == b || b->dot->ball->copy == b);
		n = n->next;
	}

	// Verify that for every dot in the primal matching, that the dot is
	// associated with a ball, and that the ball is associated with the dot. 
	n = qc->m_pr->dots->next;
	while (n != qc->m_pr->dots) {
		d = (DOT *)n->key;
		assert(d->ball != NULL);
		assert(d->ball->dot == d);
		n = n->next;
	}

	// Verify that for every ball in the dual nest, that if the ball has an
	// associated dot, that the dot's ball is either set to this ball, or the
	// dot's ball's copy is set to the ball. 
	n = qc->nest_du->ball_cdll->next;
	while (n != qc->nest_du->ball_cdll) {
		b = (BALL *)n->key;
		if (b->dot != NULL) assert(b->dot->ball == b || b->dot->ball->copy == b);
		n = n->next;
	}

	// Verify that for every dot in the dual matching, that the dot is
	// associated with a ball, and that the ball is associated with the dot. 
	n = qc->m_du->dots->next;
	while (n != qc->m_du->dots) {
		d = (DOT *)n->key;
		assert(d->ball != NULL);
		assert(d->ball->dot == d);
		n = n->next;
	}
}

/**
 * \brief Prints statistics associated with the \ref nest\s of a \ref qc
 * 
 * \param[in] 
 */
void print_nest_nums(QC *qc) {
	int pb, ps, db, ds;	
	//int pd, pl, dd, dl;
	
	CDLL_NODE *n;

	pb = 0;
	n = qc->nest_pr->ball_cdll->next;
	while (n != qc->nest_pr->ball_cdll) {
		n = n->next;
		pb++;
	}

	ps = 0;
	n = qc->nest_pr->stick_cdll->next;
	while (n != qc->nest_pr->stick_cdll) {
		n = n->next;
		ps++;
	}

	/*
	pd = 0;
	n = qc->m_pr->dots->next;
	while (n != qc->m_pr->dots) {
		n = n->next;
		pd++;
	}

	pl = 0;
	n = qc->m_pr->lines->next;
	while (n != qc->m_pr->lines) {
		n = n->next;
		pl++;
	}
	*/

	db = 0;
	n = qc->nest_du->ball_cdll->next;
	while (n != qc->nest_du->ball_cdll) {
		n = n->next;
		db++;
	}

	ds = 0;
	n = qc->nest_du->stick_cdll->next;
	while (n != qc->nest_du->stick_cdll) {
		n = n->next;
		ds++;
	}

	/*
	dd = 0;
	n = qc->m_du->dots->next;
	while (n != qc->m_du->dots) {
		n = n->next;
		dd++;
	}

	dl = 0;
	n = qc->m_du->lines->next;
	while (n != qc->m_du->lines) {
		n = n->next;
		dl++;
	}
	*/

	// fprintf(stderr, "pb: %d, ps: %d, pd: %d, pl: %d, db: %d, ds: %d, dd: %d, dl: %d\n", pb, ps, pd, pl, db, ds, dd, dl);
	fprintf(stderr, "pb: %d, ps: %d, db: %d, ds: %d\n", pb, ps, db, ds);
}

