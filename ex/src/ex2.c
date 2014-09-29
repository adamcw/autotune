#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>

#include "depolar_loss/depolar_loss.h"
#include "memory/memory.h"
#include "my_time/my_time.h"
#include "ex2.h"

//#define PRINT_ROUNDS 1

/**
 * \brief Example 2: Implements Topological Cluster State  
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
	args->t_delete = -1;
	args->max_num_X = 10000;
	args->max_num_Z = 10000;
	args->verbose = 0;
	args->screen = 0;
	args->t_out = 0;
	args->t_delay = 0;
	args->track = 0;
	args->loss = 1;
	args->p_loss = 0.0;
	args->p_depolar = 0.0;
	args->t_check_scale = 10;
	args->boot = 1;
	args->boot_num_X = 50;
	args->boot_num_Z = 50;
	args->cap_time = 0;
	args->switch_time = 3L;
	args->flush_file = 0;

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
					"-switch_time SWITCH_TIME      When to switch to the repeating layer in the boot up phase [default: %ld]\n"
					"-t_check_scale T_CHECK_SCALE  A scaling factor applied to t_check to ensure it is reasonable. [default: %d]\n"
					"-boot BOOT	                   Whether or not to have a boot up phase to calculate t_check [default: %d]\n"
					"-boot_num_X BOOT_NUM_X        The number of X changes to wait for during the boot up phase. [default: %d]\n"
					"-boot_num_Z BOOT_NUM_Z        The number of Z changes to wait for during the boot up phase. [default: %d]\n"

					"-cap_time CAP_TIME            [default: %d]\n"
					"-ems EMS                      The folder containing the error models to use. [default: %s]\n"
					"-screen SCREEN                Whether or not to output to the screen. [default: %d]\n"
					"-t_out T_OUT                  How often to output the current status. Useful at low error rates. [default: %d]\n"
					"-track TRACK                  Whether or not to use full error tracking (instead of recipe) [default: %d]\n"
					"-loss LOSS					   Base probability of loss occuring [default: %f]\n"
					"-depolar DEPOLAR			   Base probability of depolarization when interacting with a lost qubit. [default: %f]\n"
					"-flush_file FLUSH_FILE		   Whether or not to close and reopen the file after every error, flushing the writes. [default: %d]", 

			args->p, args->d, args->big_t_max, args->s0, args->s1, 
			args->t_check, args->t_delete, args->max_num_X, args->max_num_Z, args->verbose, 
			args->switch_time, args->t_check_scale, args->boot, args->boot_num_X, args->boot_num_Z, 
			args->cap_time, args->ems, args->screen, args->t_out, args->track, args->p_loss, args->p_depolar, args->flush_file);

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
		else if (!strcmp(argv[j], "-t_delay")) {
			args->t_delay = atoi(argv[++j]);
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
		else if (!strcmp(argv[j], "-flush_file")) {
			args->flush_file = atoi(argv[++j]);
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
		else if (!strcmp(argv[j], "-track")) {
			args->track = 1;
		}
		else if (!strcmp(argv[j], "-loss")) {
			args->p_loss = atof(argv[++j]);
			args->loss = 1;
		}
		else if (!strcmp(argv[j], "-depolar")) {
			args->p_depolar = atof(argv[++j]);
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

	// Set t_delete to be 5 times the code distance by default
	if (args->t_delete == -1) {
		args->t_delete = args->d * 5;
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
	SC_DPL_QC *sc_dpl_qc;
	DPL_QC *dpl_qc;
	QC *qc; 

	big_t = 0;

	// Initialise the Surface Code Depolarizing Quantum Computer
	sc_dpl_qc = create_sc_dpl_qc(args, recipe, args->loss);
	dpl_qc = sc_dpl_qc->dpl_qc;
	qc = dpl_qc->qc;

	if (args->s0 == 0) {
		args->s0 = qc->s0;
	}
	if (args->s1 == 0) {
		args->s1 = qc->s1;
	}

	// Let QC know about the boundaries
	qc->track = FALSE;
	qc->m_pr->t_delay = args->t_delay;
	qc->m_du->t_delay = args->t_delay;
	
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

	dpl_qc->p = 0;
	dpl_qc->p_loss = 0;
	measure_stabilizers(sc_dpl_qc);
	qc_convert_nests(qc, FALSE);
	big_t++;
	dpl_qc->p = args->p;
	dpl_qc->p_loss = args->p_loss;

	while (TRUE) {
		measure_stabilizers(sc_dpl_qc); 
		qc_convert_nests(qc, FALSE);
		qc_mwpm(qc, FALSE);
		correct_mts(sc_dpl_qc);
		qc_trim_nests(qc, qc->unfinalized_big_t - 20);
		m_time_delete(qc->m_pr);	
		m_time_delete(qc->m_du);

		if (big_t > 0 && (big_t % args->t_check) == 0) {
			test_correct(sc_dpl_qc, &args->out_raw);	

			// If we have gotten to 50 time checks and have yet to find a change, 
			// double t_check or we'll be here all year.
			if (sc_dpl_qc->num_X_changes == 0 && sc_dpl_qc->num_Z_changes == 0 && sc_dpl_qc->num_checks >= 50) {
				args->t_check <<= 1;
				fprintf(args->out_raw, "Doubling tcheck: %d\n", args->t_check);
				sc_dpl_qc->num_checks = 0;
			}

			if (sc_dpl_qc->num_X_changes >= args->boot_num_X || sc_dpl_qc->num_Z_changes >= args->boot_num_Z) {
				changes = (sc_dpl_qc->num_X_changes >= sc_dpl_qc->num_Z_changes) ? sc_dpl_qc->num_X_changes : sc_dpl_qc->num_Z_changes; 
				fprintf(args->out_raw, "Calculating tcheck: %d %d %d %d\n", args->t_check, sc_dpl_qc->num_checks, changes, args->t_check_scale);

				args->t_check = args->t_check * sc_dpl_qc->num_checks / changes / args->t_check_scale;

				if (args->t_check < 1) { 
					args->t_check = 1;
				}
				
				recipe->in_cycle = 0;
				break;
			}
		}

		// If we are verbose, then print out the current SC_DPL_QC statistics
		if (args->verbose) {
			print_stats(sc_dpl_qc, &args->out_raw);
		}

		// Ensure that the qc has been properly advanced in big_t
		//assert(qc->big_t == big_t+1);
		
		big_t++;
	}

	// Free the sc_dpl_qc now that the t_check has been calculated
	free_sc_dpl_qc(sc_dpl_qc);
}

/**
 * \brief Generates a \ref recipe based on provided \ref args
 * 
 * \param[in] args The \ref args to generate the \ref recipe
 * \param[out] recipe The \ref recipe to be generated 
 */
void generate_recipe(ARGS *args, RECIPE_ADV *recipe) {
	SC_DPL_QC *sc_dpl_qc;

	sc_dpl_qc = create_sc_dpl_qc(args, NULL, FALSE);

	do {
		measure_stabilizers(sc_dpl_qc);
	} while (qc_boot_up_adv(sc_dpl_qc->dpl_qc->qc, recipe, args->switch_time, 12) != DONE);

	free_sc_dpl_qc(sc_dpl_qc);
}



/**
 * \brief Simulates the given \ref recipe based on provided \ref args
 * 
 * \param[in] args The \ref args to the simulation
 * \param[in] recipe  A generated \ref recipe to simulate
 */
void simulate_recipe(ARGS *args, RECIPE_ADV *recipe) {
	SC_DPL_QC *sc_dpl_qc;
	DPL_QC *dpl_qc;
	QC *qc;

	long int big_t;

	if (args->track) {
		sc_dpl_qc = create_sc_dpl_qc(args, NULL, args->loss);
	} else {
		sc_dpl_qc = create_sc_dpl_qc(args, recipe, args->loss);
	}
	dpl_qc = sc_dpl_qc->dpl_qc;
	qc = dpl_qc->qc;

	qc->track = args->track;
	qc->m_pr->t_delay = args->t_delay;
	qc->m_du->t_delay = args->t_delay;

	// Run the simulation
	// Now that the bootup process is complete, simulate as usual
	fprintf(args->out_raw, 
		"p = %g\n"
		"d = %d\n"
		"big_t_max = %ld\n"
		"new t_check = %d\n"
		"t_delete = %d\n"
		"loss = %d\n"
		"p_loss = %g\n"
		"p_depolar = %g\n"
		"max_num_X = %d\n"
		"max_num_Z = %d\n"
		"verbose = %d\n"
		"s0 = %d;\n"
		"s1 = %d;\n"
		"-s0 %d -s1 %d\n",
	args->p, args->d, args->big_t_max, args->t_check, args->t_delete, args->loss, args->p_loss, args->p_depolar,
	args->max_num_X, args->max_num_Z,
	args->verbose, qc->s0, qc->s1, qc->s0, qc->s1);

	if (args->out_raw != stdout && args->flush_file) {
		fflush(args->out_raw);
		fclose(args->out_raw);
		args->out_raw = (FILE *)fopen("out_raw", "a");
	}

	dpl_qc->p = 0;
	dpl_qc->p_loss = 0;

	big_t = 0;
	measure_stabilizers(sc_dpl_qc);
	qc_convert_nests(qc, FALSE);
	big_t++;

	dpl_qc->p = args->p;
	dpl_qc->p_loss = args->p_loss;

	while (big_t <= args->big_t_max && (sc_dpl_qc->num_X_changes < args->max_num_X || sc_dpl_qc->num_Z_changes < args->max_num_Z)) {
		measure_stabilizers(sc_dpl_qc);
		qc_convert_nests(qc, FALSE);

		//printf("%d %d\n", args->t_delay, args->t_delay);
		qc_mwpm(qc, FALSE);
		correct_mts(sc_dpl_qc);
		//printf("---\n");

		// not quite sure why following needs -2 rather than -1
		qc_trim_nests(qc, qc->unfinalized_big_t - 20);
		m_time_delete(qc->m_pr);
		m_time_delete(qc->m_du);

		if (args->cap_time > 0 && big_t == args->cap_time) {
			break;
		}

		if (big_t > 0 && big_t % args->t_check == 0) {
			test_correct(sc_dpl_qc, &args->out_raw);
		}

		if (args->verbose || (args->t_out > 0 && big_t % args->t_out == 0)) {
			print_stats(sc_dpl_qc, &args->out_raw);
		}

		//assert(qc->big_t == big_t+1);
		big_t++;
	}
	
	free_sc_dpl_qc(sc_dpl_qc);
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
 * \brief Creates a new \ref sc_dpl_qc
 * 
 * \param[in] args The args structure as filled from the command line
 * \param[in] recipe The recipe the \ref sc_dpl_qc will use for computation 
 * \param[in] has_loss Whether or not the sc_dpl_qc should have loss
 */
SC_DPL_QC *create_sc_dpl_qc(ARGS *args, RECIPE_ADV *recipe, int has_loss) {
	SC_DPL_QC *sc_dpl_qc;
	DPL_QC *dpl_qc;
	QC *qc;
	int i, j, is, js, d, n;
	QUBIT ***q_arr, ***q_arr2;
	SET *bdy_s1_pr, *bdy_s2_pr, *bdy_s1_du, *bdy_s2_du; // spatial boundaries
	SET *bdy_t1_pr, *bdy_t2_pr, *bdy_t1_du, *bdy_t2_du; // temporal boundaries
	SYNDROME ***syn_arr_pr, ***syn_arr_du;
	SET ****set_arr_pr, ****set_arr_du;
	int num_meas_left_1, num_meas_left_0;
	SET *bdy;
	SET_POS *set_pos;

	sc_dpl_qc = (SC_DPL_QC *)my_malloc(sizeof(SC_DPL_QC));
	sc_dpl_qc->flush_file = args->flush_file;

	sc_dpl_qc->d = d = args->d;
	sc_dpl_qc->n = n = 2 * d - 1;
	sc_dpl_qc->ems = args->ems;

	sc_dpl_qc->dpl_qc = dpl_qc = dpl_create_dpl_qc(args->s0, args->s1, DE_HT_FACTOR * n * n, 
		STICK_HT_FACTOR * d * d, args->p, args->p_loss, args->p_depolar, args->t_delete, recipe, args->ems, has_loss);
	qc = dpl_qc->qc;
	
	sc_dpl_qc->frame = (int **)my_2d_calloc(n, n, sizeof(int));

	sc_dpl_qc->q_arr = q_arr = (QUBIT ***)my_2d_calloc(n, n, sizeof(QUBIT *));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			q_arr[i][j] = qc_create_and_insert_qubit(qc, i, j, 0, QUBIT_HT_SIZE);
		}
	}
	
	sc_dpl_qc->q_arr2 = q_arr2 = (QUBIT ***)my_2d_calloc(n, n, sizeof(QUBIT *));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			q_arr2[i][j] = qc_create_and_insert_qubit(qc, i, j, 1, QUBIT_HT_SIZE);
		}
	}

	// the location of boundaries is irrelevant, but a sensible choice is good for debugging
	sc_dpl_qc->bdy_s1_pr = bdy_s1_pr = qc_create_set_adv(qc, PRIMAL_BOUNDARY, -1, 0, 0, NULL);
	sc_dpl_qc->bdy_s2_pr = bdy_s2_pr = qc_create_set_adv(qc, PRIMAL_BOUNDARY, -2, 0, 0, NULL);
	sc_dpl_qc->bdy_t1_pr = bdy_t1_pr = qc_create_set_adv(qc, PRIMAL_BOUNDARY, -1, 1, 0, NULL);
	sc_dpl_qc->bdy_t2_pr = bdy_t2_pr = qc_create_set_adv(qc, PRIMAL_BOUNDARY, -2, 1, 0, NULL);

	sc_dpl_qc->bdy_s1_du = bdy_s1_du = qc_create_set_adv(qc, DUAL_BOUNDARY, -1, 0, 0, NULL);
	sc_dpl_qc->bdy_s2_du = bdy_s2_du = qc_create_set_adv(qc, DUAL_BOUNDARY, -2, 0, 0, NULL);
	sc_dpl_qc->bdy_t1_du = bdy_t1_du = qc_create_set_adv(qc, DUAL_BOUNDARY, -1, 1, 0, NULL);
	sc_dpl_qc->bdy_t2_du = bdy_t2_du = qc_create_set_adv(qc, DUAL_BOUNDARY, -2, 1, 0, NULL);

	sc_dpl_qc->boundaries = (BALL **)my_malloc(8 * sizeof(BALL *));
	sc_dpl_qc->boundaries[0] = bdy_s1_pr->ball;
	sc_dpl_qc->boundaries[1] = bdy_s2_pr->ball;
	sc_dpl_qc->boundaries[2] = bdy_t1_pr->ball;
	sc_dpl_qc->boundaries[3] = bdy_t2_pr->ball;
	sc_dpl_qc->boundaries[4] = bdy_s1_du->ball;
	sc_dpl_qc->boundaries[5] = bdy_s2_du->ball;
	sc_dpl_qc->boundaries[6] = bdy_t1_du->ball;
	sc_dpl_qc->boundaries[7] = bdy_t2_du->ball;

	qc->get_boundary = get_boundary;
	qc->boundaries = (void *)sc_dpl_qc->boundaries; 

	sc_dpl_qc->syn_arr_pr = syn_arr_pr = (SYNDROME ***)my_2d_calloc(d, d - 1, sizeof(SYNDROME *));
	sc_dpl_qc->set_arr_pr = set_arr_pr = (SET ****)my_3d_calloc(d, d - 1, 2, sizeof(SET *));
	for (is = 0; is < d; is++) {
		num_meas_left_1 = 1;
		if (is == 0 || is == d - 1) {
			num_meas_left_0 = 5;
		}
		else {
			num_meas_left_0 = 6;
		}

		for (js = 0; js < d - 1; js++) {
			// current
			set_arr_pr[is][js][1] = qc_create_set_adv(qc, PRIMAL, is << 1, (js << 1) + 1, num_meas_left_1, bdy_t1_pr);
			set_pos = qc_create_set_pos(is, js, 1);
			cdll_insert_head(set_arr_pr[is][js][1]->set_pos_cdll, set_pos, NULL);

			syn_arr_pr[is][js] = qc_create_syndrome();
			qc_insert_syndrome(qc, syn_arr_pr[is][js]);
			qc_associate_syndrome(set_arr_pr[is][js][1], syn_arr_pr[is][js]);
			qc_insert_set(qc, set_arr_pr[is][js][1]);

			// future
			set_arr_pr[is][js][0] = qc_create_set_adv(qc, PRIMAL, is << 1, (js << 1) + 1, num_meas_left_0, bdy_t1_pr);
			set_pos = qc_create_set_pos(is, js, 0);
			cdll_insert_head(set_arr_pr[is][js][0]->set_pos_cdll, set_pos, NULL);

			qc_insert_set(qc, set_arr_pr[is][js][0]);
		}
	}
	
	sc_dpl_qc->syn_arr_du = syn_arr_du = (SYNDROME ***)my_2d_calloc(d - 1, d, sizeof(SYNDROME *));
	sc_dpl_qc->set_arr_du = set_arr_du = (SET ****)my_3d_calloc(d - 1, d, 2, sizeof(SET *));
	for (is = 0; is < d - 1; is++) {
		if (is == 0) {
			bdy = bdy_s1_du;
		}
		else if (is == d - 2) {
			bdy = bdy_s2_du;
		}
		else {
			bdy = NULL;
		}

		for (js = 0; js < d; js++) {
			if (js == 0 || js == d - 1) {
				num_meas_left_1 = 4;
				num_meas_left_0 = 5;
			}
			else {
				num_meas_left_1 = 5;
				num_meas_left_0 = 6;
			}

			set_arr_du[is][js][1] = qc_create_set_adv(qc, DUAL, (is << 1) + 1, js << 1, num_meas_left_1, bdy);
			set_arr_du[is][js][0] = qc_create_set_adv(qc, DUAL, (is << 1) + 1, js << 1, num_meas_left_0, bdy);
			
			set_pos = qc_create_set_pos(is, js, 1);
			cdll_insert_head(set_arr_du[is][js][1]->set_pos_cdll, set_pos, NULL);
			set_pos = qc_create_set_pos(is, js, 0);
			cdll_insert_head(set_arr_du[is][js][0]->set_pos_cdll, set_pos, NULL);

			syn_arr_du[is][js] = qc_create_syndrome();
			qc_insert_syndrome(qc, syn_arr_du[is][js]);
			qc_associate_syndrome(set_arr_du[is][js][1], syn_arr_du[is][js]);
			qc_insert_set(qc, set_arr_du[is][js][1]);
			qc_insert_set(qc, set_arr_du[is][js][0]);
		}
	}
	
	sc_dpl_qc->t1 = double_time();

	sc_dpl_qc->num_checks = 0;
	sc_dpl_qc->last_X_check = 0;
	sc_dpl_qc->last_Z_check = 0;
	sc_dpl_qc->num_X_changes = 0;
	sc_dpl_qc->num_Z_changes = 0;
	sc_dpl_qc->bdy_meas_pr = 1;
	sc_dpl_qc->bdy_edges_pr = 0;
	sc_dpl_qc->bdy_meas_du = 1;
	sc_dpl_qc->bdy_edges_du = 0;

	return sc_dpl_qc;
}

/**
 * \brief Copies a \ref sc_dpl_qc
 * 
 * \param[in] sc_dpl_qc The \ref sc_dpl_qc to be copied
 *
 * \return The newly copied \ref sc_dpl_qc
 */
SC_DPL_QC *copy_sc_dpl_qc(SC_DPL_QC *sc_dpl_qc) {
	SC_DPL_QC *sc_dpl_qc2;
	int i, j;
	int is, js;
	int d, n;

	d = sc_dpl_qc->d;
	n = sc_dpl_qc->n;
	
	sc_dpl_qc2 = (SC_DPL_QC *)my_malloc(sizeof(SC_DPL_QC));
	
	// Copy values from sc_dpl_qc to sc_dpl_qc2
	*sc_dpl_qc2 = *sc_dpl_qc;

	// Copy the dpl_qc
	sc_dpl_qc2->dpl_qc = dpl_copy_dpl_qc(sc_dpl_qc->dpl_qc);

	// Copy the frame
	sc_dpl_qc2->frame = (int **)my_2d_calloc(n, n, sizeof(int));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			sc_dpl_qc2->frame[i][j] = sc_dpl_qc->frame[i][j];
		}
	}

	// Copy the qubits
	sc_dpl_qc2->q_arr = (QUBIT ***)my_2d_calloc(n, n, sizeof(QUBIT *));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			sc_dpl_qc2->q_arr[i][j] = sc_dpl_qc->q_arr[i][j]->copy;
		}
	}
	sc_dpl_qc2->q_arr2 = (QUBIT ***)my_2d_calloc(n, n, sizeof(QUBIT *));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			sc_dpl_qc2->q_arr2[i][j] = sc_dpl_qc->q_arr2[i][j]->copy;
		}
	}

	// Copy the PRIMAL syndrome and set arrays 
	sc_dpl_qc2->syn_arr_pr = (SYNDROME ***)my_2d_calloc(d, d - 1, sizeof(SYNDROME *));
	sc_dpl_qc2->set_arr_pr = (SET ****)my_3d_calloc(d, d - 1, 2, sizeof(SET *));
	for (is = 0; is < d; is++) {
		for (js = 0; js < d - 1; js++) {
			sc_dpl_qc2->syn_arr_pr[is][js] = sc_dpl_qc->syn_arr_pr[is][js]->copy;
			
			if (!qc_is_boundary_set(sc_dpl_qc->set_arr_pr[is][js][0])) {
				sc_dpl_qc2->set_arr_pr[is][js][0] = sc_dpl_qc->set_arr_pr[is][js][0]->copy;
			} else {
				// Do not copy boundaries. The copy of a boundary is NULL.
				sc_dpl_qc2->set_arr_pr[is][js][0] = sc_dpl_qc->set_arr_pr[is][js][0];
			}
			if (!qc_is_boundary_set(sc_dpl_qc->set_arr_pr[is][js][1])) {
				sc_dpl_qc2->set_arr_pr[is][js][1] = sc_dpl_qc->set_arr_pr[is][js][1]->copy;
			} else {
				// Do not copy boundaries. The copy of a boundary is NULL.
				sc_dpl_qc2->set_arr_pr[is][js][1] = sc_dpl_qc->set_arr_pr[is][js][1];
			}
		}
	}
	
	// Copy the DUAL syndrom and set arrays
	sc_dpl_qc2->syn_arr_du = (SYNDROME ***)my_2d_calloc(d - 1, d, sizeof(SYNDROME *));
	sc_dpl_qc2->set_arr_du = (SET ****)my_3d_calloc(d - 1, d, 2, sizeof(SET *));
	for (is = 0; is < d - 1; is++) {
		for (js = 0; js < d; js++) {
			sc_dpl_qc2->syn_arr_du[is][js] = sc_dpl_qc->syn_arr_du[is][js]->copy;

			if (!qc_is_boundary_set(sc_dpl_qc->set_arr_du[is][js][0])) {
				sc_dpl_qc2->set_arr_du[is][js][0] = sc_dpl_qc->set_arr_du[is][js][0]->copy;
			} else {
				// Do not copy boundaries. The copy of a boundary is NULL.
				sc_dpl_qc2->set_arr_du[is][js][0] = sc_dpl_qc->set_arr_du[is][js][0];
			}
			if (!qc_is_boundary_set(sc_dpl_qc->set_arr_du[is][js][1])) {
				sc_dpl_qc2->set_arr_du[is][js][1] = sc_dpl_qc->set_arr_du[is][js][1]->copy;
			} else {
				// Do not copy boundaries. The copy of a boundary is NULL.
				sc_dpl_qc2->set_arr_du[is][js][1] = sc_dpl_qc->set_arr_du[is][js][1];
			}
		}
	}

	return sc_dpl_qc2;
}

/**
 * \brief Frees a \ref sc_dpl_qc
 * 
 * \param[in] sc_dpl_qc The \ref sc_dpl_qc to be freed
 */
void free_sc_dpl_qc(SC_DPL_QC *sc_dpl_qc) {
	int is, js, d, n;

	d = sc_dpl_qc->d;
	n = sc_dpl_qc->n;

	// Free the \ref dpl_qc
	dpl_free_dpl_qc(sc_dpl_qc->dpl_qc);
	
	// Free the Pauli frame
	my_2d_free(n, (void **)sc_dpl_qc->frame);

	// Free the arrays of qubits (the qubits themselves are owned by the qc)
	my_2d_free(n, (void **)sc_dpl_qc->q_arr);
	my_2d_free(n, (void **)sc_dpl_qc->q_arr2);
	
	// Free the primal syndrome array and the syndromes
	for (is = 0; is < d; is++) {
		for (js = 0; js < d - 1; js++) {
			qc_free_syndrome(sc_dpl_qc->syn_arr_pr[is][js]);
		}
	}
	my_2d_free(d, (void **)sc_dpl_qc->syn_arr_pr);
	
	// Free the dual syndrome array and the syndromes
	for (is = 0; is < d - 1; is++) {
		for (js = 0; js < d; js++) {
			qc_free_syndrome(sc_dpl_qc->syn_arr_du[is][js]);
		}
	}
	my_2d_free(d - 1, (void **)sc_dpl_qc->syn_arr_du);

	// Free the primal and dual set arrays
	my_3d_free(d, d - 1, (void ***)sc_dpl_qc->set_arr_pr);
	my_3d_free(d - 1, d, (void ***)sc_dpl_qc->set_arr_du);

	// Free the sets associated with the boundaries

	qc_free_bdy(sc_dpl_qc->bdy_s1_pr);
	qc_free_bdy(sc_dpl_qc->bdy_s2_pr);
	qc_free_bdy(sc_dpl_qc->bdy_t1_pr);
	qc_free_bdy(sc_dpl_qc->bdy_t2_pr);
	qc_free_bdy(sc_dpl_qc->bdy_s1_du);
	qc_free_bdy(sc_dpl_qc->bdy_s2_du);
	qc_free_bdy(sc_dpl_qc->bdy_t1_du);
	qc_free_bdy(sc_dpl_qc->bdy_t2_du);

	// Free the boundaries array
	free(sc_dpl_qc->boundaries);

	free(sc_dpl_qc);
}

/**
 * \brief Frees a copied \ref sc_dpl_qc
 * 
 * \param[in] sc_dpl_qc The \ref sc_dpl_qc copy to be freed
 */
void free_sc_dpl_qc_copy(SC_DPL_QC *sc_dpl_qc) {
	/*
	 * \remark No need to free the boundaries, these are not copied
	 */
	int is, js, d, n;

	d = sc_dpl_qc->d;
	n = sc_dpl_qc->n;

	// Free the copied dpl_qc
	dpl_free_dpl_qc_copy(sc_dpl_qc->dpl_qc);

	// Free the Pauli frame and the qubit arrays (The qubits are stored in qc)
	my_2d_free(n, (void **)sc_dpl_qc->frame);
	my_2d_free(n, (void **)sc_dpl_qc->q_arr);
	my_2d_free(n, (void **)sc_dpl_qc->q_arr2);
	
	// Free the primal syndrome array and the syndromes
	for (is = 0; is < d; is++) {
		for (js = 0; js < d - 1; js++) {
			qc_free_syndrome(sc_dpl_qc->syn_arr_pr[is][js]);
		}
	}
	my_2d_free(d, (void **)sc_dpl_qc->syn_arr_pr);
	
	// Free the dual syndrome array and the syndromes
	for (is = 0; is < d - 1; is++) {
		for (js = 0; js < d; js++) {
			qc_free_syndrome(sc_dpl_qc->syn_arr_du[is][js]);
		}
	}
	my_2d_free(d - 1, (void **)sc_dpl_qc->syn_arr_du);

	// Free the primal and dual set arrays
	my_3d_free(d, d - 1, (void ***)sc_dpl_qc->set_arr_pr);
	my_3d_free(d - 1, d, (void ***)sc_dpl_qc->set_arr_du);

	free(sc_dpl_qc);
}

/**
 * \brief 
 * 
 * \param[in] 
 */

void measurement(SC_DPL_QC *sc_dpl_qc, int i, int j, QUBIT ***q_arr, SET ****set_arr, SET *set1, SET *set2) {
	DPL_QC *dpl_qc;
	SET *setk, *setr;
	int mres, prod, sets_are_eq;
	CDLL_NODE *node;
	SET_POS *sp;

	dpl_qc = sc_dpl_qc->dpl_qc;

	// If a merging occurs due to loss, there is no way to
	// determine which set was kept, so we need to determine this
	// before measuring the qubit.
	setk = (qc_is_boundary_set(set1) || set1->big_t > set2->big_t) ? set1 : set2;
	setr = (setk == set1) ? set2 : set1;
	
	sets_are_eq = (setk == setr) ? 1 : 0;

	dpl_H(dpl_qc, q_arr[i][j]);
	mres = dpl_meas_Z(dpl_qc, q_arr[i][j], set1, set2);

	if (mres != 1 && (setk->type == PRIMAL_BOUNDARY || setk->type == DUAL_BOUNDARY)) {
		prod = 1;

		if (mres == 0) {
			prod = setk->ball->mp;
			setk->ball->mp = 1;
		} else if (!sets_are_eq) {
			prod = mres;
		}

		if (prod != 1 && setk->i == -1) {
			if (setk->type == PRIMAL_BOUNDARY) {
				sc_dpl_qc->bdy_meas_pr *= prod;
			} else {
				sc_dpl_qc->bdy_meas_du *= prod;
			}
		}
	}

	// If the measurement result was a loss, then a set merging has
	// occured. Therefore, we will need to update the array of
	// sets to point to the newly merged set.
	if (mres == 0) {
		node = setk->set_pos_cdll->next;
		while (node != setk->set_pos_cdll) {
			sp = (SET_POS *)node->key;
			assert(setk != NULL);
			if (set_arr[sp->i][sp->j][sp->lay] == setr) {
				set_arr[sp->i][sp->j][sp->lay] = setk;
			}
			node = node->next;
		}
	}
}

/**
 * \brief 
 * 
 * \param[in] 
 */

void final_measurement(SC_DPL_QC *sc_dpl_qc, SET ****set_arr, SYNDROME ***syn_arr, QUBIT ***q_arr, int i, int j, int lay1, int lay2, int type, SET *bdy, int num_meas_left) {
	DPL_QC *dpl_qc;
	QC *qc;

	int mres, prod;
	int is;
	int js;
	int sets_are_eq;
	SET *set1, *set2, *setk, *setr;
	SET_POS *sp, *set_pos;
	CDLL_NODE *node, *n;

	dpl_qc = sc_dpl_qc->dpl_qc;
	qc = dpl_qc->qc;


	is = i >> 1;
	js = j >> 1;

	//printf("-- %d, %d, %d --\n", is, js, lay1);

	set1 = set_arr[is][js][lay1];
	set2 = set_arr[is][js][lay2];

	setk = (qc_is_boundary_set(set1) || set1->t > set2->t) ? set1 : set2;
	setr = (setk == set1) ? set2 : set1;
	
	sets_are_eq = (setk == setr) ? 1 : 0;

	/*
	printf("X Keep Set: "); qc_print_set(setk);
	printf("X Remove Set: "); qc_print_set(setr);
	//*/

	// Perform the measurement
	dpl_H(dpl_qc, q_arr[i][j]);
	mres = dpl_meas_Z(dpl_qc, q_arr[i][j], set1, set2);

	// This condition is only true with loss enabled
	if (mres != 1 && (setk->type == PRIMAL_BOUNDARY || setk->type == DUAL_BOUNDARY)) {
		prod = 1;

		if (mres == 0) {
			prod = setk->ball->mp;
			setk->ball->mp = 1;
		} else if (!sets_are_eq) {
			prod = mres;
		}

		if (prod != 1 && setk->i == -1) {
			//printf("FOOO: (%d, %d): %d, %d, %ld\n", i, j, mres, prod, sc_dpl_qc->bdy_meas_pr);
			//printf("\t"); qc_print_set(set1);
			//printf("\t"); qc_print_set(set2);

			if (setk->type == PRIMAL_BOUNDARY) {
				//printf("FPRIMAL\n");
				sc_dpl_qc->bdy_meas_pr *= prod;
			} else {
				//printf("FDUAL\n");
				sc_dpl_qc->bdy_meas_du *= prod;
			}
		}
	}

	// Update any of the set positions that point to the removed set to point
	// to the new set. Note that due to the asynchronous application of
	// measurements not all positions in set_arr will point to setr or setk,
	// they may have already advanced to a new set.
	if (mres == 0) {
		node = setk->set_pos_cdll->next;
		while (node != setk->set_pos_cdll) {
			sp = (SET_POS *)node->key;
			assert(setk != NULL);
			if (set_arr[sp->i][sp->j][sp->lay] == setr) {
				set_arr[sp->i][sp->j][sp->lay] = setk;
			}
			node = node->next;
		}
	}

	// We now need to move the syndrome from the bottom layer
	// (Layer 1) to the top layer (Layer 2). 
	//
	// If we have merged due to loss, and the first set is the boundary, then
	// we do not want to move the syndrome into the new set. This is due to the
	// new set no longer having a syndrome cdll, as it was merged into the
	// boundary's. 
	if (!(mres == 0 && qc_is_boundary_set(setk))) { 
		qc_unassociate_syndrome(syn_arr[is][js]);
		qc_associate_syndrome(set_arr[is][js][lay2], syn_arr[is][js]);
	}

	//validate_set_arrays(sc_dpl_qc->d, sc_dpl_qc->set_arr_pr, sc_dpl_qc->set_arr_du);
	//printf("XX: "); qc_print_set(set_arr[is][js][lay1]);
	//printf("QQ: "); qc_print_set(set_arr[is][js][lay2]);

	// We want to remove and free the set_pos that was being used by the set at
	// this current position. Without this, we would leak memory.
	n = set_arr[is][js][lay1]->set_pos_cdll->next;
	while (n != set_arr[is][js][lay1]->set_pos_cdll) {
		sp = (SET_POS *)n->key;
		n = n->next;
		if (sp->i == is && sp->j == js && sp->lay == lay1) {
			cdll_delete_node(n->prev, free);
		}
	}

	//printf("AA: "); qc_print_set(set_arr[is][js][lay1]);
	//printf("BB: "); qc_print_set(set_arr[is][js][lay2]);

	set_arr[is][js][lay1] = qc_create_set_adv(qc, type, i, j, num_meas_left, bdy);
	qc_insert_set(qc, set_arr[is][js][lay1]);
	set_pos = qc_create_set_pos(is, js, lay1);
	cdll_insert_head(set_arr[is][js][lay1]->set_pos_cdll, set_pos, NULL);

	//printf("CC: "); qc_print_set(set_arr[is][js][lay1]);
	//printf("DD: "); qc_print_set(set_arr[is][js][lay2]);

	// Reverse the order of the sets to ensure that the new set is in layer 2
	set1 = set_arr[is][js][lay1];
	set_arr[is][js][lay1] = set_arr[is][js][lay2];
	set_arr[is][js][lay2] = set1;
	qc_swap_set_layer(set_arr[is][js][lay1], lay1);
	qc_swap_set_layer(set_arr[is][js][lay2], lay2);
	
	set1 = set_arr[is][js][lay1];
	set2 = set_arr[is][js][lay2];

	//printf("YY: "); qc_print_set(set1);
	//printf("ZZ: "); qc_print_set(set2);

	// We want to delete any and all set_pos that have been put into a boundary
	// set, however we want to keep the head node of the cdll. Hence we do
	// this, rather than using cdll_free().
	if (qc_is_boundary_set(setk)) {
		n = setk->set_pos_cdll->next;
		while (n != setk->set_pos_cdll) {
			sp = (SET_POS *)n->key;
			n = n->next;
			cdll_delete_node(n->prev, free);
		}
	}

	//printf("\n");

	//validate_set_arrays(sc_dpl_qc->d, sc_dpl_qc->set_arr_pr, sc_dpl_qc->set_arr_du);
}

//*
void validate_set_arrays(int d, SET ****set_arr_pr, SET ****set_arr_du) {
	//*
	SET *temp_set;
	CDLL_NODE *temp_cdll;
	SET_POS *temp_set_pos;
	int i, j, lay, found;

	printf("Validate Set Arrays:\n");
	for (i = 0; i < d; i++) {
		for (j = 0; j < d - 1; j++) {
			for (lay = 0; lay < 1; lay++) {
				temp_set = set_arr_pr[i][j][lay];
				if (temp_set) {
					temp_cdll = temp_set->set_pos_cdll->next;
					found = FALSE;
					while (temp_cdll != temp_set->set_pos_cdll) {
						temp_set_pos = (SET_POS *)temp_cdll->key;
						if (temp_set_pos->i == i && temp_set_pos->j == j && temp_set_pos->lay == lay) {
							found = TRUE;
						}
						temp_cdll = temp_cdll->next;
					}
					//printf("%2d,%d ", temp_set->i, lay);
					if (!qc_is_boundary_set(temp_set) && !found) {
						qc_print_set(temp_set);
						assert(found);
						printf("\n\n\n\n\n\n");
					}
				}
			}
			for (lay = 1; lay < 2; lay++) {
				temp_set = set_arr_pr[i][j][lay];
				if (temp_set) {
					temp_cdll = temp_set->set_pos_cdll->next;
					found = FALSE;
					while (temp_cdll != temp_set->set_pos_cdll) {
						temp_set_pos = (SET_POS *)temp_cdll->key;
						if (temp_set_pos->i == i && temp_set_pos->j == j && temp_set_pos->lay == lay) {
							found = TRUE;
						}
						temp_cdll = temp_cdll->next;
					}
					//printf("%2d,%d ", temp_set->i, lay);
					if (!qc_is_boundary_set(temp_set) && !found) {
						qc_print_set(temp_set);
						assert(found);
						printf("\n\n\n\n\n\n");
					}
				}
			}
		}
		//printf("\n");
	}

	for (i = 0; i < d - 1; i++) {
		for (j = 0; j < d ; j++) {
			for (lay = 0; lay < 2; lay++) {
				temp_set = set_arr_du[i][j][lay];
				if (temp_set) {
					temp_cdll = temp_set->set_pos_cdll->next;
					found = FALSE;
					while (temp_cdll != temp_set->set_pos_cdll) {
						temp_set_pos = (SET_POS *)temp_cdll->key;
						if (temp_set_pos->i == i && temp_set_pos->j == j && temp_set_pos->lay == lay) {
							found = TRUE;
						}
						temp_cdll = temp_cdll->next;
					}
					if (!qc_is_boundary_set(temp_set) && !found) {
						qc_print_set(temp_set);
						assert(found);
						printf("\n\n\n\n\n\n");
					}
				}
			}
		}
	}
	//*/
}
//*/

/**
 * \brief Measures the stabilizers of the topological cluster state
 *
 *
 * The timing found here is outlined in the paper _Topological Code Autotune_, 
* DOI: 10.1103/PhysRevX.2.041003, arXiv: http://arxiv.org/abs/1202.6111
* It details the simulation of a topological 3D cluster-state using only
* two layers of qubits with nearest neighbour interactions. There two
* layers can physically be interwoven, however for simplicity of
* code and simulation they are kept seperate. 

* In this simulation Hadamard is considered to take 0 time, essentially
* allowing initialisation into the + state. This is done to retain a
* simpler 6 round measurement process, while having no effect on the
* logical error rate obtained.

* The association of qubits as either primal or dual is arbitary. 

* The pattern of qubits chosen to have actions performed in each layer
* during each round was chosen so that each qubit is only ever acted on by
* one other qubit at once or is being initialised or measured. These
* patterns have been labelled A, B, C to try and make things a little
* clearer. 
 * 
 * \param[in] sc_dpl_qc The \ref sc_dpl_qc
 * \param[in] big_t The current big_t of simulation 
 */
void measure_stabilizers(SC_DPL_QC *sc_dpl_qc) {
	int i, j, is, js, lay1_pr, lay1_du, lay2_pr, lay2_du, d, n, inc;
	long int t;
	DPL_QC *dpl_qc;
	QC *qc;
	QUBIT ***q_arr, ***q_arr2;
	SET *bdy_s1_pr, *bdy_s2_pr; // primal boundaries
	SET *bdy_s1_du, *bdy_s2_du; // dual boundaries
	SYNDROME ***syn_arr_pr, ***syn_arr_du;
	SET ****set_arr_pr, ****set_arr_du;
	SET *set1, *set2;
	int ***mts;
	int num_meas_left;
	SET *bdy;

	d = sc_dpl_qc->d;
	n = sc_dpl_qc->n;

	mts = (int ***)my_3d_calloc(n, n, 2, sizeof(int));

	dpl_qc = sc_dpl_qc->dpl_qc;
	qc = dpl_qc->qc;

	// Setup nice variable names
	q_arr = sc_dpl_qc->q_arr;
	q_arr2 = sc_dpl_qc->q_arr2;
	bdy_s1_pr = sc_dpl_qc->bdy_s1_pr;
	bdy_s2_pr = sc_dpl_qc->bdy_s2_pr;
	bdy_s1_du = sc_dpl_qc->bdy_s1_du;
	bdy_s2_du = sc_dpl_qc->bdy_s2_du;
	syn_arr_pr = sc_dpl_qc->syn_arr_pr;
	syn_arr_du = sc_dpl_qc->syn_arr_du;
	set_arr_pr = sc_dpl_qc->set_arr_pr;
	set_arr_du = sc_dpl_qc->set_arr_du;

	//validate_set_arrays(d, set_arr_pr, set_arr_du);

	lay1_pr = lay1_du = 1;
	lay2_pr = lay2_du = 0;

	// Get the duration of the identity gate of the same length as the cZ gate,
	// and increment the time of the simulation.
	t = q_arr[0][0]->t;
	inc = qc->ems[dpl_qc->iden_cZ_id]->duration;
	t += inc;

	//
	// Round 1
	//
	#ifdef PRINT_ROUNDS
	printf("Round 1\n");
	#endif
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			//
			// Layer 1
			//

			// A
			if ((i % 4 == 0 && j % 4 == 1) || 
				(i % 4 == 2 && j % 4 == 3)) {
				// Perform a cZ gate between the qubit and the qubit to the
				// south. For qubits with no southern neighbour, perform the
				// identity gate.
				if (i < n - 1) {
					dpl_cZ(dpl_qc, q_arr[i][j], q_arr[i + 1][j]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr[i][j]);
				}
			}

			// B
			else if ((i % 4 == 0 && j % 4 == 3) || 
					 (i % 4 == 2 && j % 4 == 1)) {
				assert(j > 0);
				// Perform a cZ gate between the qubit and the qubit to the
				// west. These qubits will always have a western neighbour.
				dpl_cZ(dpl_qc, q_arr[i][j], q_arr[i][j - 1]);
			}
			
			// C (Dual Wall Measurement):
			else if ((i + j) % 4 == 0) {
				is = i >> 1;
				js = j >> 1;

				// We must define the two sets which the qubit being
				// measured belongs to. 

				// If i is odd, then the measurement will be part of two
				// cluster states sharing an eastern-western wall.
				if (i % 2 == 1) {
					set1 = set_arr_du[is][js][lay1_du];
					set2 = set_arr_du[is][js + 1][lay1_du];
				}

				// If we are on the northern boundary, then the qubit is
				// connected to the first dual boundary.
				else if (is == 0) {
					set1 = bdy_s1_du;
					set2 = set_arr_du[is][js][lay1_du];
				}

				// If we are on the southern boundary, then the qubit is
				// connected to the second dual boundary.
				else if (is == d - 1) {
					set1 = set_arr_du[is - 1][js][lay1_du];
					set2 = bdy_s2_du;
				}

				// Otherwise, the measurement is part of two cluster states
				// sharing a northen-southern wall.
				else {
					set1 = set_arr_du[is - 1][js][lay1_du];
					set2 = set_arr_du[is][js][lay1_du];
				}

				// Dual Measurement
				//printf("Round 1 Dual Meas\n");
				//qc_print_qubit(q_arr[i][j]);
				measurement(sc_dpl_qc, i, j, q_arr, set_arr_du, set1, set2);
				//printf("End: Round 1 Dual Meas\n");
			}

			//
			// Layer 2
			//
		
			// A
			if ((i + j) % 4 == 2) {
				// Initialise the qubits 
				dpl_init_Z(dpl_qc, q_arr2[i][j]);
				dpl_H(dpl_qc, q_arr2[i][j]);
			}

			// B
			else if ((i % 4 == 1 && j % 4 == 2) || 
					 (i % 4 == 3 && j % 4 == 0)) {
				// Perform a cZ gate between the qubit and the qubit to the
				// east. For qubits with no eastern neighbour, perform the
				// identity gate.
				if (j < n - 1) {
					dpl_cZ(dpl_qc, q_arr2[i][j], q_arr2[i][j + 1]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr2[i][j]);
				}
			}
			
			// C
			else if ((i % 4 == 1 && j % 4 == 0) || 
					 (i % 4 == 3 && j % 4 == 2)) {
				assert(i > 0);
				// Perform a cZ gate between the qubit and the qubit to the
				// west. These qubits will always have a northern neighbour.
				dpl_cZ(dpl_qc, q_arr2[i][j], q_arr2[i - 1][j]);
			}
		}
	}



	//validate_set_arrays(d, set_arr_pr, set_arr_du);
	
	// Advance the boundary, performing identity gates on remaining qubits.
	// Increment the time by inc
	advance_boundary(sc_dpl_qc, t);
	t += inc;

	//
	// Round 2
	//
	#ifdef PRINT_ROUNDS
	printf("Round 2\n");
	#endif
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			//
			// Layer 1
			// 
			
			// A* (Final Primal Measurement)
			if (i % 2 == 0 && j % 2 == 1) {
				if (j == 1) {
					bdy = bdy_s1_pr;
				}
				else if (j == n - 2) {
					bdy = bdy_s2_pr;
				}
				else {
					bdy = NULL;
				}
				num_meas_left = (i == 0 || i == n - 1) ? 5 : 6;
				//printf("Start: Round 2 Final Measurement Primal\n");
				final_measurement(sc_dpl_qc, set_arr_pr, syn_arr_pr, q_arr, i, j, lay1_pr, lay2_pr, PRIMAL, bdy, num_meas_left);
				//printf("End: Round 2 Final Measurement Primal\n");
			}
			
			// B*
			else if ((i + j) % 4 == 2) {
				dpl_cZ(dpl_qc, q_arr[i][j], q_arr2[i][j]);
			}
			
			// C
			else if ((i + j) % 4 == 0) {
				dpl_init_Z(dpl_qc, q_arr[i][j]);
				dpl_H(dpl_qc, q_arr[i][j]);
			}

			//
			// Layer 2
			//

			// B
			if ((i % 4 == 1 && j % 4 == 2) || 
				(i % 4 == 3 && j % 4 == 0)) {
				if (i < n - 1) {
					dpl_cZ(dpl_qc, q_arr2[i][j], q_arr2[i + 1][j]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr2[i][j]);
				}
			}
			
			// C
			else if ((i % 4 == 1 && j % 4 == 0) ||
					 (i % 4 == 3 && j % 4 == 2)) {
				if (j > 0) {
					dpl_cZ(dpl_qc, q_arr2[i][j], q_arr2[i][j - 1]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr2[i][j]);
				}
			}
		}
	}
	//validate_set_arrays(d, set_arr_pr, set_arr_du);

	advance_boundary(sc_dpl_qc, t);
	t += inc;

	// Round 3
	#ifdef PRINT_ROUNDS
	printf("Round 3\n");
	#endif
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			//
			// Layer 1
			//
			
			// A*
			if (i % 2 == 0 && j % 2 == 1) {
				dpl_init_Z(dpl_qc, q_arr[i][j]);
				dpl_H(dpl_qc, q_arr[i][j]);
			}

			// B*
			else if ((i + j) % 4 == 2) {
				is = i >> 1;
				js = j >> 1;

				if (i % 2 == 1) {
					set1 = set_arr_du[is][js][lay1_du];
					set2 = set_arr_du[is][js + 1][lay1_du];
				}
				else if (is == 0) {
					set1 = bdy_s1_du;
					set2 = set_arr_du[is][js][lay1_du];
				}
				else if (is == d - 1) {
					set1 = set_arr_du[is - 1][js][lay1_du];
					set2 = bdy_s2_du;
				}
				else {
					set1 = set_arr_du[is - 1][js][lay1_du];
					set2 = set_arr_du[is][js][lay1_du];
				}

				// Dual Measurement
				//printf("Start: Round 3 Dual Measurement\n");
				//qc_print_qubit(q_arr[i][j]);
				measurement(sc_dpl_qc, i, j, q_arr, set_arr_du, set1, set2);
				//printf("End: Round 3 Dual Measurement\n");
			}
			// C
			else if ((i + j) % 4 == 0) {
				dpl_cZ(dpl_qc, q_arr[i][j], q_arr2[i][j]);
			}

			//
			// Layer 2
			//

			// B
			if ((i % 4 == 1 && j % 4 == 2) || 
				(i % 4 == 3 && j % 4 == 0)) {
				assert(i > 0);
				dpl_cZ(dpl_qc, q_arr2[i][j], q_arr2[i - 1][j]);
			}

			// C
			else if ((i % 4 == 1 && j % 4 == 0) || 
					 (i % 4 == 3 && j % 4 == 2)) {
				if (j < n - 1) {
					dpl_cZ(dpl_qc, q_arr2[i][j], q_arr2[i][j + 1]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr2[i][j]);
				}
			}
		}
	}
	//validate_set_arrays(d, set_arr_pr, set_arr_du);

	advance_boundary(sc_dpl_qc, t);
	t += inc;

	// 
	// Round 4
	//
	#ifdef PRINT_ROUNDS
	printf("Round 4\n");
	#endif
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			// 
			// Layer 1
			// 

			// A
			if ((i % 4 == 0 && j % 4 == 1) || 
				(i % 4 == 2 && j % 4 == 3)) {
				if (i > 0) {
					dpl_cZ(dpl_qc, q_arr[i][j], q_arr[i-1][j]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr[i][j]);
				}
			}

			// B*
			else if ((i + j) % 4 == 2) {
				dpl_init_Z(dpl_qc, q_arr[i][j]);
				dpl_H(dpl_qc, q_arr[i][j]);
			}

			// C*
			else if ((i % 4 == 2 && j % 4 == 1) || 
					 (i % 4 == 0 && j % 4 == 3)) {
				if (j < n - 1) {
					dpl_cZ(dpl_qc, q_arr[i][j], q_arr[i][j + 1]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr[i][j]);
				}
			}

			//
			// Layer 2
			//
			
			// A (Primal Wall Measurement):
			if ((i + j) % 4 == 0) {
				is = i >> 1;
				js = j >> 1;

				if (j % 2 == 1) {
					set1 = set_arr_pr[is][js][lay1_pr];
					set2 = set_arr_pr[is + 1][js][lay1_pr];
				}
				else if (js == 0) {
					set1 = bdy_s1_pr;
					set2 = set_arr_pr[is][js][lay1_pr];
				}
				else if (js == d - 1) {
					set1 = set_arr_pr[is][js - 1][lay1_pr];
					set2 = bdy_s2_pr;
				}
				else {
					set1 = set_arr_pr[is][js - 1][lay1_pr];
					set2 = set_arr_pr[is][js][lay1_pr];
				}

				// Primal Measurement
				//printf("Start: Round 4 Primal Measurement\n");
				measurement(sc_dpl_qc, i, j, q_arr2, set_arr_pr, set1, set2);
				//printf("End: Round 4 Primal Measurement\n");
			}

			// B
			else if ((i % 4 == 1 && j % 4 == 0) || 
					 (i % 4 == 3 && j % 4 == 2)) {
				if (i < n - 1) {
					dpl_cZ(dpl_qc, q_arr2[i][j], q_arr2[i + 1][j]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr2[i][j]);
				}
			}

			// C
			else if ((i % 4 == 1 && j % 4 == 2) ||
					 (i % 4 == 3 && j % 4 == 0)) {
				if (j > 0) {
					dpl_cZ(dpl_qc, q_arr2[i][j], q_arr2[i][j - 1]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr2[i][j]);
				}
			}
		}
	}
	//validate_set_arrays(d, set_arr_pr, set_arr_du);

	advance_boundary(sc_dpl_qc, t);
	t += inc;

	//
	// Round 5
	//
	#ifdef PRINT_ROUNDS
	printf("Round 5\n");
	#endif
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			//
			// Layer 1
			//

			// A
			if ((i % 4 == 0 && j % 4 == 1) ||
					 (i % 4 == 2 && j % 4 == 3)) {
				if (j > 0) {
					dpl_cZ(dpl_qc, q_arr[i][j], q_arr[i][j - 1]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr[i][j]);
				}
			}

			// B
			else if ((i % 4 == 0 && j % 4 == 3) || 
					 (i % 4 == 2 && j % 4 == 1)) {
				if (i < n - 1) {
					dpl_cZ(dpl_qc, q_arr[i][j], q_arr[i + 1][j]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr[i][j]);
				}
			}

			// C
			else if ((i + j) % 4 == 2) {
				dpl_cZ(dpl_qc, q_arr[i][j], q_arr2[i][j]);
			}

			//
			// Layer 2
			//

			// A* (Final Dual Cell Measurement)
			if (i % 2 == 1 && j % 2 == 0) {
				is = i >> 1;
				js = j >> 1;

				if (is == 0) {
					bdy = bdy_s1_du;
				}
				else if (is == d - 2) {
					bdy = bdy_s2_du;
				}
				else {
					bdy = NULL;
				}
				num_meas_left = (js == 0 || js == d-1) ? 5 : 6;

				//printf("Start: Round 5 Final Measurement Dual\n");
				final_measurement(sc_dpl_qc, set_arr_du, syn_arr_du, q_arr2, i, j, lay1_du, lay2_du, DUAL, bdy, num_meas_left);
				//printf("End: Round 5 Final Measurement Dual\n");
			}

			// A
			else if ((i + j) % 4 == 0) {
				dpl_init_Z(dpl_qc, q_arr2[i][j]);
				dpl_H(dpl_qc, q_arr2[i][j]);
			}
		}
	}
	//validate_set_arrays(d, set_arr_pr, set_arr_du);

	advance_boundary(sc_dpl_qc, t);
	t += inc;

	//
	// Round 6
	//
	#ifdef PRINT_ROUNDS
	printf("Round 6\n");
	#endif
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			//
			// Layer 1
			//
			
			// A
			if ((i % 4 == 0 && j % 4 == 1) ||
				(i % 4 == 2 && j % 4 == 3)) {
				if (j < n - 1) {
					dpl_cZ(dpl_qc, q_arr[i][j], q_arr[i][j + 1]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr[i][j]);
				}
			}

			// B
			else if ((i % 4 == 2 && j % 4 == 1) ||
					 (i % 4 == 0 && j % 4 == 3)) {
				if (i > 0) {
					dpl_cZ(dpl_qc, q_arr[i][j], q_arr[i - 1][j]);
				}
				else {
					dpl_iden_cZ(dpl_qc, q_arr[i][j]);
				}
			}

			// C
			else if ((i + j) % 4 == 0) {
				dpl_cZ(dpl_qc, q_arr[i][j], q_arr2[i][j]);
			}

			//
			// Layer 2
			//

			// A*
			if (i % 2 == 1 && j % 2 == 0) {
				dpl_init_Z(dpl_qc, q_arr2[i][j]);
				dpl_H(dpl_qc, q_arr2[i][j]);
			}

			// A (Primal Wall Measurement)
			else if ((i + j) % 4 == 2) {
				is = i >> 1;
				js = j >> 1;

				if (j % 2 == 1) {
					set1 = set_arr_pr[is][js][lay1_pr];
					set2 = set_arr_pr[is + 1][js][lay1_pr];
				}
				else if (js == 0) {
					set1 = bdy_s1_pr;
					set2 = set_arr_pr[is][js][lay1_pr];
				}
				else if (js == d - 1) {
					set1 = set_arr_pr[is][js - 1][lay1_pr];
					set2 = bdy_s2_pr;
				}
				else {
					set1 = set_arr_pr[is][js - 1][lay1_pr];
					set2 = set_arr_pr[is][js][lay1_pr];
				}

				// Primal Measurement
				//printf("Start: Round 6 Primal Measurement\n");
				measurement(sc_dpl_qc, i, j, q_arr2, set_arr_pr, set1, set2);
				//printf("End: Round 6 Primal Measurement\n");
			}
		}
	}
	//validate_set_arrays(d, set_arr_pr, set_arr_du);

	advance_boundary(sc_dpl_qc, t);
	
	//print_all_es(sc_dpl_qc->n, sc_dpl_qc->q_arr);
	//print_all_es(sc_dpl_qc->n, sc_dpl_qc->q_arr2);

	my_3d_free(n, n, (void ***)mts);
}

/**
 * \brief Advances the boundary
 * 
 * \param[in] sc_dpl_qc The \ref sc_dpl_qc
 * \param[in] t The current t 
 */
void advance_boundary(SC_DPL_QC *sc_dpl_qc, long int t) {
	int i, j, n;
	DPL_QC *dpl_qc;
	QUBIT ***q_arr, ***q_arr2;

	n = sc_dpl_qc->n;
	dpl_qc = sc_dpl_qc->dpl_qc;
	q_arr = sc_dpl_qc->q_arr;
	q_arr2 = sc_dpl_qc->q_arr2;

	for (j = 0; j < n; j++) {
		if (q_arr[0][j]->t < t) {
			dpl_iden_cZ(dpl_qc, q_arr[0][j]);
			assert(q_arr[0][j]->t == t);
		}
		if (q_arr[n - 1][j]->t < t) {
			dpl_iden_cZ(dpl_qc, q_arr[n - 1][j]);
			assert(q_arr[n - 1][j]->t == t);
		}
	}

	for (j = 0; j < n; j += 2) {
		if (q_arr2[0][j]->t < t) {
			dpl_iden_cZ(dpl_qc, q_arr2[0][j]);
			assert(q_arr2[0][j]->t == t);
		}

		if (q_arr2[n - 1][j]->t < t) {
			dpl_iden_cZ(dpl_qc, q_arr2[n - 1][j]);
			assert(q_arr2[n - 1][j]->t == t);
		}
	}

	for (i = 0; i < n; i += 2) {
		if (q_arr[i][0]->t < t) {
			dpl_iden_cZ(dpl_qc, q_arr[i][0]);
			assert(q_arr[i][0]->t == t);
		}

		if (q_arr[i][n - 1]->t < t) {
			dpl_iden_cZ(dpl_qc, q_arr[i][n - 1]);
			assert(q_arr[i][n - 1]->t == t);
		}
	}

	for (i = 0; i < n; i++) {
		if (q_arr2[i][0]->t < t) {
			dpl_iden_cZ(dpl_qc, q_arr2[i][0]);
			assert(q_arr2[i][0]->t == t);
		}

		if (q_arr2[i][n - 1]->t < t) {
			dpl_iden_cZ(dpl_qc, q_arr2[i][n - 1]);
			assert(q_arr2[i][n - 1]->t == t);
		}
	}
}

/**
 * \brief Processes an augmented edge 
 * 
 * \param[in] ae The \ref ae to process
 * \param[in] sc_dpl_qc The \ref sc_dpl_qc to augment the edge on
 * \param[in] type The type of edge (PRIMAL, DUAL) 
 */
void process_aug_edge(AUG_EDGE *ae, SC_DPL_QC *sc_dpl_qc, int type) {
	VERTEX *v;
	v = ae->vb;

	if (v->v_num != -1) {
		return;
	}

	if (type == PRIMAL) {
		sc_dpl_qc->bdy_edges_pr++;
	} else {
		sc_dpl_qc->bdy_edges_du++;
	}
}

/**
 * \brief Process augmented edges from a \ref matching
 * 
 * \param[in] m The \ref matching containing the augmented edges
 * \param[in] sc_dpl_qc The \ref sc_dpl_qc to augment the edges on
 * \param[in] type The type of edges (PRIMAL, DUAL) 
 */
void process_aug_edges(MATCHING *m, SC_DPL_QC *sc_dpl_qc, int type) {
	AUG_EDGE *ae;

	ae = m_get_aug_edge(m);
	while (ae != NULL) {
		process_aug_edge(ae, sc_dpl_qc, type);
		m_delete_aug_edge(ae);
		ae = m_get_aug_edge(m);
	}
}

/**
 * \brief Corrects the \ref measurements in a ref sc_dpl_qc
 * 
 * Processes the augmented edges from the primal and dual matchings
 *
 * \param[in] sc_dpl_qc The \ref sc_dpl_qc to correct 
 */
void correct_mts(SC_DPL_QC *sc_dpl_qc) {
	QC *qc;

	qc = sc_dpl_qc->dpl_qc->qc;

	process_aug_edges(qc->m_pr, sc_dpl_qc, PRIMAL);
	process_aug_edges(qc->m_du, sc_dpl_qc, DUAL);
}

/**
 * \brief Performs three rounds of perfect stabilizer measurement, then
 * determines if there has been a logical X and/or Z error.
 * 
 * \param[in] sc_dpl_qc The \ref sc_dpl_qc to test
 * \param[out] out A file pointer to where the output stats of the test should
 * be sent.
 */
void test_correct(SC_DPL_QC *sc_dpl_qc, FILE **out) {
	SC_DPL_QC *sc_dpl_qc2;
	DPL_QC *dpl_qc2;
	QC *qc2;
	int count, t_delay;

	//printf("\ntc: %ld\n", sc_dpl_qc->dpl_qc->qc->big_t);

	sc_dpl_qc2 = copy_sc_dpl_qc(sc_dpl_qc);
	dpl_qc2 = sc_dpl_qc2->dpl_qc;
	qc2 = dpl_qc2->qc;
	dpl_save_rng();

	qc2->perfect_gates = true;
	dpl_qc2->has_loss =	false;

	t_delay = qc2->m_pr->t_delay;
	qc2->m_pr->undo_flag = TRUE;
	qc2->m_du->undo_flag = TRUE;
	qc2->m_pr->t_delay = 0;
	qc2->m_du->t_delay = 0;

	if (qc2->track) {
		measure_stabilizers(sc_dpl_qc2);
		measure_stabilizers(sc_dpl_qc2);
		qc_finalize_nests(qc2, qc2->big_t);
		qc_convert_nests(qc2, TRUE);
	} else {
		measure_stabilizers(sc_dpl_qc2);
		measure_stabilizers(sc_dpl_qc2);

		// The recipe has already finalized the nest by its own nature, however
		// this will update qc->unfinalized_big_t to current time so that
		// qc_convert_nests will insert all dots and lines into matching. 
		qc_finalize_nests(qc2, qc2->big_t);

		// The recipe has already finalized the dots and lines at this stage,
		// this will simply insert them into the matching.  
		qc_convert_nests(qc2, TRUE);
	}

	qc2->big_t += 2;

	qc_mwpm(qc2, TRUE);
	correct_mts(sc_dpl_qc2);

	sc_dpl_qc->num_checks++;

	if (sc_dpl_qc2->bdy_meas_pr == 1) {
		count = sc_dpl_qc2->bdy_edges_pr;
	}
	else {
		count = sc_dpl_qc2->bdy_edges_pr + 1;
	}
	if (qc2->logical_loss_err_pr != qc2->last_logical_loss_err_pr) {
		count += qc2->logical_loss_err_pr - qc2->last_logical_loss_err_pr;
		sc_dpl_qc->dpl_qc->qc->logical_loss_err_pr = qc2->logical_loss_err_pr;
		sc_dpl_qc->dpl_qc->qc->last_logical_loss_err_pr = qc2->logical_loss_err_pr;
	}
	if (count % 2 != sc_dpl_qc->last_Z_check) {
		sc_dpl_qc->last_Z_check = count%2;
		sc_dpl_qc->num_Z_changes++;
		print_stats(sc_dpl_qc, out);
	}

	if (sc_dpl_qc2->bdy_meas_du == 1) {
		count = sc_dpl_qc2->bdy_edges_du;
	}
	else {
		count = sc_dpl_qc2->bdy_edges_du + 1;
	}
	if (qc2->logical_loss_err_du != qc2->last_logical_loss_err_du) {
		count += qc2->logical_loss_err_du - qc2->last_logical_loss_err_du;
		sc_dpl_qc->dpl_qc->qc->logical_loss_err_du = qc2->logical_loss_err_du;
		sc_dpl_qc->dpl_qc->qc->last_logical_loss_err_du = qc2->logical_loss_err_du;
	}
	if (count % 2 != sc_dpl_qc->last_X_check) {
		sc_dpl_qc->last_X_check = count%2;
		sc_dpl_qc->num_X_changes++;
		print_stats(sc_dpl_qc, out);
	}

	qc_undo_mwpm(qc2);
	free_sc_dpl_qc_copy(sc_dpl_qc2);

	sc_dpl_qc->dpl_qc->qc->m_pr->undo_flag = FALSE;
	sc_dpl_qc->dpl_qc->qc->m_du->undo_flag = FALSE;
	sc_dpl_qc->dpl_qc->qc->m_pr->t_delay = t_delay;
	sc_dpl_qc->dpl_qc->qc->m_du->t_delay = t_delay;
	dpl_restore_rng();

	//printf("etc\n\n");
}

/**
 * \brief Outputs the current statistics of the simulation including: Current
 * time, current big_t, the parity of X and Z in the last check, the number of
 * checks, and the number of X and Z logical changes.
 * 
 * \param[in] sc_dpl_qc The \ref sc_dpl_qc to print the statistics of
 * \param[out] out A file pointer to where the output stats of the test should
 * be sent.
 */
void print_stats(SC_DPL_QC *sc_dpl_qc, FILE **out) {
	double t2;

	t2 = double_time();
	fprintf(*out, ">>> %.2f | %ld | Last: %d %d | Checks: %d %d | Loss: %d %d | Changes: %d %d\n", 
		t2 - sc_dpl_qc->t1, 
		sc_dpl_qc->dpl_qc->qc->big_t, 
		sc_dpl_qc->last_X_check, 
		sc_dpl_qc->last_Z_check, 
		sc_dpl_qc->num_checks, 
		sc_dpl_qc->num_checks, 
		sc_dpl_qc->dpl_qc->qc->logical_loss_err_pr,
		sc_dpl_qc->dpl_qc->qc->logical_loss_err_du,
		sc_dpl_qc->num_X_changes, 
		sc_dpl_qc->num_Z_changes
	);

	if (*out != stdout && sc_dpl_qc->flush_file) {
		fflush(*out);
		fclose(*out);
		*out = (FILE *)fopen("out_raw", "a");
	}
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
