// examples of how to use the Polyestimate library

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "polyestimate/polyestimate.h"
#include "my_time/my_time.h"

int main(int argc, char **argv) {
	int i; // loop variable
	int d; // code distance
	int type; // either X or Z depending on the desired type of logical error rate
	char fname[MAX_STRING_LENGTH]; // path to and name of database file
	char p_dir[MAX_STRING_LENGTH]; // path to directory containing depolarizing error rates in a file p.json
	char gates_dir[MAX_STRING_LENGTH]; // path to directory containing detailed gate error models in a file gates.json
	char em_dir[MAX_STRING_LENGTH]; // path to directory containing detailed gate error models in files of the form *_em
	DATABASE *db; // structure containing the loaded database
	double p_cnot, p_H, p_init_Z, p_meas_Z, p_iden_cnot, p_iden_H, p_iden_init_Z, p_iden_meas_Z; // individual gate error rates
	double p; // global gate error rate, overrides individual data error rates if set to positive value
	double p_L; // logical error rate
	double p_arr[8]; // array containing all of the individual data error rates
	FORM *f; // structure containing formulae

	// default path to and name of database file
	strcpy(fname, "../data/stats.txt");

	// default data directories
	strcpy(p_dir, "../data/p");
	strcpy(gates_dir, "../data/gates");
	strcpy(em_dir, "../data/em");

	// default logical error type
	type = X;

	// default code distance, can be changed here or by a command line argument
	d = 7;

	// default gate error rates, user can change these default values if they only desire a specific logical error rate or
	// alternatively can change these values using command line arguments or by loading a file with the desired values as shown
	// in example 3
	p = -1;
	p_cnot = 0.001;
	p_H = 0.001;
	p_init_Z = 0.001;
	p_meas_Z = 0.001;
	p_iden_cnot = 0.001;
	p_iden_H = 0.001;
	p_iden_init_Z = 0.001;
	p_iden_meas_Z = 0.001;

	for (i=1; i<argc; i++) {
      if (!strcmp(argv[i], "-fname")) {
			strcpy(fname, argv[++i]);
      }
      else if (!strcmp(argv[i], "-p_dir")) {
			strcpy(p_dir, argv[++i]);
      }
      else if (!strcmp(argv[i], "-gates_dir")) {
			strcpy(gates_dir, argv[++i]);
      }
      else if (!strcmp(argv[i], "-em_dir")) {
			strcpy(em_dir, argv[++i]);
      }
		else if (!strcmp(argv[i], "-d")) {
			d = atoi(argv[++i]);
		}
		else if (!strcmp(argv[i], "-p")) {
			p = atof(argv[++i]);
		}
		else if (!strcmp(argv[i], "-p_cnot")) {
			p_cnot = atof(argv[++i]);
		}
		else if (!strcmp(argv[i], "-p_H")) {
			p_H = atof(argv[++i]);
		}
		else if (!strcmp(argv[i], "-p_init_Z")) {
			p_init_Z = atof(argv[++i]);
		}
		else if (!strcmp(argv[i], "-p_meas_Z")) {
			p_meas_Z = atof(argv[++i]);
		}
		else if (!strcmp(argv[i], "-p_iden_cnot")) {
			p_iden_cnot = atof(argv[++i]);
		}
		else if (!strcmp(argv[i], "-p_iden_H")) {
			p_iden_H = atof(argv[++i]);
		}
		else if (!strcmp(argv[i], "-p_iden_init_Z")) {
			p_iden_init_Z = atof(argv[++i]);
		}
		else if (!strcmp(argv[i], "-p_iden_meas_Z")) {
			p_iden_meas_Z = atof(argv[++i]);
		}
      else {
         fprintf(stderr, "Unknown switch: %s\n", argv[i]);
         exit(EXIT_FAILURE);
      }
   }

	p_arr[0] = p_cnot;
	p_arr[1] = p_H;
	p_arr[2] = p_init_Z;
	p_arr[3] = p_meas_Z;
	p_arr[4] = p_iden_cnot;
	p_arr[5] = p_iden_H;
	p_arr[6] = p_iden_init_Z;
	p_arr[7] = p_iden_meas_Z;

	if (p > 0) {
		for (i=0; i<8; i++) p_arr[i] = p;
	}

	// first step is always to load the data base
	db = pe_load_database(fname);


	// EXAMPLE 1: p_arr input

	// calculate formulae from the array of 8 gate error rates
	f = pe_calculate_formulae_from_8ps(p_arr, db);

	// calculate and print desired logical error rate from the formulae
	p_L = pe_calculate_logical_error(d, type, f);
	printf("p_L = %g\n", p_L);

	// free the structure containing the formulae if they are no longer needed
	pe_free_formulae(f);

	// Notes
	// - a different set of data error rates necessitates calculation of a different set of formulae
	// - straightforward to modify the above to generate data for an entire graph rather than a single d, type, p_arr set


	// EXAMPLE 2: p input, same as p_arr input with all entries equal to p
	
	if (p <= 0) p = 0.001;

	// calculate formulae from p alone
	f = pe_calculate_formulae_from_p(p, db);

	// calculate and print desired logical error rate from the formulae
	p_L = pe_calculate_logical_error(d, type, f);
	printf("p_L = %g\n", p_L);

	// free the structure containing the formulae if they are no longer needed
	pe_free_formulae(f);


	// EXAMPLE 3: p.json file input

	// calculate formulae from p.json in p_dir
	f = pe_calculate_formulae_from_p_json(p_dir, db);

	// calculate and print desired logical error rate from the formulae
	p_L = pe_calculate_logical_error(d, type, f);
	printf("p_L = %g\n", p_L);

	// free the structure containing the formulae if they are no longer needed
	pe_free_formulae(f);

	// Notes
	// - an example p.json file is given in the ../data/p directory, it is essentially just an alternative way to specify the
	//   individual gate error rates


	// EXAMPLE 4: gates.json file input, enables the user to specify detailed Pauli error models for each gate

	// calculate formulae from gates.json in gates_dir
	f = pe_calculate_formulae_from_gates_json(gates_dir, db);

	// calculate and print desired logical error rate from the formulae
	p_L = pe_calculate_logical_error(d, type, f);
	printf("p_L = %g\n", p_L);

	// free the structure containing the formulae if they are no longer needed
	pe_free_formulae(f);

	// Notes
	// - the format for gates.json is fairly lengthy, however an example is included in the ../data/gates directory
	// - IMPORTANT: gates.json uses the convention I=0, X=1, Y=2, Z=3,  which is different to the convention used in the next example
	// - global_error sets a baseline error rate
	// - rel_error is set for each gate, and the total probability of error of a gate is global_error * rel_error
	// - num_qubits is the number of qubits touched by each gate
	// - errors is a list of possible errors following the associated gate
	// - pauli_indices marks the beginning of a description of an error, an array of num_qubits errors of the form I=0, X=1, Y=2, Z=3
	// - prob is an integer characterising the relative probability of its associated error, namely prob_i/(Sum prob_i). For example
	//   if all values of prob are equal, the relative probabilities of each error will be equal. If a simulation generates
	//   absolute probabilities of specific errors, as is more likely, one simply uniformly multiplies these absolute probabilities
	//   by a sufficiently large power of 10 to obtain sufficiently large integers to accurately characterise these relative probabilities.


	// EXAMPLE 5: *_em files used as input, enables the user to specify detailed Pauli error models for each gate, p must be specified

	// calculate formulae from *_em files
	f = pe_calculate_formulae_from_ems(p, em_dir, db);

	// calculate and print desired logical error rate from the formulae
	p_L = pe_calculate_logical_error(d, type, f);
	printf("p_L = %g\n", p_L);

	// free the structure containing the formulae if they are no longer needed
	pe_free_formulae(f);

	// Notes
	// - using the definitions in the previous example, the chosen value of p is global_error
	// - the first number in each file is num_qubits
	// - the second number is rel_error
	// - the third number is the number of different possible errors
	// - the following lines contain numbers corresponding to prob (integer as described in EXAMPLE 4) then Pauli error terms
	// - IMPORTANT: *_em files use the convention I=0, X=1, Z=2, Y=3, note the reversal of Z and Y
	// - the final number is the gate duration in arbitrary units (units can be arbitrary provide the same unit is used in all files
	//   and a sufficiently small unit is used such that the duration is an integer


	// final step is to free the memory associated with the database
	pe_free_database(db);

	return EXIT_SUCCESS;
}
