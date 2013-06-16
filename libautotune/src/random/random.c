#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <limits.h>
#include "random.h"
#include "../memory/memory.h"

/**
 * \file 
 * Example of use
 *
 * \code
 * RNG *rng;
 * int s0, s1;
 *
 * unsigned int rand_ui;
 * double rand_d;
 *
 * // Zeros will seed using current time
 * s0 = s1 = 0;
 *
 * // Create and initialise RNG
 * rng = rng_create();
 * rng_randoms(rng, &s0, &s1);
 *
 * // Get some random numbers
 * rand_ui = rng_randomui(rng); // Returns an unsigned integer
 * rand_d = rng_randomd(rng); // Returns a double
 * 
 * // Done, now we can free the RNG
 * rng_free(rng);
 * \endcode
 */

/**
 * \brief Creates a new random number generator
 *
 * Creates a new random number generator, implementing the Mersenne Twister MT19937.
 *
 * \returns A random number generator
 */
RNG *rng_create() {
	RNG *rng;

	rng = (RNG *)my_malloc(sizeof(RNG));
	
	rng->mt = (unsigned int *)my_malloc(MERS_N * sizeof(unsigned int));
	rng->mt2 = (unsigned int *)my_malloc(MERS_N * sizeof(unsigned int));
	rng->mti = 0;
	rng->mti2 = 0;

	return rng;
}

/**
 * \brief Frees a random number generator
 *
 * \param[in] rng The random number generator to be freed
 */
void rng_free(RNG *rng) {
	free(rng->mt);
	free(rng->mt2);
	free(rng);
}

/**
 * \brief [Private] Initialises the Mersenne Twister initial state
 *
 * Intialises the random number generator state. This function should not need
 * to be called outside of this library, see rng_randoms.
 *
 * \param[in,out] rng The random number generator to be initialised
 * \param[in] seed The seed to initialise the randome number generator
 *
 * \sa rng_randoms
 */
void rng_init(RNG *rng, int seed) {
	const unsigned int factor = 1812433253UL;
	rng->mt[0]= seed;
	for (rng->mti = 1; rng->mti < MERS_N; rng->mti++) {
		rng->mt[rng->mti] = (factor * (rng->mt[rng->mti-1] ^ (rng->mt[rng->mti-1] >> 30)) + rng->mti);
	}
}

/**
 * \brief [Private] Performs the Mersenne Twisting operation to obtain the
 * initial state 
 *
 * Initialises the random number generator state. This function should not need
 * to be called outside of this library, see rng_randoms.
 *
 * \param[in,out] rng The random number generator to be initialised
 * \param[in] seeds An array of seeds to be used
 * \param[in] num_seeds The number of seeds that have been passed in
 *
 * \sa rng_randoms
 */
void rng_init_by_array(RNG *rng, int const seeds[], int num_seeds) {
	// Seed by more than 32 bits
	int i, j, k;

	// Initialize
	rng_init(rng, 19650218);

	if (num_seeds <= 0) return;

	// Randomize mt[] using whole seeds[] array
	i = 1;  j = 0;
	k = (MERS_N > num_seeds ? MERS_N : num_seeds);
	for (; k; k--) {
		rng->mt[i] = (rng->mt[i] ^ ((rng->mt[i-1] ^ (rng->mt[i-1] >> 30)) * 1664525UL)) + (unsigned int)seeds[j] + j;
		i++; j++;
		if (i >= MERS_N) {
			rng->mt[0] = rng->mt[MERS_N-1]; 
			i = 1;
		}
		if (j >= num_seeds) { 
			j = 0;
		}
	}
	for (k = MERS_N-1; k; k--) {
		rng->mt[i] = (rng->mt[i] ^ ((rng->mt[i-1] ^ (rng->mt[i-1] >> 30)) * 1566083941UL)) - i;
		if (++i >= MERS_N) {
			rng->mt[0] = rng->mt[MERS_N-1]; 
			i = 1;
		}
	}
	rng->mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

	// Randomize some more
	rng->mti = 0;
	for (i = 0; i <= MERS_N; i++) {
		rng_randomui(rng);
	}
}

/**
 * \brief Initialises the random number generator by invoking the Mersenne
 * Twister initalisation. 
 *
 * Initialises the random number generator by invoking the Mersenne
 * Twister initalisation. Passing in no seeds will generate new seeds using the
 * current time.
 *
 * \param[in,out] rng The random number generator to be initialised
 * \param[in] s0 The first seed
 * \param[in] s1 The second seed 
 */
void rng_randoms(RNG *rng, int *s0, int *s1) {
	struct timeval t1970;
	int seeds[2];

	if (*s0 == 0 && *s1 == 0) {
		// microseconds since 1/1/1970 used to initialise random number generator
		gettimeofday(&t1970, NULL);
		*s0 = seeds[0] = t1970.tv_sec;
		*s1 = seeds[1] = t1970.tv_usec;
	}
	else {
		seeds[0] = *s0;
		seeds[1] = *s1;
	}

	// seed Mersenne Twister random number generator
	rng_init_by_array(rng, seeds, 2);
}

/** 
 * \brief Generates a random unsigned integer using the Mersenne Twister
 * algorithm
 *
 * \param[in] rng The random number generator used to generate the number. This
 * must first be initialised using rng_randoms.
 *
 * \sa rng_randoms
 *
 * \return A pseudo-random unsigned integer
 */
unsigned int rng_randomui(RNG *rng) {
	// Generate 32 random bits
	unsigned int y;

	if (rng->mti >= MERS_N) {
		// Generate MERS_N words at one time
		const unsigned int LOWER_MASK = (1LU << MERS_R) - 1;  // Lower MERS_R bits
		const unsigned int UPPER_MASK = 0xFFFFFFFF << MERS_R; // Upper (32 - MERS_R) bits
		static const unsigned int mag01[2] = {0, MERS_A};

		int kk;
		for (kk=0; kk < MERS_N-MERS_M; kk++) {	 
			y = (rng->mt[kk] & UPPER_MASK) | (rng->mt[kk+1] & LOWER_MASK);
			rng->mt[kk] = rng->mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

		for (; kk < MERS_N-1; kk++) {	 
			y = (rng->mt[kk] & UPPER_MASK) | (rng->mt[kk+1] & LOWER_MASK);
			rng->mt[kk] = rng->mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}		

		y = (rng->mt[MERS_N-1] & UPPER_MASK) | (rng->mt[0] & LOWER_MASK);
		rng->mt[MERS_N-1] = rng->mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
		rng->mti = 0;
	}
	y = rng->mt[rng->mti++];

	// Tempering (May be omitted):
	y ^=  y >> MERS_U;
	y ^= (y << MERS_S) & MERS_B;
	y ^= (y << MERS_T) & MERS_C;
	y ^=  y >> MERS_L;

	return y;
}

/**
 * \brief Generates a random double using the Mersenne Twister algorithm
 * between 0 and 1
 *
 * Generates a random double using the Mersenne Twister algorithm. An alias for
 * generating an unsigned integer and casting it to a double and dividing by
 * UINT_MAX.
 *
 * \param[in] rng The random number generator used to generate the number. This
 * must first be initialised using rng_randoms.
 *
 * \sa rng_randoms
 *
 * \return A pseudo-random double between 0 and 1
 */
double rng_randomd(RNG *rng) {
	return (double)rng_randomui(rng)/UINT_MAX;
}

/**
 * \brief Saves the current internal state of a random number generator
 *
 * Saves the current internal state of a random number generator. Use this to
 * save an intermediate state. To restore a random number generator from
 * scratch, initialise it with the same initial seed.
 *
 * \param[in] rng The random number generator used to generate the number. This
 * must first be initialised using rng_randoms.
 *
 * \sa rng_randoms
 */
void rng_save(RNG *rng) {
	int i;

	for (i = 0; i < MERS_N; i++) rng->mt2[i] = rng->mt[i];
	rng->mti2 = rng->mti;
}

/**
 * \brief Restores a saved internal state of a random number generator
 *
 * Restores a saved internal state of a random number generator. The state must
 * first have been saved using rng_save.
 *
 * \param[in] rng The random number generator used to generate the number. This
 * must first be initialised using rng_randoms.
 *
 * \sa rng_save
 */
void rng_restore(RNG *rng) {
	int i;

	for (i = 0; i < MERS_N; i++) rng->mt[i] = rng->mt2[i];
	rng->mti = rng->mti2;
}





//
// Cautionary functionality
//
// Use the above functions instead, the below exists
// to maintain backward compatibility with previous code.
// Only use in the case it is otherwise infeasible to pass
// the RNG around.
//

// State vector and backup
static unsigned int mt[MERS_N];	
static int mti;					

// Backup for state vector and index
static unsigned int mt2[MERS_N];  
static int mti2;						

void Init0(int seed) {
	// Seed generator
	const unsigned int factor = 1812433253UL;
	mt[0]= seed;
	for (mti=1; mti < MERS_N; mti++) {
		mt[mti] = (factor * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
	}
}

void RandomInitByArray(int const seeds[], int NumSeeds) {
	// Seed by more than 32 bits
	int i, j, k;

	// Initialize
	Init0(19650218);

	if (NumSeeds <= 0) return;

	// Randomize mt[] using whole seeds[] array
	i = 1;  j = 0;
	k = (MERS_N > NumSeeds ? MERS_N : NumSeeds);
	for (; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + (unsigned int)seeds[j] + j;
		i++; j++;
		if (i >= MERS_N) {
			mt[0] = mt[MERS_N-1]; 
			i = 1;
		}
		if (j >= NumSeeds) { 
			j = 0;
		}
	}
	for (k = MERS_N-1; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
		if (++i >= MERS_N) {
			mt[0] = mt[MERS_N-1]; 
			i = 1;
		}
	}
	mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

	// Randomize some more
	mti = 0;
	for (i = 0; i <= MERS_N; i++) randomui();
}

void randoms(int *s0, int *s1) {
	struct timeval t1970;
	int seeds[2];

	if (*s0 == 0 && *s1 == 0) {
		// microseconds since 1/1/1970 used to initialise random number generator
		gettimeofday(&t1970, NULL);
		*s0 = seeds[0] = t1970.tv_sec;
		*s1 = seeds[1] = t1970.tv_usec;
	}
	else {
		seeds[0] = *s0;
		seeds[1] = *s1;
	}

	// seed Mersenne Twister random number generator
	RandomInitByArray(seeds, 2);
}

unsigned int randomui() {
	// Generate 32 random bits
	unsigned int y;

	if (mti >= MERS_N) {
		// Generate MERS_N words at one time
		const unsigned int LOWER_MASK = (1LU << MERS_R) - 1;  // Lower MERS_R bits
		const unsigned int UPPER_MASK = 0xFFFFFFFF << MERS_R; // Upper (32 - MERS_R) bits
		static const unsigned int mag01[2] = {0, MERS_A};

		int kk;
		for (kk=0; kk < MERS_N-MERS_M; kk++) {	 
			y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

		for (; kk < MERS_N-1; kk++) {	 
			y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}		

		y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
		mti = 0;
	}
	y = mt[mti++];

	// Tempering (May be omitted):
	y ^=  y >> MERS_U;
	y ^= (y << MERS_S) & MERS_B;
	y ^= (y << MERS_T) & MERS_C;
	y ^=  y >> MERS_L;

	return y;
}

double randomd() {
	return (double)randomui()/UINT_MAX;
}

void save_rand() {
	int i;

	for (i=0; i<MERS_N; i++) mt2[i] = mt[i];
	mti2 = mti;
}

void restore_rand() {
	int i;

	for (i=0; i<MERS_N; i++) mt[i] = mt2[i];
	mti = mti2;
}

void save_rand_to_file(const char *filename) {
	FILE *fp;
	char line[256];
	int i;

	fp = fopen(filename, "w+");
	if (fp == NULL) {
		printf("Could not open file: %s\n", filename);
		exit(1);
	}
	for (i=0; i<MERS_N; i++) {
		sprintf(line, "%u\n", mt[i]); 
		fputs(line, fp);
	}
	sprintf(line, "%d\n", mti); 
	fputs(line, fp);
	
	fclose(fp);
}

void restore_rand_from_file(const char *filename) {
	FILE *fp;
	int i;
	unsigned int ui;

	fp = fopen(filename, "r+");
	if (fp == NULL) {
		printf("Could not open file: %s\n", filename);
		exit(1);
	}
	for (i=0; i<MERS_N; i++) {
		fscanf(fp, "%u\n", &ui);
		mt[i] = ui;
	}
	fscanf(fp, "%d\n", &i);
	mti = i;

	fclose(fp);
}
