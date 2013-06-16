#ifndef RANDOM_H
#define RANDOM_H

#define MERS_N	624
#define MERS_M	397
#define MERS_R	31
#define MERS_U	11
#define MERS_S	7
#define MERS_T	15
#define MERS_L	18
#define MERS_A	0x9908B0DF
#define MERS_B	0x9D2C5680
#define MERS_C	0xEFC60000

/**
 * \brief A random number generator
 */
typedef struct rng RNG;
struct rng {	
    /*@{*/
    /** The current internal state */
    unsigned int *mt;
    /** The saved internal state, if any */
    unsigned int *mt2;
    /*@}*/
    
    /*@{*/
    /** The current state index */
	int mti;
    /** The saved state index, if any */
    int mti2;
    /*@}*/
};

RNG *rng_create();
void rng_free(RNG *rng);
void rng_init(RNG *rng, int seed); 
void rng_init_by_array(RNG *rng, int const seeds[], int num_seeds); 
void rng_randoms(RNG *rng, int *s0, int *s1); 
unsigned int rng_randomui(RNG *rng);
double rng_randomd(RNG *rng);
void rng_save(RNG *rng);
void rng_restore(RNG *rng);

//
// Cautionary functionality
//
// Use the above functions instead, the below exists
// to maintain backward compatibility with previous code.
// Only use in the case it is otherwise infeasible to pass
// the RNG around.
//

void Init0(int seed);
void RandomInitByArray(int const seeds[], int NumSeeds);
void randoms(int *s0, int *s1);
unsigned int randomui();
double randomd();
void save_rand();
void restore_rand();
void save_rand_to_file(const char *filename); 
void restore_rand_from_file(const char *filename);

#endif
