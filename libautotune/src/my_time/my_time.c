#include <stdlib.h>
#include <sys/time.h>
#include "my_time.h"

/**
 * \brief Gets the time in seconds
 *
 * \return The time in seconds.
 */
double double_time() {
	struct timeval t1970;

	gettimeofday(&t1970, NULL);

	return (double)t1970.tv_sec + (double)t1970.tv_usec/1000000;
}
