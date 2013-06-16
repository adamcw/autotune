#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "memory.h"
#if defined(__APPLE__) && defined(__MACH__)
	#include <mach/mach.h>
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	#include <unistd.h>
	#include <sys/resource.h>
#endif

/**
 * \brief Get current program memory usage in kilobytes.
 *
 * \return A size_t representing the current program memory usage in kilobytes.
 */
int my_mem() {
#if defined(__APPLE__) && defined(__MACH__)
	struct mach_task_basic_info info;
	mach_msg_type_number_t infoCount;
	
	infoCount = MACH_TASK_BASIC_INFO_COUNT;
	if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount ) != KERN_SUCCESS) {
		return (int)0;
	}
	return (int)(info.resident_size / 1024L);
#else 
	long rss;;
	FILE* fp;
	
	rss = 0L;
	fp = NULL;

	if ((fp = fopen( "/proc/self/statm", "r" )) == NULL) {
		return (int)0L;
	}

	if (fscanf(fp, "%*s%ld", &rss) != 1) {
		fclose(fp);
		return (int)0L;
	}

	fclose(fp);

	return (int)((size_t)rss * (size_t)sysconf(_SC_PAGESIZE) / 1024L);
#endif
}

/**
 * \brief A safe implementation of malloc
 *
 * A safe implementation of malloc. Will exit if no memory available 
 *
 * \param[in] size The size of memory to be allocated.
 *
 * \return A pointer to the allocated memory block.
 */
void *my_malloc(size_t size) {
	void *ptr;

	if (size > 0) {
		ptr = malloc(size);
		if (ptr == NULL) {
			exit(1);
		}
		return ptr;
	}

	return NULL;
}

/**
 * \brief A safe implementation of realloc
 *
 * A safe implementation of realloc. Will exit if no memory available 
 *
 * \param[in] ptr The pointer to the memory block to be reallocated
 * \param[in] size the size the memory block should be resized to
 *
 * \return A pointer to the new reallocated memory block.
 */
void *my_realloc(void *ptr, size_t size) {
	void *ptr2;

	ptr2 = realloc(ptr, size);
	if (ptr2 == NULL && size > 0) {
		exit(1);
	}
	return ptr2;
}

/**
 * \brief A safe implementation of calloc
 *
 * A safe implementation of calloc. Will exit if no memory available 
 *
 * \param[in] n The number of elements to allocate
 * \param[in] size The size of each element
 *
 * \return A pointer to the new calloced memory block.
 */
void *my_calloc(int n, int size)
{
	void *ptr;

	if (n > 0 && size > 0) {
		ptr = calloc(n, size);
		if (ptr == NULL) {
			exit(1);
		}
		return ptr;
	}

	return NULL;
}

/**
 * \brief Allocates and initialises a 2D array of memory
 *
 * \param[in] n1 The number of elements in the ith dimension
 * \param[in] n2 The number of elements in the jth dimension
 * \param[in] size The size of each element in the array
 *
 * \return A pointer to the newly allocated memory block
 */
void *my_2d_calloc(int n1, int n2, int size)		  
{
	int i;
	void **pptr;

	pptr = (void **)my_calloc(n1, sizeof(void *));
		
	for (i = 0; i < n1; i++) {
		pptr[i] = my_calloc(n2, size);
	}

	return pptr;
}

/**
 * \brief Allocates and initialises a 3D array of memory
 *
 * \param[in] n1 The number of elements in the ith dimension
 * \param[in] n2 The number of elements in the jth dimension
 * \param[in] n3 The number of elements in the kth dimension
 * \param[in] size The size of each element in the array
 *
 * \return A pointer to the newly allocated memory block
 */
void *my_3d_calloc(int n1, int n2, int n3, int size)
{
	int i;
	void ***ppptr;

	ppptr = (void ***)my_calloc(n1, sizeof(void **));
		
	for (i = 0; i < n1; i++) {
		ppptr[i] = (void **)my_2d_calloc(n2, n3, size);
	}

	return ppptr;
}

/**
 * \brief Frees a 2D array of memory
 *
 * \param[in] n1 The size of the ith dimension
 * \param[in] pptr The block of memory to be freed
 */
void my_2d_free(int n1, void **pptr)
{
	int i;

	for (i = 0; i < n1; i++) {
		free(pptr[i]);
	}
	  
	free(pptr);
}

/**
 * \brief Frees a 2D array of memory
 *
 * \param[in] n1 The size of the ith dimension 
 * \param[in] n2 The size of the jth dimension
 * \param[in] ppptr The block of memory to be freed
 */
void my_3d_free(int n1, int n2, void ***ppptr)
{
	int i;

	for (i = 0; i < n1; i++) {
		my_2d_free(n2, ppptr[i]);
	}

	free(ppptr);
}
