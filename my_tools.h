#ifndef MY_TOOLS_H
#define MY_TOOLS_H

#include <stdlib.h>
#include <stdio.h>
#include <error.h>
#include <string.h>
#include <math.h>

// array access: make accessing a specific value in the array slightly easier
#define aaa(F,i,j,k,yy,zz) F[(i)*yy*zz+(j)*zz+(k)] 

#define MY_NUM_FIELDS 6
#define MY_NUM_MATGRID 6

// very useful constants
static double EPSILON_0 = 8.85e-12; // in farads per meter
static double MU_0 = M_PI*4e-7; // in henrys per meter

// allocate memory
void * cust_alloc (size_t size) {
	register void *value = malloc (size); // allocate a block of memory
	if (value == NULL) // make sure we suceeded in allocating the desired memory
		error (-1, 0, "Virtual memory exhausted.");
	else
		memset (value, 0, size); 
	return value;
	/*This could be replaced with a #define statement:
	 * #define CALLOC(r,s,t) if(((r)=calloc(s,t)) == NULL){error(-1, 0, "Virtual memory exhausted.");} //r is pointer to array, s is the size of the array, t is the byte size of an array element
	 * this might be slightly faster and have less overhead, or totally pointless :)
	 */
}

/*
// reallocate memory
void * cust_realloc (void *ptr, size_t size) {
	register void *value = realloc (ptr, size);
	if (value == 0)
		error (-1, 0, "Virtual memory exhausted.");
	return value;
}

// retrieve parameter from argv
double get_param (char *arg) {
	char *tailptr;
	double parameter;
	parameter = (double) (strtod (arg, &tailptr));
	// check to make sure we got a double
	if ( memcmp (arg, tailptr, strlen(arg)) == 0 )
		error (-1, 0, "Could not interpret '%s' as a double", arg);
	return parameter;
}

// swap two double pointers
void my_swap (double **p1, double **p2) {
	double *temp;
	temp = *p1;
	*p1 = *p2;
	*p2 = temp;
	return;
}
*/

#endif
