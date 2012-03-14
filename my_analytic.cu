// use this to compute stuff about the grid
#ifndef MY_ANALYTIC_CU
#define MY_ANALYTIC_CU

#include "my_grid.cu"
#include "my_analytic_kernel.cu"
// #include <cutil.h>

// this is the big function which will calculate all the analytical stuff we want (or possible)
float analytic_calculate (my_grid g, dim3 grid_size, dim3 block_size, int analysis_type) {
	// allocate space in global memory to hold the result for each thread block
	// just do oen analytic calculation for now
	float *h_answers, *d_answers;
	int memsize = sizeof (float) * grid_size.x * grid_size.y;
	(cudaMalloc((void**) &d_answers, memsize)); // allocate memory in global (device) memory
	h_answers = (float*) cust_alloc (memsize); // allocate memory in host memory (RAM)

	// calculate
	analytic_energy <<<grid_size, block_size>>> (g, d_answers, grid_size.x, analysis_type);

	// transfer stuff back to host
	(cudaMemcpy(h_answers, d_answers, memsize, cudaMemcpyDeviceToHost) );

	// sum up everything to get just one answer for each analytic
	int cnt;
	float sum = 0;
	for (cnt=0 ; cnt < grid_size.x * grid_size.y ; cnt++) 
		sum = sum + h_answers [cnt];
	// printf ("%e\n", sum);
	return sum;
}

#endif
