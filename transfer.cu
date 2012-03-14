#ifndef TRANSFER_CU
#define TRANSFER_CU

#include <stdio.h>
#include "my_grid.cu"
// #include <cutil.h>

__device__ float *d_field[6];

// copy an array into device memory
void* my_copy2device (void *h_data, unsigned int memsize) {
	// allocate device memory
	void *d_data;
	(cudaMalloc((void**) &d_data, memsize));

	// copy node structure to device
	(cudaMemcpy(d_data, h_data, memsize, cudaMemcpyHostToDevice) );
	return d_data;
}


// Creates a copy of the grid structure which can be used by the kernel functions which run on the device. 
// Basically a structure whose pointers point to arrays in device global memory, instead of host memory.
my_grid grid2device (my_grid g) {
	int cnt;
	for (cnt=0 ; cnt<6 ; cnt++) {
		g.field[cnt] = (float*) my_copy2device (g.field[cnt], sizeof (float) * (g.xx*g.yy*g.zz)); 
		g.mat[cnt] = (int*) my_copy2device (g.mat[cnt], sizeof (int) * (g.xx*g.yy*g.zz)); 
	}
	for (cnt=0 ; cnt<3 ; cnt++) {
		g.epsilon[cnt] = (float*) my_copy2device (g.epsilon[cnt], sizeof (float) * (g.xx*g.yy*g.zz)); 
	}
	g.type_index = (int*) my_copy2device (g.type_index, sizeof (int) * g.mm);
	g.params = (float*) my_copy2device (g.params, sizeof (float) * g.mm * MY_MAXMATPARAMS);
	

	return g;
}

#endif // TRANSFER_CU
