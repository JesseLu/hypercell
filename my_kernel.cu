// kernel update for CUDA
#ifndef MY_KERNEL_CU 
#define MY_KERNEL_CU 

#include "my_grid.cu"
#include "my_maxwell.cu"
#include "my_kernel_help.cu"

__global__ void kernel_E (my_grid g, float t) {
	int i, j, k;

	// initialize index values
	i = tx + bx*B_XX;
	j = ty + by*B_YY;
	
	__shared__ float Hx_a[B_XX+1][B_YY+1];
	__shared__ float Hy_a[B_XX+1][B_YY+1];
	__shared__ float Hz_a[B_XX+1][B_YY+1];
	__shared__ float Hx_b[B_XX+1][B_YY+1];
	__shared__ float Hy_b[B_XX+1][B_YY+1];

	for (k=0 ; k<g.zz; k++) {
		// load the bottom layer
		if (k==0) { // just starting out, so we need to load
			load_E (g, Hx_b, i, j, 0, 3);
			load_E (g, Hy_b, i, j, 0, 4);
			// now make them negative, for the TE symmetry case
			if (g.symz == 1) {
				make_negative (Hx_b);
				make_negative (Hy_b);
			}
			//note for TM case, not making it negative gives Ex,Ey =0  on the mirror plane	
		}
		else
			{ // obtain the bottom layer through a swap operation
			transfer_E (Hx_a, Hx_b);
			transfer_E (Hy_a, Hy_b);
			}

		// load the current layer
		load_E (g, Hx_a, i, j, k, 3);
		load_E (g, Hy_a, i, j, k, 4);
		load_E (g, Hz_a, i, j, k, 5);


		__syncthreads(); // wait for all values to be loaded in
		// update Ex
		update_E (0, Hx_a, Hx_b, Hy_a, Hy_b, Hz_a, g, i, j, k, t); // Ex
		update_E (1, Hx_a, Hx_b, Hy_a, Hy_b, Hz_a, g, i, j, k, t); // Ey
		update_E (2, Hx_a, Hx_b, Hy_a, Hy_b, Hz_a, g, i, j, k, t); // Ez

		__syncthreads(); // wait for all computation to be finished before loading in new values
	}

		
	return;
}



__global__ void kernel_H (my_grid g, float t) {
	int i, j, k;

	// initialize index values
	i = tx + bx*B_XX;
	j = ty + by*B_YY;
	
	__shared__ float Ex_a[B_XX+1][B_YY+1];
	__shared__ float Ey_a[B_XX+1][B_YY+1];
	__shared__ float Ez_a[B_XX+1][B_YY+1];
	__shared__ float Ex_b[B_XX+1][B_YY+1];
	__shared__ float Ey_b[B_XX+1][B_YY+1];

	for (k=0 ; k<g.zz-1; k++) {
		// load the bottom layer
		if (k==0) { 
			// load the current layer
			load_H (g, Ex_a, i, j, 0, 0);
			load_H (g, Ey_a, i, j, 0, 1);
			//if (g.symz == -1) {
			//	load_H (g, Ex_a, i,j,1,0);
			//	load_H (g, Ey_a, i,j,1,1);
			//	make_negative(Ex_a);
			//	make_negative(Ey_a);
			//}
		}
			
		else
			{ // obtain the the new current layer through a swap operation
			transfer_H (Ex_b, Ex_a);
			transfer_H (Ey_b, Ey_a);
			}

		// load the upper layer
		load_H (g, Ex_b, i, j, k+1, 0);
		load_H (g, Ey_b, i, j, k+1, 1);
		// and Ez for the current layer
		load_H (g, Ez_a, i, j, k, 2);


		__syncthreads(); // wait for all values to be loaded in

		update_H (3, Ex_a, Ex_b, Ey_a, Ey_b, Ez_a, g, i, j, k, t);
		update_H (4, Ex_a, Ex_b, Ey_a, Ey_b, Ez_a, g, i, j, k, t);
		update_H (5, Ex_a, Ex_b, Ey_a, Ey_b, Ez_a, g, i, j, k, t);

		__syncthreads(); // wait for all computation to be finished before loading in new values
	}

		
	return;
}


#endif // MY_KERNEL_CU
