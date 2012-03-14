#include <stdio.h>
#include "my_grid.cu"
#include "my_outs.cu"
#include "transfer.cu"
#include <cutil_inline.h>
#include <time.h>
#include <sys/times.h>

#include "my_kernel.cu"
#include "my_analytic.cu"
// #include "memload_kernels.cu"

/* 	This is the main function that runs the simulation. It takes the following parameters:
	paramater 1: -device=X
		X tells us which GPU to run the simulation on. -device=0 will run on the first GPU and -device=N will run on the N+1 GPU. If there are less that N GPUs then it will run on the last GPU.
	parameter 2: input file (e.g. test.h5)
		The input file that contains all the information needed to run the simulation
	parameter 3: output file (e.g. out.h5)
		THrow the outputs here.
*/
int main (int argc, char *argv[]) {
	my_grid g, d_g; // grid structures, d_g is the grid structure on the device
	my_out *outputs; // outputs structure, defines what simulation data we want to return to the user 

	int cntT; // counter variables
	float t; // the time in simulation units
	int percentdone, oldpdone; // so that we don't get too impatient
	oldpdone = 0;

	clock_t start, end; // get the processor time, another way to measure how quickly a simulation runs
	tms tms_temp;
	          
	printf ("\nInitialize grid\n---\n"); // read in data from the input file
	grid_initialize (&g, argv[2]); 
	int oo = outs_initialize (&outputs, argv[2], g); // oo tells us how many outputs the user has requested

    	// CUT_DEVICE_INIT(2, argv); // initialize the device

	printf ("\nTransfer data to device\n---\n"); // copy over data to device
	d_g = grid2device (g);
		
	// help define how the simulation will run on the GPU 
	dim3 threads (B_XX, B_YY);
	dim3 grid (g.xx / (B_XX), g.yy /(B_YY));
	printf ("\nThread block dimensions: [%d, %d]\n", B_XX, B_YY);
	printf ("Grid dimensions: [%d, %d]\n", g.xx/B_XX, g.yy/B_YY);

	printf("\nBegin testing an algorithm\n---\n          ");
	start = times(&tms_temp);
	for (cntT=0 ; cntT<g.tt;  cntT++) {
		t = (float)(cntT);
		kernel_E <<<grid, threads>>> (d_g, t); // update E-fields
		kernel_H <<<grid, threads>>> (d_g, t); // update H-fields
		outs_extract (outputs, oo, g, d_g, cntT, grid, threads); // extract data, if needed, to output structures. Analytical data is also calculated and extracted in this function

		percentdone = (int)(100.00*cntT/g.tt);
		if (percentdone > oldpdone) {
			printf("\b\b\b\b\b\b\b\b %2d done",percentdone);
			fflush(stdout);
			oldpdone = percentdone;
		}
	}
	end = times(&tms_temp);
	printf ("\n\n%d time units.\n", end-start);

	printf ("\nStore result in hdf  file\n---\n");
	outs_write (outputs, oo, argv[3], g);

        printf ("Freeing memory\n");
	// CUDA_SAFE_CALL(cudaFree(	));
	// write something here!
	printf("Woo hoo, all done\n");
	// CUT_EXIT(argc, argv);	
	return 0;
}

