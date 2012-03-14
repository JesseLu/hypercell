// use this to compute stuff about the grid
#ifndef MY_ANALYTIC_KERNEL_CU
#define MY_ANALYTIC_KERNEL_CU

#include "my_kernel_help.cu"

// load what we need to do the analytic stuff
__device__ void my_analytic_load (float A[6], int mat_type[6], float *p[3], float epsilon[3], my_grid g, int i, int j, int k) {
	int cntF, mat_index, cnt;
	for (cntF=0 ; cntF<3 ; cntF++) {
		A[cntF] = af (i,j,k,cntF);
		mat_index = am (i,j,k,cntF);
		mat_type[cntF] = g.type_index[mat_index]; 
		// if (cntF < 3) {
			p[cntF] = &(g.params[mat_index*4]);
			epsilon[cntF] = ae (i,j,k,cntF);
		// }
	}
	for (cntF=3 ; cntF<6 ; cntF++) {
		A[cntF] = af (i,j,k,cntF);
		mat_index = am (i,j,k,cntF);
		mat_type[cntF] = g.type_index[mat_index]; 
		// if (cntF < 3) {
			// p[cntF] = &(g.params[mat_index*4]);
			// epsilon[cntF] = ae (i,j,k,cntF);
		// }
	}
	return;
}
	
__device__ void my_analytic_kernel (float A[6], int mat_type[6], float *p[3], float epsilon[3], int analysis_type, float dt, float running_sum[B_XX][B_YY]) {
	int cnt;
	int temp;
	if (analysis_type == -1) { // global energy in dielectrics
		for (cnt=0 ; cnt<3 ; cnt++) 
			running_sum[tx][ty] += ((mat_type[0]!=6)&&(mat_type[0]!=7)&&(mat_type[1]!=6)&&(mat_type[1]!=7)&&(mat_type[2]!=6)&&(mat_type[2]!=7)) * dt * powf (A[cnt], 2.0) / epsilon[cnt] * 0.5;
		for (cnt=3 ; cnt<6 ; cnt++) 
			running_sum[tx][ty] += ((mat_type[0]!=6)&&(mat_type[0]!=7)&&(mat_type[1]!=6)&&(mat_type[1]!=7)&&(mat_type[2]!=6)&&(mat_type[2]!=7)) * powf (A[cnt], 2.0) * 0.5;
	}
	else if (analysis_type >= 0) {
		for (cnt=0 ; cnt<3 ; cnt++) {
			// temp = (int)p[cnt][0];
			if ((mat_type[cnt]==9) && (p[cnt][2] == analysis_type))
				running_sum[tx][ty] += p[cnt][1] * A[cnt] *  A[(int)p[cnt][0]] ;
		}
			// running_sum[tx][ty] += 1;
			// running_sum[tx][ty] += (mat_type[cnt]==9) * ((int)p[cnt][2]==analysis_type) * (p[cnt][0] + 1.0);
			// running_sum[tx][ty] += (mat_type[cnt]==9) * ((int)p[cnt][2]==analysis_type) * p[cnt][1] * A[cnt] * A[(int)p[cnt][0]];
	}
	else {
		// not good if we actually get here
	}


	return;
}

__global__ void analytic_energy (my_grid g, float *answers, int A_XX, int analysis_type) {

	int i, j, k; // tell us where we are on the main grid
	int cnti, cntj;
	__shared__ float running_sum[B_XX][B_YY]; // this sums for each thread
	float big_sum; // the final answer
	
	float A[6]; // stores the 6 field components that we need
	int mat_type[6]; // material types only for E-fields
	float *p[3]; // parameters for those matrial types
	float epsilon[3]; // epsilon for E values

	// initialize index values
	i = tx + bx*B_XX;
	j = ty + by*B_YY;

	running_sum[tx][ty] = 0;
	for (k=0 ; k < g.zz ; k++) {
		my_analytic_load (A, mat_type, p, epsilon, g, i, j, k);
		my_analytic_kernel (A, mat_type, p, epsilon, analysis_type, g.dt, running_sum);
		// running_sum[tx][ty] = running_sum[tx][ty] + 1;
	}

	// don't forget to syncthreads to make sure everyone's done
	__syncthreads();

	// sum up stuff
	if (tx==0 && ty==0) {
		big_sum = 0;
		for (cnti=0 ; cnti<B_XX ; cnti++)
			for (cntj=0 ; cntj<B_YY ; cntj++)
				big_sum = big_sum + running_sum[cnti][cntj];
		answers[bx + by*A_XX] = big_sum;
	}
	else
		return;
}


#endif
