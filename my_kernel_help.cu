#ifndef MY_KERNEL_HELP
#define MY_KERNEL_HELP

#define B_XX 16
#define B_YY 16

#define tx threadIdx.x
#define ty threadIdx.y

#define bx blockIdx.x
#define by blockIdx.y

#define af(i,j,k,f) g.field[f][(i) + (j)*g.xx+ (k)*g.xx*g.yy] 
#define ae(i,j,k,f) g.epsilon[f][(i) + (j)*g.xx+ (k)*g.xx*g.yy] 
#define am(i,j,k,f) g.mat[f][(i) + (j)*g.xx+ (k)*g.xx*g.yy] 

__device__ void my_field_load (float *A, int *mat_type, float **p, my_grid g, int i, int j, int k, int f) {
	*A = af (i,j,k,f);
	int mat_index = am (i,j,k,f);
	*mat_type = g.type_index[mat_index]; 
	*p = &(g.params[mat_index*4]);
	return;
}

__device__ void load_E (my_grid g, float H[17][17], int i, int j, int k, int f) {
	H[tx+1][ty+1] = af (i,j,k,f);
	// H[tx+1][ty+1] = 2;
	if (tx == 0) {
		if (bx == 0)
			H[tx][ty+1] = -1;
		else
			H[tx][ty+1] = af (i-1,j,k,f);
	// H[tx+1][ty+1] = 2;
			// H[tx][ty+1] = tx+ty;
		}
	if (ty == 0) {
		if (by == 0)
			H[tx+1][ty] = -1;
		else
			H[tx+1][ty] = af (i,j-1,k,f);
	// H[tx+1][ty+1] = 2;
		}
	return;
}

__device__ void load_H (my_grid g, float E[17][17], int i, int j, int k, int f) {
	E[tx][ty] = af (i,j,k,f);
	// E[tx][ty] = i+j;
	if (tx == B_XX-1) {
		if (bx == (g.xx/B_XX)-1)
			E[tx+1][ty] = -1;
		else
			E[tx+1][ty] = af (i+1,j,k,f);
			// E[tx+1][ty] = -1;
			// H[tx][ty+1] = tx+ty;
		}
	if (ty == B_YY-1) {
		if (by == (g.yy/B_YY)-1)
			E[tx][ty+1] = -1;
		else
			E[tx][ty+1] = af (i,j+1,k,f);
		}
	return;
}

// transfers from A to B
__device__ void transfer_E (float A[17][17], float B[17][17]) {
	B[tx+1][ty+1] = A[tx+1][ty+1];
	if (tx == 0)
		B[tx][ty+1] = A[tx][ty+1];
	if (ty == 0)
		B[tx+1][ty] = A[tx+1][ty];
	return;
}

__device__ void transfer_H (float A[17][17], float B[17][17]) {
	B[tx][ty] = A[tx][ty];
	if (tx == B_XX-1)
		B[tx+1][ty] = A[tx+1][ty];
	if (ty == B_YY-1)
		B[tx][ty+1] = A[tx][ty+1];
	return;
}
	
__device__ void make_negative(float A[17][17]) {
	A[tx+1][ty+1] = A[tx+1][ty+1] * -1;
	if (tx == 0)
		A[tx][ty+1] = A[tx][ty+1] * -1;
	if (ty == 0)
		A[tx+1][ty] = A[tx+1][ty] * -1;
	return;
}

#endif
