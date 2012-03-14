// update equations for CUDA
#ifndef MY_MAXWELL_CU
#define MY_MAXWELL_CU

#include "my_kernel_help.cu"

__device__ float current_source (float *p, float t) {
//	return sinf( p[1] * t) * expf ( -1 * powf (t * p[2] - p[3], 2.0));
	return (cosf( p[1] * t)*(-2)*(p[2]*t-p[3])*p[2]-p[1]*sinf(p[1]*t)) * expf ( -1 * powf (t * p[2] - p[3], 2.0));
}


__device__ float update_psi (float val, my_grid g, int i, int j, int k, int f) {
	// get the value	
	float psi = af (i,j,k,f);
	int mat_index = am (i,j,k,f);
	float *p = &(g.params[mat_index*4]);
	psi = p[1]*psi + p[2]*val;
	af (i,j,k,f) = psi;
	return psi;
}

__device__ void update_E (	int component,
				float Hx_a[B_XX+1][B_YY+1], 
				float Hx_b[B_XX+1][B_YY+1], 
				float Hy_a[B_XX+1][B_YY+1], 
				float Hy_b[B_XX+1][B_YY+1], 
				float Hz_a[B_XX+1][B_YY+1],
				my_grid g, int i0, int j0, int k0, float t) {
	int i, j;
	i = threadIdx.x+1;
	j = threadIdx.y+1;

	float E, epsilon;
	int mat_type;
	float *p;
	my_field_load (&E, &mat_type, &p, g, i0, j0, k0, component);

	// load epsilon
	epsilon = ae (i0,j0,k0,component);
	

	switch (mat_type) {
		case 0: // hold at zero
			E = 0;
			break;
			
		case 8: // current/soft source
			E += current_source (p, t);
		case 9: // analytic, just update it normally for now
		case 1: // simple_dielectric
			switch (component) {
				case 0: // Ex
					E += epsilon * (Hz_a[i][j] - Hz_a[i][j-1] - Hy_a[i][j] + Hy_b[i][j]);
					break;
				case 1: // Ey
					E += epsilon * (Hx_a[i][j] - Hx_b[i][j] - Hz_a[i][j] + Hz_a[i-1][j]);
					break;
				case 2: // Ez
					E += epsilon * (Hy_a[i][j] - Hy_a[i-1][j] - Hx_a[i][j] + Hx_a[i][j-1]);
					break;
				}
			break;

		case 6: // cpml_dielectric
			switch (component) {
				case 0: // Ex
					E += epsilon * (Hz_a[i][j] - Hz_a[i][j-1] - Hy_a[i][j] + Hy_b[i][j]);
					if (p[2] != 0) 
						E += epsilon * update_psi (Hz_a[i][j]-Hz_a[i][j-1], g, i0, j0+(int)p[2], k0, 0);
					if (p[3] != 0) 
						E -= epsilon * update_psi (Hy_a[i][j]-Hy_b[i][j], g, i0, j0, k0+(int)p[3], 0);
					break;
				case 1:
					E += epsilon * (Hx_a[i][j] - Hx_b[i][j] - Hz_a[i][j] + Hz_a[i-1][j]);
					if (p[1] != 0) 
						E -= epsilon * update_psi (Hz_a[i][j]-Hz_a[i-1][j], g, i0+(int)p[1], j0, k0, 1);
					if (p[3] != 0) 
						E += epsilon * update_psi (Hx_a[i][j]-Hx_b[i][j], g, i0, j0, k0+(int)p[3], 1);
					break;
				case 2:
					E += epsilon * (Hy_a[i][j] - Hy_a[i-1][j] - Hx_a[i][j] + Hx_a[i][j-1]);
					if (p[1] != 0) 
						E += epsilon * update_psi (Hy_a[i][j]-Hy_a[i-1][j], g, i0+(int)p[1], j0, k0, 2);
					if (p[2] != 0) 
						E -= epsilon * update_psi (Hx_a[i][j]-Hx_a[i][j-1], g, i0, j0+(int)p[2], k0, 2);
					break;
				}
			break;
			
		case 7: // cpml_psi (don't do anything)
			return; // this is important, because it skips the write step

		default:
			E = -1;
			break;
		}

	// write out to global memory
	af (i0,j0,k0,component) = E;
	return;
}

__device__ void update_H (	int component,
				float Ex_a[B_XX+1][B_YY+1], 
				float Ex_b[B_XX+1][B_YY+1], 
				float Ey_a[B_XX+1][B_YY+1], 
				float Ey_b[B_XX+1][B_YY+1], 
				float Ez_a[B_XX+1][B_YY+1],
				my_grid g, int i0, int j0, int k0, float t) {
	int i, j;
	i = threadIdx.x;
	j = threadIdx.y;

	float H;
	int mat_type;
	float *p;
	my_field_load (&H, &mat_type, &p, g, i0, j0, k0, component);

	switch (mat_type) {
		case 0: // hold at zero
			H = 0;
			break;
		case 8: // current/soft source
			H += current_source (p, t);
		case 1: // simple_dielectric
			switch (component) {
				case 3:
					H -= p[0] * (Ez_a[i][j+1] - Ez_a[i][j] - Ey_b[i][j] + Ey_a[i][j]);
					break;
				case 4:
					H -= p[0] * (Ex_b[i][j] - Ex_a[i][j] - Ez_a[i+1][j] + Ez_a[i][j]);
					break;
				case 5:
					H -= p[0] * (Ey_a[i+1][j] - Ey_a[i][j] - Ex_a[i][j+1] + Ex_a[i][j]);
					break;
				}
			break;
		case 6: // cpml_dielectric
			switch (component) {
				case 3:
					H -= p[0] * (Ez_a[i][j+1] - Ez_a[i][j] - Ey_b[i][j] + Ey_a[i][j]);
					if (p[2] != 0) 
						H -= p[0] * update_psi (Ez_a[i][j+1]-Ez_a[i][j], g, i0, j0+(int)p[2], k0, 3);
					if (p[3] != 0) 
						H += p[0] * update_psi (Ey_b[i][j]-Ey_a[i][j], g, i0, j0, k0+(int)p[3], 3);
					break;
				case 4:
					H -= p[0] * (Ex_b[i][j] - Ex_a[i][j] - Ez_a[i+1][j] + Ez_a[i][j]);
					if (p[1] != 0) 
						H += p[0] * update_psi (Ez_a[i+1][j]-Ez_a[i][j], g, i0+(int)p[1], j0, k0, 4);
					if (p[3] != 0) 
						H -= p[0] * update_psi (Ex_b[i][j]-Ex_a[i][j], g, i0, j0, k0+(int)p[3], 4);
					break;
				case 5:
					H -= p[0] * (Ey_a[i+1][j] - Ey_a[i][j] - Ex_a[i][j+1] + Ex_a[i][j]);
					if (p[1] != 0) 
						H -= p[0] * update_psi (Ey_a[i+1][j]-Ey_a[i][j], g, i0+(int)p[1], j0, k0, 5);
					if (p[2] != 0) 
						H += p[0] * update_psi (Ex_a[i][j+1]-Ex_a[i][j], g, i0, j0+(int)p[2], k0, 5);
					break;
				}
			break;
		case 7:
			return;
		default:
			H = -1;
			break;
		}
	af (i0,j0,k0,component) = H;
	return;
}

#endif // MY_MAXWELL_CU
