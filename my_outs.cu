#ifndef MY_OUTS_H
#define MY_OUTS_H

#include <stdlib.h>
#include <H5LT.h>
#include "my_tools.h"
#include "my_analytic.cu"
// #include <cutil.h>

// defines what kind of outputs we want for our simulation
typedef struct my_out {
	int dims; // the dimensions of the out command
	int pos[3]; // position for the out command 
	int *tshots; // the t array
	int tt; // length of tshots
	float *data[6]; // holds the data to be outputted (before we write it to disk)
	int counter; // tells us where we are in the tshots list
	int tshot_size; // how big each time shot is
} my_out;



// load everything in first, then organize it
int outs_initialize (my_out **o, char *hdf5_filename, my_grid g) {
	hid_t file_id; // file id for hdf file
	int cnt, cntA;
	int oo[2];

	printf ("Opening hdf5 file %s...", hdf5_filename);
	file_id = H5Fopen (hdf5_filename, H5F_ACC_RDONLY, H5P_DEFAULT); // open the hdf5 file
	if ( file_id <= 0 )
		error (-1, 0, "error.\n");	
	printf ("file id: %d\n", file_id);

	// get oo
	printf ("Reading in variables...\n");
	H5LTread_dataset (file_id, "/info/oo", H5T_NATIVE_INT, oo); // tells us how big the other arrays are
	printf ("(oo) = [%d, %d]\n", oo[0], oo[1]);
	
	*o = (my_out*) cust_alloc (sizeof (my_out)	* oo[0]);
	
	my_out *outputs = *o;

	
	// get the dims first
	int *dims = (int*) cust_alloc (sizeof (int) * oo[0]);
	H5LTread_dataset (file_id, "/outputs/dims", H5T_NATIVE_INT, dims);

	int *pos = (int*) cust_alloc (sizeof (int) * oo[0] * 3); // the positions
	H5LTread_dataset (file_id, "/outputs/position", H5T_NATIVE_INT, pos); 

	int *tt = (int*) cust_alloc (sizeof (int) * oo[0]); // the dimensions of the various outputs 
	H5LTread_dataset (file_id, "/outputs/tt", H5T_NATIVE_INT, tt); 

	int *tshots = (int*) cust_alloc (sizeof (int) * oo[0] * oo[1]); // the positions
	H5LTread_dataset (file_id, "/outputs/tshots", H5T_NATIVE_INT, tshots); 

	for (cnt=0 ; cnt < oo[0] ; cnt++) {
		outputs[cnt].counter = 0; // set the counter to 0
		outputs[cnt].dims = dims[cnt];
		memcpy (outputs[cnt].pos, &(pos[cnt*3]), sizeof (int)*3);
		outputs[cnt].tt = tt[cnt];
		if (tshots[cnt*oo[1]] == -1) { // special case
			outputs[cnt].tshots = (int*) cust_alloc (sizeof (int) * 1);
			outputs[cnt].tshots[0] = -1;
		}
		else { // normal case
		outputs[cnt].tshots = (int*) cust_alloc (sizeof (int) * outputs[cnt].tt);
		memcpy (outputs[cnt].tshots, &(tshots[cnt*oo[1]]), sizeof (int) * outputs[cnt].tt);
		}
		// allocate memory to store the data
		if (outputs[cnt].dims == 4) { // analytic
				outputs[cnt].data[0] = (float*) cust_alloc (sizeof (float) * 1 * outputs[cnt].tt);
		}
		else {
			for (cntA = 0 ; cntA<6 ; cntA++) {
				switch (outputs[cnt].dims) {
					case 0: outputs[cnt].tshot_size = 1; break;
					case 1: outputs[cnt].tshot_size = g.xx; break;
					case 2: outputs[cnt].tshot_size = g.xx * g.yy; break;
					case 3: outputs[cnt].tshot_size = g.xx * g.yy * g.zz; break;
					default: printf ("error in OUTS_INITIALIZE.\n"); break;
				}
				outputs[cnt].data[cntA] = (float*) cust_alloc (sizeof (float) * outputs[cnt].tshot_size * outputs[cnt].tt);
				// printf ("%d, ", tshot_size);
			}
		}
	}

for (cnt = 0 ; cnt < oo[0] ; cnt++) {
		printf ("%d: %d dims, [%d, %d, %d] : (%d) ", cnt, outputs[cnt].dims, outputs[cnt].pos[0], outputs[cnt].pos[1], outputs[cnt].pos[2], outputs[cnt].tt);
		if (outputs[cnt].tshots[0] == -1)
			printf ("%d", outputs[cnt].tshots[0]);
		else 
			for (cntA = 0 ; cntA < outputs[cnt].tt ; cntA++)
				printf ("%d ", outputs[cnt].tshots[cntA]);

		printf("\n");
	}

	printf ("Closing hdf5 file (id: %d)...\n", file_id);
	H5Fclose (file_id); //close the hdf5 file

	return oo[0];
}

// puts some strings together
void  my_stringcat (char *a, int num, char *b, char *temp) {
	sprintf (temp, "%s%d%s", a, num, b);
}

void my_dimsdet (int d, hsize_t dims[4], int tt, my_grid g) {
	// determine the dimensions to use
	switch (d) {
		case 4: // analytic
		case 0: 
			dims[0] = 1;
			dims[1] = 1;
			dims[2] = 1;
			dims[3] = tt;
			break;
		case 1: 
			dims[0] = g.xx;
			dims[1] = 1;
			dims[2] = 1;
			dims[3] = tt;
			break;
		case 2: 
			dims[0] = g.xx;
			dims[1] = g.yy;
			dims[2] = 1;
			dims[3] = tt;
			break;
		case 3: 
			dims[0] = g.xx;
			dims[1] = g.yy;
			dims[2] = g.zz;
			dims[3] = tt;
			break;
		default: printf ("error in MY_DIMSDET.\n"); break;
		}
	return;
}

void outs_write (my_out *outputs, int oo, char *hdf5_filename, my_grid g) {
	hid_t file_id;
	hid_t group_id;
	hsize_t dims[4];

	char *field_names[6] = {"/Ex", "/Ey", "/Ez", "/Hx", "/Hy", "/Hz"}; // names of the field components
	char *analytic_name = "/Analytic";

	int cnt, cntF;
	char temp[12];

	printf ("Creating hdf5 file %s...\n", hdf5_filename);
	file_id = H5Fcreate (hdf5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); // create the hdf5 file
	printf ("file id: %d\n", file_id);

	// put this in so we know how many outs there are
	hsize_t oo_dims[1] = {1};
	H5LTmake_dataset (file_id, "/oo", 1, oo_dims, H5T_NATIVE_INT, &oo); // transfer to the hdf5 file

	for (cnt=0 ; cnt < oo; cnt++) { // output all the output structures
		printf ("Open/create out group...\n");
		my_stringcat ("/out", cnt, "", temp);
		// printf ("%s\n", temp);
		group_id = H5Gcreate (file_id, temp, 0); // open the group
		
		// determine the dimensions
		my_dimsdet (outputs[cnt].dims, dims, outputs[cnt].tt, g);
		// printf ("dimensions: %d %d %d %d\n", dims[0], dims[1], dims[2], dims[3]);
		
		if (outputs[cnt].dims == 4) { // analytic
			my_stringcat ("/out", cnt, analytic_name, temp);
			H5LTmake_dataset (file_id, temp, 4, dims, H5T_NATIVE_FLOAT, outputs[cnt].data[0]); // transfer to the hdf5 file
			}
		else { // the other outs
			for (cntF=0 ; cntF<6 ; cntF++) {
				my_stringcat ("/out", cnt, field_names[cntF], temp);
				H5LTmake_dataset (file_id, temp, 4, dims, H5T_NATIVE_FLOAT, outputs[cnt].data[cntF]); // transfer to the hdf5 file
			}
		}
		H5Gclose (group_id); // close the group
	}

	printf ("Closing hdf5 file (file id: %d)...\n", file_id);
	H5Fclose (file_id);

	return;
}


void outs_extract (my_out *o, int oo, my_grid g, my_grid d_g, int cntT, dim3 grid, dim3 threads) {
	
	int cnt, cntF;
	int exportit;
	int offset_output, offset_grid;

	for (cnt=0 ; cnt<oo ; cnt++) { // go through all the outputs and see if we want to save anything from disk
		if (o[cnt].tshots[0] == -1) // always output stuff
			exportit = 1;
		else if (o[cnt].tshots[o[cnt].counter] == cntT) // we want to export this time shot
			exportit = 1;
		else 
			exportit = 0;

		// export if we need to
		if (exportit == 1) {

			// now finally transfer the data
			if (o[cnt].dims == 4) {
				// takes a little more work to extract the information this time
				o[cnt].data[0][o[cnt].counter] = analytic_calculate (d_g, grid, threads, o[cnt].pos[0]);
			}
			else { // the normal extracts
				// calculate the offsets
				offset_output = o[cnt].tshot_size * o[cnt].counter;
				offset_grid = o[cnt].pos[0] + o[cnt].pos[1]*g.xx + o[cnt].pos[2]*g.xx*g.yy;
				for (cntF=0 ; cntF<6 ; cntF++) {
					(cudaMemcpy(&(o[cnt].data[cntF][offset_output]), &(d_g.field[cntF][offset_grid]), sizeof (float) * o[cnt].tshot_size, cudaMemcpyDeviceToHost) );
				}
			}
			o[cnt].counter++;
			}
		}
	return;
}
#endif
