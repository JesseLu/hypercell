#ifndef MY_GRID_H
#define MY_GRID_H

#include <stdlib.h>
#include <H5LT.h>
#include "my_tools.h"
// #include "my_materials.h"
// this efficiently implements localized auxiliary variables 
#define MY_MAXMATPARAMS 4

// structure that holds the field values
typedef struct my_grid {
	float *field[6]; // field values (Ex, Ey, ..., Hz)
	float *epsilon[3]; // epsilon values for Ex, Ey and Ez
	int *mat[6]; // material types corresponding to those fields
	int *type_index; // store the material type for every material index
	float *params; // holds the parameters for each material 
	float dx, dt; // the constants
	int xx, yy, zz, tt, mm, T, symz;
} my_grid;


// load everything in first, then organize it
void grid_initialize (my_grid *g, char *hdf5_filename) {
	hid_t file_id; // file id for hdf file
	int cnt;
	char *matgrid_h5names[6] = {"/mat/Ex", "/mat/Ey", "/mat/Ez", "/mat/Bx", "/mat/By", "/mat/Bz"}; // names of the datasets
	char *in_h5names[6] = {"/in/Ex", "/in/Ey", "/in/Ez", "/in/Bx", "/in/By", "/in/Bz"}; // names of the datasets
	char *epsilon_h5names[6] = {"/epsilon/Ex", "/epsilon/Ey", "/epsilon/Ez"}; // names of the datasets

	// temporary variables
	int mm[2]; // size of the materials list
	

	printf ("Opening hdf5 file %s...", hdf5_filename);
	file_id = H5Fopen (hdf5_filename, H5F_ACC_RDONLY, H5P_DEFAULT); // open the hdf5 file
	if ( file_id <= 0 )
		error (-1, 0, "error.\n");	
	printf ("file id: %d\n", file_id);

	// first get dx, xx, yy, zz, and dt, tt
	printf ("Reading in info variables...\n");
	H5LTread_dataset (file_id, "/info/dx", H5T_NATIVE_FLOAT, &(g->dx));
	H5LTread_dataset (file_id, "/info/xx", H5T_NATIVE_INT, &(g->xx));
	H5LTread_dataset (file_id, "/info/yy", H5T_NATIVE_INT, &(g->yy));
	H5LTread_dataset (file_id, "/info/zz", H5T_NATIVE_INT, &(g->zz));
	H5LTread_dataset (file_id, "/info/dt", H5T_NATIVE_FLOAT, &(g->dt));
	H5LTread_dataset (file_id, "/info/tt", H5T_NATIVE_INT, &(g->tt));
	H5LTread_dataset (file_id, "/info/T", H5T_NATIVE_INT, &(g->T));
	H5LTread_dataset (file_id, "/info/symz", H5T_NATIVE_INT, &(g->symz));
	H5LTread_dataset (file_id, "/info/mm", H5T_NATIVE_INT, mm);
	g->mm = mm[0];
	printf ("Info variables: (xx, yy, zz, tt, dx, dt, T, symz, mm) = (%d, %d, %d, %d, %f, %f, %d, %d, [%d,%d])\n", g->xx, g->yy, g->zz, g->tt, g->dx, g->dt, g->T, g->symz, g->mm, mm[1]);

	// now get the material and in-field values
	printf ("Reading in material and in-field grids...\n");
	for (cnt=0 ; cnt<6 ; cnt++) {
		g->mat[cnt] = (int*) cust_alloc (sizeof (int) * g->xx * g->yy * g->zz);
		H5LTread_dataset (file_id, matgrid_h5names[cnt], H5T_NATIVE_INT, g->mat[cnt]); // transfer from hdf5 file
		g->field[cnt] = (float *) cust_alloc (sizeof (float) * g->xx * g->yy * g->zz);
		H5LTread_dataset (file_id, in_h5names[cnt], H5T_NATIVE_FLOAT, g->field[cnt]); // transfer from hdf5 file
	}

	// get epsilon values
	for (cnt=0 ; cnt<3 ; cnt++) {
		g->epsilon[cnt] = (float *) cust_alloc (sizeof (float) * g->xx * g->yy * g->zz);
		H5LTread_dataset (file_id, epsilon_h5names[cnt], H5T_NATIVE_FLOAT, g->epsilon[cnt]); // transfer from hdf5 file
	}

	// input the material list
	printf ("Reading in the material list...\n");
	// the type vector first
	g->type_index = (int*) cust_alloc (sizeof (int) * mm[0]); // allocate space for temporary storage of the material list types
	H5LTread_dataset (file_id, "/matlist/type", H5T_NATIVE_INT, g->type_index); // transfer data from hdf5 file into temporary storage
	// the params list for each material
	g->params = (float*) cust_alloc (sizeof (float) * mm[0] * mm[1]); // allocate space for temporary storage of the material list params
	H5LTread_dataset (file_id, "/matlist/params", H5T_NATIVE_FLOAT, g->params); // transfer data from hdf5 file into temporary storage

	printf ("Closing hdf5 file (id: %d)...\n", file_id);
	H5Fclose (file_id); //close the hdf5 file

	return;
}

void grid_write2hdf (my_grid *g, char *hdf5_filename) {
	hid_t file_id;
	hid_t group_id;
	hsize_t dims[3];

	char *out_h5names[6] = {"/out/Ex", "/out/Ey", "/out/Ez", "/out/Bx", "/out/By", "/out/Bz"}; // names of the datasets

	float *tempfield; // gather a field's values here
	int cnt; // counter variables

	printf ("Creating hdf5 file %s...\n", hdf5_filename);
	file_id = H5Fcreate (hdf5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); // create the hdf5 file
	printf ("file id: %d\n", file_id);
	
	// dimensions of the datasets that we want to write to the hdf5 file
	dims[0] = g->xx;
	dims[1] = g->yy;
	dims[2] = g->zz; 

	printf ("Open/create out group...\n");
	group_id = H5Gcreate (file_id, "/out", 0);

	printf ("Transferring to hdf5 file...\n");
	tempfield = (float*) cust_alloc (sizeof(float) * g->xx*g->yy*g->zz); // allocate memory for temporary storage
	for (cnt=0 ; cnt<6 ; cnt++) {
		H5LTmake_dataset (file_id, out_h5names[cnt], 3, dims, H5T_NATIVE_FLOAT, g->field[cnt]); // transfer to the hdf5 file
	}
	printf ("Freeing temporary field storage...\n");
	free (tempfield);

	printf ("Closing hdf5 file (file id: %d)...\n", file_id);
	H5Gclose (group_id);
	H5Fclose (file_id);

	return;
}


void grid_write2hdf_pointout (my_grid *g, float *po, char *hdf5_filename, int tt) {
	hid_t file_id;
	hid_t group_id;
	hsize_t dims[2];

	printf ("Creating hdf5 file %s...\n", hdf5_filename);
	file_id = H5Fcreate (hdf5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); // create the hdf5 file
	printf ("file id: %d\n", file_id);
	
	// dimensions of the datasets that we want to write to the hdf5 file
	//dims[0] = g->tt;
	dims[0] = tt;
	dims[1] = 6;
	printf ("Open/create out group...\n");
	//printf( "testout");
	group_id = H5Gcreate (file_id, "/out", 0);
/*        FILE *fpr;
        float *timeout;
        timeout = (float*) cust_alloc (sizeof (float) *g->T);
        fpr=fopen("./timeout.txt","r");
        int count = 0;
        for (int i = 0; i < g->T; i++) {
                fscanf(fpr,"%lf",&timeout[i]);
                count++;
        }
        totaltime = count;
        fclose(fpr);
*/
	printf ("Transferring to hdf5 file...\n");
	H5LTmake_dataset (file_id, "/out/pointout", 2, dims, H5T_NATIVE_FLOAT, po); // transfer to the hdf5 file

	printf ("Closing hdf5 file (file id: %d)...\n", file_id);
	H5Gclose (group_id);
	H5Fclose (file_id);

	return;
}

void grid_write2hdf_sliceout (my_grid *g, float *sl, char *hdf5_filename) {
	hid_t file_id;
	hid_t group_id;
	hsize_t dims[3];

	printf ("Creating hdf5 file %s...\n", hdf5_filename);
	file_id = H5Fcreate (hdf5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); // create the hdf5 file
	printf ("file id: %d\n", file_id);
	
	// dimensions of the datasets that we want to write to the hdf5 file
	dims[0] = g->xx;
	dims[1] = g->yy;
	dims[2] = g->tt;

	printf ("Open/create out group...\n");
	group_id = H5Gcreate (file_id, "/out", 0);

	printf ("Transferring to hdf5 file...\n");
	H5LTmake_dataset (file_id, "/out/sliceout", 3, dims, H5T_NATIVE_FLOAT, sl); // transfer to the hdf5 file

	printf ("Closing hdf5 file (file id: %d)...\n", file_id);
	H5Gclose (group_id);
	H5Fclose (file_id);

	return;
}

void grid_write4d (my_grid *g, float *go_field[6], char *hdf5_filename, int tt) {
	hid_t file_id;
	hid_t group_id;
	hsize_t dims[4];

	char *out_h5names[6] = {"/gridout/Ex", "/gridout/Ey", "/gridout/Ez", "/gridout/Hx", "/gridout/Hy", "/gridout/Hz"}; // names of the datasets

	int cnt; // counter variables

	printf ("Creating hdf5 file %s...\n", hdf5_filename);
	file_id = H5Fcreate (hdf5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); // create the hdf5 file
	printf ("file id: %d\n", file_id);
	
	// dimensions of the datasets that we want to write to the hdf5 file
	dims[0] = g->xx;
	dims[1] = g->yy;
	dims[2] = g->zz; 
	dims[3] = tt; 

	printf ("Open/create out group...\n");
	group_id = H5Gcreate (file_id, "/gridout", 0);

	printf ("Transferring to hdf5 file...\n");
	for (cnt=0 ; cnt<6 ; cnt++) {
		H5LTmake_dataset (file_id, out_h5names[cnt], 4, dims, H5T_NATIVE_FLOAT, go_field[cnt]); // transfer to the hdf5 file
	}
	printf ("Freeing temporary field storage...\n");

	printf ("Closing hdf5 file (file id: %d)...\n", file_id);
	H5Gclose (group_id);
	H5Fclose (file_id);

	return;
}

#endif
