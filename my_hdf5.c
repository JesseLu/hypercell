// Functions to create, input data into and close an hdf5 file, which will be used to store the simulation data

#ifndef MY_HDF5_H
#define MY_HDF5_H

#define MY_GROUPSIZE_GUESS 6 // just a guess for how large a group will be
#define MY_MAXMEM 100000000 // 1e8 bytes

#include "hdf5.h"
#include "my_tools.h"

// structure for holding the necessary identifiers for our hdf5 file and its objects
typedef struct my_hdf5_ids {
	hid_t file, group, *dset; // identifiers for the file and currently opened group, datasets and dataspaces
	int dnum; // number of dataspaces and datasets in the currently opened group

	int rank; // number of dimensions
	hsize_t *dims; // stores the sizes of the dimensions 
	int *counter; // tell us when we need to extend the dataset

	float *temp_storage; // store stuff here (in host RAM) before we transfer it over to disk
} my_hdf5_ids;

// create an hdf5 file
void hdf5_openfile (my_hdf5_ids *id, char *filename) {
	id->file = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	return;
}

// close the hdf5 file
void hdf5_closefile (my_hdf5_ids *id) {
	herr_t status;
	status = H5Fclose (id->file);
	return;
}

// create a group and its accompanying datasets and dataspaces
// each dataspace will be "chunked" only in the last dimension
void hdf5_opengroup (my_hdf5_ids *id, char *groupname, int dset_num, char **dset_name, int rank0, int xx, int yy, int zz) {
	herr_t status; // not really useful yet
	hid_t dspace; // temporarily holds the dataspace id
	int cnt; // generic counter variable

	// figure out what lastdimsize should be, i.e. how big a 4D 'chunk' should be
	int lastdimsize = MY_MAXMEM / (xx*yy*zz);
	// copy the dimension size information
	id->rank = rank0;
	id->dims = (hsize_t*) cust_alloc (sizeof (hsize_t) * id->rank);
	id->dims[0] = (hsize_t) xx;
	id->dims[1] = (hsize_t) yy;
	id->dims[2] = (hsize_t) zz;
	id->dims[3] = (hsize_t) lastdimsize;

	// create the group
	id->group = H5Gcreate (id->file, groupname, MY_GROUPSIZE_GUESS);

	// allocate space in the id struct to hold all the necessary dataset ids
	id->dnum = dset_num;
	id->dset = (hid_t*) cust_alloc (sizeof (hid_t) * id->dnum);

	// create the parameters used to enable chunking
	hid_t cparms = H5Pcreate (H5P_DATASET_CREATE); 
	status = H5Pset_chunk (cparms, id->rank, id->dims); // set the chunk dims to dims

	// create the maxdims array
	hsize_t *maxdims = (hsize_t*) cust_alloc ( sizeof(hsize_t) * id->rank );
	memcpy (maxdims, id->dims, sizeof (hsize_t) * id->rank);
	maxdims [id->rank-1] = H5S_UNLIMITED;

	// initiate values that keep track of when we need to extend the dataset
	id->counter = (int*) cust_alloc (sizeof (int) * dset_num);

	// create the datasets and dataspaces
	for (cnt=0 ; cnt < id->dnum ; cnt++) {
		// create a dataspace
		dspace = H5Screate_simple (id->rank, id->dims, maxdims);

		// create a dataset containing the dataspace
		id->dset[cnt] = H5Dcreate (id->group, dset_name[cnt], H5T_IEEE_F64BE, dspace, cparms);

		// close the dataspace
		H5Sclose (dspace);
	}
		
	return;
}

// close all datasets and dataspaces within a group, and then close the group itself
void hdf5_closegroup (my_hdf5_ids *id) {
	// close all dataspaces
	int cnt;
	herr_t status;

	for (cnt=0 ; cnt < id->dnum ; cnt++)
		status = H5Dclose (id->dset[cnt]);

	// close the group
	status = H5Gclose (id->group);

	// return the arrays that were used back to memory
	free (id->dims);
	free(id->counter);

	return;
}

// write data to dataset ind
void hdf5_write (my_hdf5_ids *id, float *data, int ind) {
	// hdf5 variables
	hid_t dataspace, filespace;
	herr_t status;

	// construct the offset array
	hsize_t *offset = (hsize_t*) cust_alloc (sizeof (hsize_t) * id->rank);
	offset[id->rank-1] = id->counter[ind];

	// construct the dimension array for the hyperslab
	hsize_t *hs_dims = (hsize_t*) cust_alloc (sizeof (hsize_t) * id->rank);
	memcpy ( hs_dims, id->dims, sizeof (hsize_t) * id->rank);
	hs_dims[id->rank-1] = 1;

	// update the dims
	id->dims[id->rank-1] = id->counter[ind]+1;
	
	// extend the dataset, if needed
	status = H5Dextend (id->dset[ind], id->dims);
	
	// create a new dataspace and get some kind of filespace id (?)
	dataspace = H5Screate_simple (id->rank, hs_dims, NULL); 
    	filespace = H5Dget_space (id->dset[ind]);
	
	// select and write to the appropriate hyperslab
    	status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL, hs_dims, NULL);  
	status = H5Dwrite (id->dset[ind], H5T_NATIVE_FLOAT, dataspace, filespace, H5P_DEFAULT, data);

	H5Sclose (dataspace); // close the dataspace
	
	id->counter[ind]++; // increment counter 

	return;
}

#endif
