#include "ch_frb_io.hpp"
#include "ch_frb_io_internals.hpp"

//
// Minimal C++ wrappers for the libhdf5 C library.
//
// FIXME (minor): There are error paths where a resource doesn't get freed.
// For example, if read_dataset() fails, then the dataset_id returned by H5Dopen()
// probably won't get closed.
//


// The libhdf5 API changed between version 1.6 and 1.8.  
// This block of macros ensures that the code below will work with either version.
#if H5_VERS_MINOR == 6
# define H5Acreate1 H5Acreate
# define H5Dopen1 H5Dopen
# define H5Dcreate1 H5Dcreate
# define H5Gopen1 H5Gopen
#endif

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


hdf5_file::hdf5_file(const string &filename_, bool write, bool clobber)
{
    this->filename = filename_;

    if (write) {
	// FIXME is there something better in libhdf5?
	if (!clobber && file_exists(filename))
	    throw runtime_error(filename + ": file already exists, and clobber=false was specified when creating file");
	this->file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
	this->file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    if (file_id < 0)
	throw runtime_error(filename + ": couldn't open file");
}


hdf5_file::~hdf5_file()
{
    H5Fclose(file_id);
}


hdf5_group::hdf5_group(const hdf5_file &f, const string &group_name_)
{
    this->filename = f.filename;
    this->group_name = group_name_;
    this->group_id = H5Gopen1(f.file_id, group_name.c_str());

    if (group_id < 0)
	throw runtime_error(filename + ": couldn't open group'" + group_name);
}


hdf5_group::~hdf5_group()
{
    H5Gclose(group_id);
}


void hdf5_group::_get_attribute_shape(const string &attr_name, hid_t attr_id, vector<hsize_t> &shape) const
{
    hid_t space_id = H5Aget_space(attr_id);
    if (space_id < 0)
	throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": get_space() failed?!");

    int ndims = H5Sget_simple_extent_ndims(space_id);
    if (ndims < 0)
	throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": get_ndims() failed?!");

    shape.resize(ndims, 0);

    int err = H5Sget_simple_extent_dims(space_id, &shape[0], NULL);
    if (err < 0)
	throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": get_extent_dims() failed?!");

    H5Sclose(space_id);
}


void hdf5_group::get_attribute_shape(const string &attr_name, vector<hsize_t> &shape) const
{
    hid_t attr_id = H5Aopen_name(this->group_id, attr_name.c_str());
    if (attr_id < 0)
	throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": attribute not found");

    _get_attribute_shape(attr_name, attr_id, shape);
    
    H5Aclose(attr_id);    
}


void hdf5_group::_read_attribute(const string &attr_name, hid_t hdf5_type, void *out, const vector<hsize_t> &expected_shape) const
{
    hid_t attr_id = H5Aopen_name(this->group_id, attr_name.c_str());
    if (attr_id < 0)
	throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": attribute not found");

    vector<hsize_t> shape;
    _get_attribute_shape(attr_name, attr_id, shape);

    if (shape.size() != expected_shape.size())
	throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": attribute shape in file didn't match expected shape");

    for (unsigned int i = 0; i < shape.size(); i++) {
	if (shape[i] != expected_shape[i])
	    throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": attribute shape in file didn't match expected shape");
    }
    
    int err = H5Aread(attr_id, hdf5_type, out);
    if (err < 0)
	throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": attribute read failed?!");

    H5Aclose(attr_id);
}


void hdf5_group::_write_attribute(const string &attr_name, hid_t hdf5_type, const void *data, const vector<hsize_t> &shape)
{
    H5S_class_t space_type = (shape.size() > 0) ? H5S_SIMPLE : H5S_SCALAR;

    hid_t space_id = H5Screate(space_type);
    if (space_id < 0)
	throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": H5Screate() failed?!");

    if (shape.size() > 0) {
	if (H5Sset_extent_simple(space_id, shape.size(), &shape[0], &shape[0]) < 0)
	    throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": set_extent() failed?!");
    }

    hid_t attr_id = H5Acreate1(group_id, attr_name.c_str(), hdf5_type, space_id, H5P_DEFAULT);
    if (attr_id < 0)
	throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": couldn't create attribute");

    if (H5Awrite(attr_id, hdf5_type, data) < 0)
	throw runtime_error(filename + ":" + group_name + "/" + attr_name + ": couldn't write attribute");

    H5Aclose(attr_id);
    H5Sclose(space_id);
}


void hdf5_group::_get_dataset_shape(const string &dataset_name, hid_t dataset_id, vector<hsize_t> &shape) const
{
    hid_t space_id = H5Dget_space(dataset_id);
    if (space_id < 0)
	throw runtime_error(filename + ": couldn't open dataspace in dataset '" + dataset_name + "'?!");

    int ndims = H5Sget_simple_extent_ndims(space_id);
    if (ndims < 0)
	throw runtime_error(filename + ": couldn't get dimensions of dataset '" + dataset_name + "'?!");

    shape.resize(ndims, 0);
    
    int err = H5Sget_simple_extent_dims(space_id, &shape[0], NULL);
    if (err < 0)
	throw runtime_error(filename + ": couldn't get dimensions of dataset '" + dataset_name + "'?!");

    H5Sclose(space_id);
}


ssize_t hdf5_group::get_dataset_slen(const string &dataset_name) const
{
    hid_t dataset_id = H5Dopen(this->group_id, dataset_name.c_str(), H5P_DEFAULT);
    if (dataset_id < 0)
	throw runtime_error(filename + ": dataset '" + dataset_name + "' not found");

    hid_t filetype = H5Dget_type(dataset_id);
    if (filetype < 0)
	throw runtime_error(filename + "/" + dataset_name + ": H5Dget_type() failed");

    // FIXME should add check that dataset is a string type!

    size_t size = H5Tget_size(filetype);
    if (size == 0)
	throw runtime_error(filename + "/" + dataset_name + ": H5Tget_size() failed");

    return size + 1;
}


void hdf5_group::get_dataset_shape(const string &dataset_name, vector<hsize_t> &shape) const
{
    hid_t dataset_id = H5Dopen1(this->group_id, dataset_name.c_str());
    if (dataset_id < 0)
	throw runtime_error(filename + ": dataset '" + dataset_name + "' not found");

    this->_get_dataset_shape(dataset_name, dataset_id, shape);

    H5Dclose(dataset_id);
}


void hdf5_group::_read_dataset(const string &dataset_name, hid_t hdf5_type, void *out, const vector<hsize_t> &expected_shape) const
{    
    hid_t dataset_id = H5Dopen1(this->group_id, dataset_name.c_str());
    if (dataset_id < 0)
	throw runtime_error(filename + ": dataset '" + dataset_name + "' not found");

    vector<hsize_t> shape;
    this->_get_dataset_shape(dataset_name, dataset_id, shape);

    if (shape.size() != expected_shape.size())
	throw runtime_error(filename + ": dataset '" + dataset_name + "' is a " + to_string(shape.size()) + "-d array, expected " + to_string(expected_shape.size()) + "-d array");

    for (unsigned int i = 0; i < shape.size(); i++) {
	if (shape[i] != expected_shape[i])
	    throw runtime_error(filename + ": dataset '" + dataset_name + "' has shape " + vstr(shape) + ", expected shape " + vstr(expected_shape));
    }

    int err = H5Dread(dataset_id, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, out);
    if (err < 0)
	throw runtime_error(filename + ": error reading dataset '" + dataset_name + "'");

    if (dataset_id >= 0)
	H5Dclose(dataset_id);
}


void hdf5_group::read_string_dataset(const string &dataset_name, char *out, const vector<hsize_t> &expected_shape, ssize_t slen) const
{
    if (slen < 2)
	throw runtime_error(filename + "/" + dataset_name + ": expected slen >= 2 in hdf5_group::read_string_dataset()");

    hid_t memtype = H5Tcopy(H5T_C_S1);
    if (memtype < 0)
	throw runtime_error(filename + ": H5Tcopy(H5T_C_S1) failed?!");

    herr_t status = H5Tset_size(memtype, slen);
    if (status < 0)
	throw runtime_error(filename + ": ");

    memset(out, 0, prod(expected_shape) * slen);
    this->_read_dataset(dataset_name, memtype, reinterpret_cast<void *> (out), expected_shape);
}


void hdf5_group::_write_dataset(const string &dataset_name, hid_t hdf5_type, const void *data, const vector<hsize_t> &shape)
{
    hid_t space_id = H5Screate(H5S_SIMPLE);
    if (space_id < 0)
	throw runtime_error(filename + ": couldn't create dataspace for dataset '" + dataset_name + "'?!");

    int ret = H5Sset_extent_simple(space_id, shape.size(), &shape[0], &shape[0]);
    if (ret < 0)
	throw runtime_error(filename + ": couldn't set extents in dataset '" + dataset_name + "'?!");

    hid_t dataset_id = H5Dcreate1(group_id, dataset_name.c_str(), hdf5_type, space_id, H5P_DEFAULT);
    if (dataset_id < 0)
	throw runtime_error(filename + ": couldn't create dataset '" + dataset_name + "'");

    ret = H5Dwrite(dataset_id, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (ret < 0)
	throw runtime_error(filename + ": error writing dataset '" + dataset_name + "'");

    if (dataset_id >= 0)
	H5Dclose(dataset_id);
    if (space_id >= 0)
	H5Sclose(space_id);
}


bool hdf5_group::has_attribute(const string &attr_name) const
{
    hid_t attr_id = H5Aopen_name(this->group_id, attr_name.c_str());

    if (attr_id >= 0) {
        H5Aclose(attr_id);
        return true;
    }

    // FIXME: check that error is due to non-existence of attribute with given name, rather than some other problem.
    return false;
}


bool hdf5_group::has_dataset(const string &dataset_name) const
{
    hid_t dataset_id = H5Dopen1(this->group_id, dataset_name.c_str());

    if (dataset_id >= 0) {
        H5Dclose(dataset_id);
        return true;
    }

    // FIXME: check that error is due to non-existence of dataset with given name, rather than some other problem.
    return false;
}


}  // namespace ch_frb_io