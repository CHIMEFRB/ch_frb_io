#include <iostream>
#include "ch_frb_io_internals.hpp"

//
// Minimal C++ wrappers for the libhdf5 C library.
//
// FIXME (minor): There are error paths where a resource doesn't get freed.
// For example, if read_dataset() fails, then the dataset_id returned by H5Dopen()
// probably won't get closed.
//

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


//
// The 'prop_id' argument should be a property list id created with H5Pcreate().
//
// The 'bitshuffle' argument has the following meaning:
//   0 = no compression
//   1 = try to compress, but if plugin fails then just write uncompressed data instead
//   2 = try to compress, but if plugin fails then print a warning and write uncompressed data instead
//   3 = compression mandatory
//
static void set_bitshuffle(const string &name, hid_t prop_id, int bitshuffle)
{
    if (bitshuffle < 0 || bitshuffle > 3)
	throw runtime_error(name + ": bad value of 'bitshuffle' argument (expected 0 <= bitshuffle <= 3)");

    if (bitshuffle == 0)
	return;

    H5Z_filter_t filter_id = 32008;
    vector<unsigned int> cd_values = { 0, 2 };  // trailing "2" means combine with LZ4 compression

    // H5Pset_filter() only records the filter in the property list.  It does not check
    // that the filter exists, and returns success whether or not the plugin is
    // installed -- a missing plugin is not discovered until H5Dcreate() tries to apply
    // it.  So ask H5Zfilter_avail() directly, and pass H5Z_FLAG_OPTIONAL for the
    // non-mandatory levels, which tells HDF5 to drop the filter and write the dataset
    // uncompressed instead of failing to create it.
    bool avail = (H5Zfilter_avail(filter_id) > 0);

    if (!avail && (bitshuffle == 3))
	throw runtime_error(name + ": Fatal: couldn't load bitshuffle plugin (HDF5 filter 32008).  See README for instructions for compiling bitshuffle, and make sure $HDF5_PLUGIN_PATH is set");
    if (!avail && (bitshuffle == 2))
	cerr << (name + ": warning: couldn't load bitshuffle plugin, data will be written uncompressed\n");

    unsigned int flag = (bitshuffle == 3) ? H5Z_FLAG_MANDATORY : H5Z_FLAG_OPTIONAL;
    herr_t status = H5Pset_filter(prop_id, filter_id, flag, cd_values.size(), &cd_values[0]);

    if (status < 0)
	throw runtime_error(name + ": H5Pset_filter() failed for the bitshuffle filter");
}


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


hdf5_group::hdf5_group(const hdf5_file &f, const string &group_name_, bool create)
{
    this->filename = f.filename;
    this->group_name = group_name_;
    this->group_id = create ? H5Gcreate1(f.file_id, group_name.c_str(), 0) : H5Gopen1(f.file_id, group_name.c_str());

    if (group_id < 0)
	throw runtime_error(filename + ": couldn't open group '" + group_name);
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


void hdf5_group::_check_dataset_shape(const string &dataset_name, hid_t dataset_id, const vector<hsize_t> &expected_shape) const
{
    vector<hsize_t> shape;
    this->_get_dataset_shape(dataset_name, dataset_id, shape);

    if (shape.size() != expected_shape.size())
	throw runtime_error(filename + ": dataset '" + dataset_name + "' is a " + to_string(shape.size()) + "-d array, expected " + to_string(expected_shape.size()) + "-d array");

    for (unsigned int i = 0; i < shape.size(); i++) {
	if (shape[i] != expected_shape[i])
	    throw runtime_error(filename + ": dataset '" + dataset_name + "' has shape " + vstr(shape) + ", expected shape " + vstr(expected_shape));
    }
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

    this->_check_dataset_shape(dataset_name, dataset_id, expected_shape);

    int err = H5Dread(dataset_id, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, out);
    if (err < 0)
	throw runtime_error(filename + ": error reading dataset '" + dataset_name + "'");

    if (dataset_id >= 0)
	H5Dclose(dataset_id);
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


//
// References: 
//   https://www.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/C/H5T/h5ex_t_string.c
//   https://www.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/C/H5T/h5ex_t_vlstring.c
//
void hdf5_group::read_string_dataset(const std::string &dataset_name, std::vector<std::string> &data, const std::vector<hsize_t> &expected_shape) const
{
    hid_t memtype = H5Tcopy(H5T_C_S1);
    if (memtype < 0)
	throw runtime_error(filename + ": H5Tcopy() failed?!");

    hid_t dataset_id = H5Dopen1(this->group_id, dataset_name.c_str());
    if (dataset_id < 0)
	throw runtime_error(filename + ": dataset '" + dataset_name + "' not found");

    this->_check_dataset_shape(dataset_name, dataset_id, expected_shape);    

    hid_t space_id = H5Dget_space(dataset_id);
    if (space_id < 0)
	throw runtime_error(filename + ": H5Dget_space() failed on string-valued dataset '" + dataset_name + "'");

    hid_t datatype_id = H5Dget_type(dataset_id);
    if (datatype_id < 0)
	throw runtime_error(filename + ": H5Dget_type() failed on string-valued dataset '" + dataset_name + "'");

    htri_t is_variable = H5Tis_variable_str(datatype_id);
    if (is_variable < 0)
	throw runtime_error(filename + ": H5Tis_variable_str() failed on string-valued dataset '" + dataset_name + "'");

    hsize_t nstrings = prod(expected_shape);
    data = vector<string> (nstrings);
    
    if (is_variable) {
	vector<char *> c_strings(nstrings, nullptr);
    
	herr_t status = H5Tset_size(memtype, H5T_VARIABLE);
	if (status < 0)
	    throw runtime_error(filename + ": H5Tset_size() failed?!");

	status = H5Dread(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &c_strings[0]);
	if (status < 0)
	    throw runtime_error(filename + ": couldn't read string-valued dataset '" + dataset_name + "'");

	for (hsize_t i = 0; i < nstrings; i++)
	    data[i] = std::string(c_strings[i]);
	
	H5Dvlen_reclaim(memtype, space_id, H5P_DEFAULT, &c_strings[0]);
    }
    else {
	size_t size = H5Tget_size(datatype_id);
	if (size == 0)
	    throw runtime_error(filename + ": H5Tget_size() failed on string-valued dataset '" + dataset_name + "'");

	herr_t status = H5Tset_size(memtype, size+1);
	if (status < 0)
	    throw runtime_error(filename + ": ");
	
	vector<char> buf(nstrings * (size+1), 0);
	this->_read_dataset(dataset_name, memtype, reinterpret_cast<void *> (&buf[0]), expected_shape);

	for (hsize_t i = 0; i < nstrings; i++)
	    data[i] = std::string(&buf[i*(size+1)]);
    }

    if (memtype >= 0)
	H5Tclose(memtype);
    if (dataset_id >= 0)
	H5Dclose(dataset_id);
    if (space_id >= 0)
	H5Sclose(space_id);
    if (datatype_id >= 0)
	H5Tclose(datatype_id);
}


//
// References: 
//   https://www.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/C/H5T/h5ex_t_string.c
//   https://www.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/C/H5T/h5ex_t_vlstring.c
//
// FIXME should have a switch to write fixed-length strings.
//
void hdf5_group::write_string_dataset(const string &dataset_name, const vector<string> &data, const vector<hsize_t> &shape)
{
    if (data.size() != prod(shape))
	throw runtime_error(filename + ": write_string_dataset() called with wrong data.size()");

    vector<const char *> cstr_array(data.size());
    for (unsigned int i = 0; i < data.size(); i++)
	cstr_array[i] = data[i].c_str();

    hid_t space_id = H5Screate_simple(shape.size(), &shape[0], NULL);
    if (space_id < 0)
	throw runtime_error(filename + ": couldn't create dataspace for dataset '" + dataset_name + "'?!");

    hid_t datatype_id = H5Tcopy(H5T_C_S1);
    if (datatype_id < 0)
	throw runtime_error(filename + ": couldn't create datatype for dataset '" + dataset_name + "'?!");

    herr_t status = H5Tset_size(datatype_id, H5T_VARIABLE);
    if (status < 0)
	throw runtime_error(filename + ": H5Tset_size() failed when creating dataset '" + dataset_name + "'?!");    

    hid_t proplist_id = H5Pcreate(H5P_DATASET_CREATE);
    if (proplist_id < 0)
	throw runtime_error(filename + ": H5Pcreate() failed when creating dataset '" + dataset_name + "'?!");    

    hid_t dataset_id = H5Dcreate1(group_id, dataset_name.c_str(), datatype_id, space_id, proplist_id);
    if (dataset_id < 0)
	throw runtime_error(filename + ": couldn't create string-valued dataset '" + dataset_name + "'");

    status = H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cstr_array[0]);
    if (status < 0)
	throw runtime_error(filename + ": couldn't write string-valued dataset '" + dataset_name + "'");

    if (dataset_id >= 0)
	H5Dclose(dataset_id);
    if (proplist_id >= 0)
	H5Pclose(proplist_id);
    if (datatype_id >= 0)
	H5Tclose(datatype_id);
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


_hdf5_extendable_dataset::_hdf5_extendable_dataset(const hdf5_group &g, const std::string &dataset_name_, 
						   const std::vector<hsize_t> &chunk_shape, int axis_, hid_t type_, int bitshuffle) :
    filename(g.filename),
    group_name(g.group_name),
    dataset_name(dataset_name_),
    full_name(g.filename + ": " + group_name + "/" + dataset_name),
    curr_shape(chunk_shape),
    axis(axis_),
    type(type_)
{
    int ndim = chunk_shape.size();

    if (ndim < 1)
	throw runtime_error("hdf5_extendable_dataset: attempt to create zero-dimensional dataset");
    if ((axis < 0) || (axis >= ndim))
	throw runtime_error("hdf5_extendable_dataset: axis is out of range");    

    for (int i = 0; i < ndim; i++) {
	if (chunk_shape[i] <= 0)
	    throw runtime_error("hdf5_extendable_dataset: chunk_shape entry is <= 0");
    }

    vector<hsize_t> max_shape = chunk_shape;
    curr_shape[axis] = 0;
    max_shape[axis] = H5S_UNLIMITED;

    hid_t space_id = H5Screate_simple(ndim, &curr_shape[0], &max_shape[0]);
    if (space_id < 0)
	throw runtime_error(full_name + ": H5Screate_simple() failed?!");

    hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
    if (prop_id < 0)
	throw runtime_error(full_name + ": H5Pcreate() failed?!");

    herr_t err = H5Pset_chunk(prop_id, ndim, &chunk_shape[0]);
    if (err < 0)
	throw runtime_error(full_name + ": H5Pset_chunk() failed?!");

    // Note: the h5py source code refers to a "resize glitch" in chunked HDF5 files,
    // and includes the following line of python code after set_chunk():
    //	    dcpl.set_fill_time(h5d.FILL_TIME_ALLOC)  # prevent resize glitch
    //
    // I don't know what's being referred to here, but I'm guessing it was a temporary
    // bug in some version of libhdf5 which got fixed in later versions.  I decided not
    // to worry about it for now, but left this comment in case the issue somehow
    // becomes relevant later.

    set_bitshuffle(full_name, prop_id, bitshuffle);

    this->dataset_id = H5Dcreate2(g.group_id, dataset_name.c_str(), type, space_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
    if (dataset_id < 0)
	throw runtime_error(full_name + ": couldn't create dataset");

    H5Pclose(prop_id);
    H5Sclose(space_id);
}


_hdf5_extendable_dataset::~_hdf5_extendable_dataset()
{
    H5Dclose(dataset_id);
}


void _hdf5_extendable_dataset::write(const void *data, const vector<hsize_t> &shape)
{
    if (shape.size() != curr_shape.size())
	throw runtime_error(full_name + ": array argument to write() has wrong number of dimensions");

    for (int i = 0; i < (int)shape.size(); i++) {
	if ((i == axis) && (shape[i] <= 0))
	    throw runtime_error(full_name + ": array argument to write() has length-zero axis");
	if ((i != axis) && (shape[i] != curr_shape[i]))
	    throw runtime_error(full_name + ": array argument to write() has dimension mismatched to chunk dimension");
    }

    std::vector<hsize_t> offsets(shape.size(), 0);
    std::vector<hsize_t> new_shape = curr_shape;
    new_shape[axis] += shape[axis];
    offsets[axis] = curr_shape[axis];

    herr_t status = H5Dset_extent(dataset_id, &new_shape[0]);
    if (status < 0)
	throw runtime_error(full_name + ": H5Dset_extent() failed?!");

    hid_t file_space_id = H5Dget_space(dataset_id);
    if (file_space_id < 0)
	throw runtime_error(full_name + ": H5Dget_space() failed?!");
	
    status = H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, &offsets[0], NULL, &shape[0], NULL);
    if (status < 0)
	throw runtime_error(full_name + ": H5Sselect_hyperslab() failed?!");	

    hid_t mem_space_id = H5Screate_simple(shape.size(), &shape[0], NULL);
    if (mem_space_id < 0)
	throw runtime_error(full_name + ": H5Screate_simple() failed?!");

    status = H5Dwrite(dataset_id, type, mem_space_id, file_space_id, H5P_DEFAULT, data);
    if (status < 0)
	throw runtime_error(full_name + ": write to extendable dataset failed");

    this->curr_shape = new_shape;

    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);
}


}  // namespace ch_frb_io
