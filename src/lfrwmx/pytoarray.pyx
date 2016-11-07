

cdef public api tonumpyarray(double* data, long long size):
    if not (data and size >= 0): raise ValueError
    cdef np.npy_intp dims = size
    #NOTE: it doesn't take ownership of `data`. You must free `data` yourself
    return np.PyArray_SimpleNewFromData(1, &dims, np.NPY_DOUBLE, <void*>data)


from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF

cdef class ArrayWrapper:
    cdef void* data_ptr
    cdef int size

    cdef set_data(self, int size, void* data_ptr):
        """ Set the data of the array
        This cannot be done in the constructor as it must recieve C-level
        arguments.
        Parameters:
        -----------
        size: int
            Length of the array.
        data_ptr: void*
            Pointer to the data            
        """
        self.data_ptr = data_ptr
        self.size = size

    def __array__(self):
        """ Here we use the __array__ method, that is called when numpy
            tries to get an array from the object."""
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(1, shape,
                                               np.NPY_INT, self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)