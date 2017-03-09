
// First include
//
#include "Python/PythonHeader.hpp"
// External include

#include <numpy/ndarrayobject.h>
#include <numpy/npy_3kcompat.h>

#include "Python/PythonWrapper.hpp"

namespace QuICC{

class PythonWrapperNumpy: public PythonWrapper{

public:

	static void init();

    /*
     * @brief Fill a full matrix (Eigen::MatrixXd) with data from PyObject
     */
    static void getMatrix(Matrix& rMatrix, PyObject* pMat);

    /*
     * @brief Fill a vector (Eigen::VectorXd) with data from PyObject
     */
    static void getVector(Array& rArray, PyObject* pVec);

    /*
     * @brief Fill a complex vector (Eigen::VectorXcd) with data from PyObjct
     */
    static void getVector(ArrayZ& rArrayZ, PyObject* pVec);

protected:
    /**
     * @brief Constructor
     */
    PythonWrapperNumpy();

    /**
     * @brief Destructor
     */
    ~PythonWrapperNumpy();

private:
    /**
     * @brief Python module object
     */
    static PyObject* mpModule;

    /**
     * @brief Python function object
     */
    static PyObject* mpFunc;

    /**
     * @brief Python class object
     */
    static PyObject* mpClass;

    /**
     * @brief Python class method object
     */
    static PyObject* mpMethod;
};

PyObject* PythonWrapperNumpy::mpModule = NULL;

PyObject* PythonWrapperNumpy::mpFunc = NULL;

PyObject* PythonWrapperNumpy::mpClass = NULL;

PyObject* PythonWrapperNumpy::mpMethod = NULL;

PythonWrapperNumpy::PythonWrapperNumpy():PythonWrapper(){

}

void PythonWrapperNumpy::init(){
	PythonWrapper::init();
}


void PythonWrapperNumpy::getMatrix(Matrix& rMatrix, PyObject* pMat)
{
	   // TODO: some precondition-checking on the number of dimensions and the size
	   // get the size of Python matrix
	   long* dims;
	   dims = PyArray_DIMS(pMat);
	   //dims = ((PyArrayObject_fields *)pMat)->dimensions;

	   // resize the matrix to correct size and assign the pointer
	   // note that that numpy default storage is RowMajor whereas Eigen 3.3.1 is ColumnMajor
	   //rMatrix = Eigen::Map<Eigen::Matrix<MHDFloat,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >((double*)PyArray_DATA(pMat),dims[0],dims[1]);
	   rMatrix = Eigen::Map<Matrix>((double*)PyArray_DATA(pMat),dims[0],dims[1]);
	   //rMatrix = Eigen::Map<Matrix>( ( double*) ((PyArrayObject_fields *)pMat)->data, dims[0], dims[1]);

}

void PythonWrapperNumpy::getVector(Array& rVector, PyObject* pVec)
{

	   // get the  lenght of the vector
	   int len = PyArray_DIM(pVec,0);

	   // resize the eigen::vector and assign the pointer
	   rVector = Eigen::Map<Eigen::VectorXd>((double*)PyArray_DATA(pVec),len);
}

void PythonWrapperNumpy::getVector(ArrayZ& rVector, PyObject* pVec)
{

	   // get the  lenght of the vector
	   int len = PyArray_DIM(pVec,0);

	   // resize the eigen::vector and assign the pointer
	   rVector = Eigen::Map<Eigen::VectorXcd>((std::complex<double>*)PyArray_DATA(pVec),len);
}


}
