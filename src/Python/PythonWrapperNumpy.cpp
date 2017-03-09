/**
 * @file PythonWrapperNumpy.cpp
 * @brief Source of Python interpreter wrapper
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */

// Configuration includes
//

// System includes
//

// External includes
//
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>
#include <numpy/npy_3kcompat.h>


// Class include
//
#include "Python/PythonWrapperNumpy.hpp"

// Project includes
//
#include "Python/PythonWrapper.hpp"
#include "Python/PythonConfig.hpp"
#include "Python/PythonCoordinator.hpp"
#include "Exceptions/Exception.hpp"
#include "IoTools/HumanToId.hpp"

namespace QuICC{

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
