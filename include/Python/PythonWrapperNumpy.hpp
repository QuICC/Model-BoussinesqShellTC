
// First include
//
#include "Python/PythonHeader.hpp"
// External include

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
}
