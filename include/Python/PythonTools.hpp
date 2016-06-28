/**
 * @file PythonTools.hpp
 * @brief Tools to work with Python objects
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PYTHONTOOLS_HPP
#define PYTHONTOOLS_HPP

// First include
//
#include "Python/PythonHeader.hpp"

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Tools to work with Python objects
    */
   class PythonTools
   {
      public:
         /**
          * @brief Make a tuple from integer array
          */
         static PyObject* makeTuple(const ArrayI& val);

         /**
          * @brief Make a tuple from vector of double
          */
         static PyObject* makeTuple(const std::vector<MHDFloat>& val);

         /**
          * @brief Make a list from vector of int
          */
         static PyObject* makeList(const std::vector<int>& val);

         /**
          * @brief Make a list of list from vector of vector of int
          */
         static PyObject* makeList(const std::vector<std::vector<int> >& val);

         /**
          * @brief Make a dictionary
          */
         static PyObject* makeDict(const std::vector<std::string>& key, const std::vector<MHDFloat>& val);

         /**
          * @brief Make a dictionary
          */
         static PyObject* makeDict(const std::map<std::string,int>& map);

         /**
          * @brief Make a dictionary
          */
         static PyObject* makeDict(const std::map<NonDimensional::Id,MHDFloat>& map);

         /**
          * @brief Get data from list
          */
         static void getList(std::vector<bool> &rList, PyObject *pList);

         /**
          * @brief Get data from list
          */
         static void getList(std::vector<NonDimensional::Id> &rList, PyObject *pList);

         /**
          * @brief Get data from list
          */
         static void getList(std::vector<PhysicalNames::Id> &rList, PyObject *pList);

         /**
          * @brief Get data from list
          */
         static void getList(std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> > &rList, PyObject *pList);

         /**
          * @brief Get data from dict
          */
         static void getDict(std::map<NonDimensional::Id,MHDFloat> &rMap, PyObject *pDict, const bool replace);

         /**
          * @brief Fill sparse matrix with data from Python call
          */
         static void fillMatrix(SparseMatrix& rMatrix, PyObject* pPyMat);

         /**
          * @brief Fill sparse matrix with data from Python call
          */
         static void fillMatrix(DecoupledZSparse& rMatrix, PyObject* pPyMat);
  
         /**
          * @brief Cleanup wrapper without finalize
          */
         static void cleanup();
  
         /**
          * @brief Finalise the Python interpreter
          */
         static void finalize();
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         PythonTools();

         /**
          * @brief Destructor
          */
         ~PythonTools();
   };

}

#endif // PYTHONTOOLS_HPP
