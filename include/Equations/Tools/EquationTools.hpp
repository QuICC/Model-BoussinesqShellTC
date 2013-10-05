/** 
 * @file EquationTools.hpp
 * @brief Implementation of equation tools
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONTOOLS_HPP
#define EQUATIONTOOLS_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "Equations/Tools/EquationConditions.hpp"
#include "Equations/Tools/EquationSorters.hpp"

namespace GeoMHDiSCC {

namespace Equations {

/**
 * @brief Implementation of equation tools
 */
namespace Tools {

   /**
    * @brief  
    *
    * @tparam T                  Template for scalar or vector equation
    * @param rEqs                Vector of shared equations
    * @param rPrognosticRange    Range of prognostic equations
    * @param rDiagnosticRange    Range of diagnostic equations
    * @param rTrivialRange       Range of trivial equations
    */
   template <typename T> void sortByType(typename std::vector<T>& rEqs, std::pair<typename std::vector<T>::iterator, typename std::vector<T>::iterator>& rPrognosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rDiagnosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rTrivialRange);

   /**
    * @brief Identify the solver indexes and set it on the equations
    *
    * @param scalarRange   Range of scalar equations of same type
    * @param vectorRange   Range of vector equations of same type
    *
    * \mhdBug This implementation is very BAD! Needs a cleanup as soon as possible
    */
   void identifySolver(const std::pair<std::vector<SharedIScalarEquation>::iterator,std::vector<SharedIScalarEquation>::iterator>& scalarRange, const std::pair<std::vector<SharedIVectorEquation>::iterator,std::vector<SharedIVectorEquation>::iterator>& vectorRange);

   /**
    * @brief Create a small sparse matrix to select matrix block through Kronecker product
    *
    * @param nBlocks Number of blocks (ie. nBlocks x nBlocks matrix)
    * @param row     Row of the required block
    * @param col     Column of the required block
    */
   SparseMatrix makeBlockMatrix(const int nBlocks, const int row, const int col);

//
// Implementation follows
//

   template <typename T> void sortByType(typename std::vector<T>& rEqs, std::pair<typename std::vector<T>::iterator, typename std::vector<T>::iterator>& rPrognosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rDiagnosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rTrivialRange)
   {
      // Sort equations
      std::stable_sort(rEqs.begin(), rEqs.end(), Sorters::EquationType());

      // Initialise empty scalar ranges
      rPrognosticRange = std::make_pair(rEqs.end(), rEqs.end());
      rDiagnosticRange = std::make_pair(rEqs.end(), rEqs.end());
      rTrivialRange = std::make_pair(rEqs.end(), rEqs.end());

      // Determine the ranges for the different types
      typename std::vector<T>::iterator  eqIt;

      eqIt = std::find_if(rEqs.begin(), rEqs.end(), Conditions::IsPrognostic());
      if(eqIt != rEqs.end())
      {
         rPrognosticRange = std::equal_range(rEqs.begin(), rEqs.end(), *eqIt, Sorters::EquationType());
      }

      eqIt = std::find_if(rEqs.begin(), rEqs.end(), Conditions::IsDiagnostic());
      if(eqIt != rEqs.end())
      {
         rDiagnosticRange = std::equal_range(rEqs.begin(), rEqs.end(), *eqIt, Sorters::EquationType());
      }

      eqIt = std::find_if(rEqs.begin(), rEqs.end(), Conditions::IsTrivial());
      if(eqIt != rEqs.end())
      {
         rTrivialRange = std::equal_range(rEqs.begin(), rEqs.end(), *eqIt, Sorters::EquationType());
      }
   }
}
}
}

#endif // EQUATIONTOOLS_HPP
