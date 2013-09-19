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

namespace GeoMHDiSCC {

namespace Equations {

namespace Sorters {

   /**
    * @brief Sorting functor by equation type
    */
   class EquationType
   {
      public:
         /**
          * @brief Sort scalar equations by equation type
          *
          * @param eqA Left scalar equation
          * @param eqB Right scalar equation
          */
         bool operator()(SharedIScalarEquation eqA, SharedIScalarEquation eqB);

         /**
          * @brief Sort vector equations by equation type 
          *
          * @param eqA Left vector equation
          * @param eqB Right vector equation
          */
         bool operator()(SharedIVectorEquation eqA, SharedIVectorEquation eqB);

      private:
         /**
          * @brief Convert enum equation type to integer for scalar equation
          *
          * @param eqA Scalar equation
          */
         int computeEquationType(SharedIScalarEquation eqA);

         /**
          * @brief Convert enum equation type to integer for vector equation
          *
          * @param eqA Vector equation
          */
         int computeEquationType(SharedIVectorEquation eqA);

   };

}

namespace Conditions {

   /**
    * @brief Condition functor for prognostic equation
    */
   struct IsPrognostic
   {
      /**
       * @brief Check scalar equation
       *
       * @param eqA Scalar equation
       */
      bool operator()(SharedIScalarEquation eqA);

      /**
       * @brief Check vector equation
       *
       * @param eqA vector equation
       */
      bool operator()(SharedIVectorEquation eqA);
   };

   /**
    * @brief Condition functor for diagnostic equation
    */
   struct IsDiagnostic
   {
      /**
       * @brief Check scalar equation
       *
       * @param eqA Scalar equation
       */
      bool operator()(SharedIScalarEquation eqA);

      /**
       * @brief Check vector equation
       *
       * @param eqA vector equation
       */
      bool operator()(SharedIVectorEquation eqA);
   };

   /**
    * @brief Condition functor for trivial equation
    */
   struct IsTrivial
   {
      /**
       * @brief Check scalar equation
       *
       * @param eqA Scalar equation
       */
      bool operator()(SharedIScalarEquation eqA);

      /**
       * @brief Check vector equation
       *
       * @param eqA vector equation
       */
      bool operator()(SharedIVectorEquation eqA);
   };

   /**
    * @brief Condition functor for wrapper
    */
   struct IsWrapper
   {
      /**
       * @brief Check scalar equation
       *
       * @param eqA Scalar equation
       */
      bool operator()(SharedIScalarEquation eqA);

      /**
       * @brief Check vector equation
       *
       * @param eqA vector equation
       */
      bool operator()(SharedIVectorEquation eqA);
   };

}

   /**
    * @brief Implementation of equation tools
    */
   class EquationTools
   {
      public:
         /**
          * @brief  
          *
          * @tparam T                  Template for scalar or vector equation
          * @param rEqs                Vector of shared equations
          * @param rPrognosticRange    Range of prognostic equations
          * @param rDiagnosticRange    Range of diagnostic equations
          * @param rTrivialRange       Range of trivial equations
          */
         template <typename T> static void sortByType(typename std::vector<T>& rEqs, std::pair<typename std::vector<T>::iterator, typename std::vector<T>::iterator>& rPrognosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rDiagnosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rTrivialRange);

         /**
          * @brief Identify the solver indexes and set it on the equations
          *
          * @param scalarRange   Range of scalar equations of same type
          * @param vectorRange   Range of vector equations of same type
          *
          * \mhdBug This implementation is very BAD! Needs a cleanup as soon as possible
          */
         static void identifySolver(const std::pair<std::vector<SharedIScalarEquation>::iterator,std::vector<SharedIScalarEquation>::iterator>& scalarRange, const std::pair<std::vector<SharedIVectorEquation>::iterator,std::vector<SharedIVectorEquation>::iterator>& vectorRange);

      protected:

      private:
         /**
          * @brief Constructor
          */
         EquationTools();

         /**
          * @brief Destructor
          */
         ~EquationTools();
   };

   template <typename T> void EquationTools::sortByType(typename std::vector<T>& rEqs, std::pair<typename std::vector<T>::iterator, typename std::vector<T>::iterator>& rPrognosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rDiagnosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rTrivialRange)
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

#endif // EQUATIONTOOLS_HPP
