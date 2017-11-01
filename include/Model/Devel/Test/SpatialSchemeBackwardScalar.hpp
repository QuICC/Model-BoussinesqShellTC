/**
 * @file SpatialSchemeBackwardScalar.hpp
 * @brief Implementation of a test equation for the spatial schemes with exact known spectral space scalar solution
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTSPATIALCHEMEBACKWARDSCALAR_HPP
#define TESTSPATIALCHEMEBACKWARDSCALAR_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IScalarEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Test {

   /**
    * @brief Implementation of a test equation for the spatial schemes with exact known spectral space scalar solution
    */
   class SpatialSchemeBackwardScalar: public IScalarEquation
   {
      public:
         /**
          * @brief Enums for the avaialable exact solutions
          */
         enum SolutionTypeId {
            ZERO,
            CONSTANT,
            EXACT,
            FULL
         };

         /**
          * @brief Constructor
          */
         SpatialSchemeBackwardScalar();

         /**
          * @brief Destructor
          */
         virtual ~SpatialSchemeBackwardScalar();

         /**
          * @brief Overloaded init
          */
         virtual void init();

         /**
          * @brief Set the unknown name and requirements 
          *
          * @param name Name to use for unknown
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set the solution type id
          *
          * @param id ID of the exact solution
          */
         void setSolutionType(const SpatialSchemeBackwardScalar::SolutionTypeId id);

         /**
          * @brief Set the exact solution as nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

         /**
          * @brief Compute the error to exact solution for scalar value 
          */
         MHDFloat computeScalarError() const;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

      private:
         /**
          * @brief Compute single point of solution 
          *
          * @param i Fastest index
          * @param j Medium index
          * @param k Slow index
          */
         Datatypes::SpectralScalarType::PointType scalarPoint(const int i, const int j, const int k) const;

         /**
          * @brief Type of the exact solution
          */
         SolutionTypeId mTypeId;
   };

   /// Typedef for a shared SpatialSchemeBackwardScalar
   typedef SharedPtrMacro<SpatialSchemeBackwardScalar> SharedSpatialSchemeBackwardScalar;

}
}
}

#endif // TESTSPATIALCHEMEBACKWARDSCALAR_HPP