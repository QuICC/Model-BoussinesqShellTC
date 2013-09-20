/**
 * @file TestWFTForwardScalar.hpp
 * @brief Implementation of a test equation for the WFT scheme with exact known physical space scalar solution
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTWFTFORWARDSCALAR_HPP
#define TESTWFTFORWARDSCALAR_HPP

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

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of a test equation for the WFT scheme with exact known physical space scalar solution
    */
   class TestWFTForwardScalar: public IScalarEquation
   {
      public:
         /**
          * @brief Enums for the avaialable exact solutions
          */
         enum SolutionTypeId {
            ZERO,
            CONSTANT,
            EXACT
         };

         /**
          * @brief Constructor
          */
         TestWFTForwardScalar();

         /**
          * @brief Destructor
          */
         virtual ~TestWFTForwardScalar();

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
         void setSolutionType(const TestWFTForwardScalar::SolutionTypeId id);

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

         /**
          * @brief Compute the error to exact solution for gradient value 
          */
         MHDFloat computeGradientError() const;

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
          * @param z    Z coordinate
          * @param th   Theta coordinate
          * @param x    X coordinate
          *
          * @return 
          */
         MHDFloat scalarPoint(const MHDFloat z, const MHDFloat th, const MHDFloat x) const;

         /**
          * @brief Compute single point of gradient of solution 
          *
          * @param z       Z coordinate
          * @param th      Theta coordinate
          * @param x       X coordinate
          * @param compId  ID of component
          *
          * @return 
          */
         MHDFloat gradientPoint(const MHDFloat z, const MHDFloat th, const MHDFloat x, const FieldComponents::Physical::Id compId) const;

         /**
          * @brief Type of the exact solution
          */
         SolutionTypeId mTypeId;
   };

   /// Typedef for a shared TestWFTForwardScalar
   typedef SharedPtrMacro<TestWFTForwardScalar> SharedTestWFTForwardScalar;

}
}

#endif // TESTWFTFORWARDSCALAR_HPP
