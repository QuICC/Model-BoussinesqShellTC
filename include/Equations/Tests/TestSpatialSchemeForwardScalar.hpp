/**
 * @file TestSpatialSchemeForwardScalar.hpp
 * @brief Implementation of a test equation for the spatial schemes with exact known physical space scalar solution
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTSPATIALCHEMEFORWARDSCALAR_HPP
#define TESTSPATIALCHEMEFORWARDSCALAR_HPP

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
    * @brief Implementation of a test equation for the spatial schemes with exact known physical space scalar solution
    */
   class TestSpatialSchemeForwardScalar: public IScalarEquation
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
          *
          * @param pyName     Python script name
          */
         TestSpatialSchemeForwardScalar(const std::string& pyName);

         /**
          * @brief Destructor
          */
         virtual ~TestSpatialSchemeForwardScalar();

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
         void setSolutionType(const TestSpatialSchemeForwardScalar::SolutionTypeId id);

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
         Datatypes::PhysicalScalarType::PointType scalarPoint(const MHDFloat z, const MHDFloat th, const MHDFloat x) const;

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
         Datatypes::PhysicalScalarType::PointType gradientPoint(const MHDFloat z, const MHDFloat th, const MHDFloat x, const FieldComponents::Physical::Id compId) const;

         /**
          * @brief Type of the exact solution
          */
         SolutionTypeId mTypeId;
   };

   /// Typedef for a shared TestSpatialSchemeForwardScalar
   typedef SharedPtrMacro<TestSpatialSchemeForwardScalar> SharedTestSpatialSchemeForwardScalar;

}
}

#endif // TESTSPATIALCHEMEFORWARDSCALAR_HPP
