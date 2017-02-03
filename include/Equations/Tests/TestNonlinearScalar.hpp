/**
 * @file TestNonlinearScalar.hpp
 * @brief Implementation of a scalar nonlinear test equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTNONLINEARSCALAR_HPP
#define TESTNONLINEARSCALAR_HPP

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
#include "Equations/IScalarEquation.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of a scalar nonlinear test equation
    */
   class TestNonlinearScalar: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         TestNonlinearScalar(SharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~TestNonlinearScalar();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);
         
         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

         /**
          * @brief Compute the source term
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual Datatypes::SpectralScalarType::PointType sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

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
   };

   /// Typedef for a shared TestNonlinearScalar
   typedef SharedPtrMacro<TestNonlinearScalar> SharedTestNonlinearScalar;

}
}

#endif // TESTNONLINEARSCALAR_HPP
