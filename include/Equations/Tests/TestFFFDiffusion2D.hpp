/**
 * @file TestFFFDiffusion2D.hpp
 * @brief Implementation of the FFF test equation for 2D diffusion (within 3D model) 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTFFFDIFFUSION2D_HPP
#define TESTFFFDIFFUSION2D_HPP

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
    * @brief Implementation of the FFF test equation for 2D diffusion (within 3D model)
    */
   class TestFFFDiffusion2D: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         TestFFFDiffusion2D(SharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~TestFFFDiffusion2D();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Initialise the boundary condition objects
          */
         virtual void createBoundaries(FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief Generic operator row dispatcher
          */
         virtual DecoupledZSparse operatorRow(const OperatorRowId opId, FieldComponents::Spectral::Id comp, const int matIdx) const;

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

   /// Typedef for a shared TestFFFDiffusion2D
   typedef SharedPtrMacro<TestFFFDiffusion2D> SharedTestFFFDiffusion2D;

   /**
    * @brief Get the time matrix block
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @parm eigs     Wave number k
    */
   void timeBlock(const TestFFFDiffusion2D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs);

   /**
    * @brief Get the linear matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @parm eigs     Wave number k
    */
   void linearBlock(const TestFFFDiffusion2D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs);

   /**
    * @brief Get the boundary condition matrix block on given field
    *
    * @param eq      Equation to work on
    * @param fieldId Physical ID of the field
    * @parm eigs     Wave number k
    */
   void boundaryBlock(const TestFFFDiffusion2D& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, std::vector<MHDFloat>& coeffs, std::vector<Boundary::BCIndex>& bcIdx);

}
}

#endif // TESTFFFDIFFUSION2D_HPP