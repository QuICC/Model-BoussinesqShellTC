/**
 * @file TestTFTBidiffusion3D.hpp
 * @brief Implementation of the TFT test equation for 3D bi-diffusion (bilaplacian) 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTTFTBIDIFFUSION3D_HPP
#define TESTTFTBIDIFFUSION3D_HPP

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
    * @brief Implementation of the TFT test equation for 3D bi-diffusion (bilaplacian)
    */
   class TestTFTBidiffusion3D: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         TestTFTBidiffusion3D(SharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~TestTFTBidiffusion3D();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

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

   /// Typedef for a shared TestTFTBidiffusion3D
   typedef SharedPtrMacro<TestTFTBidiffusion3D> SharedTestTFTBidiffusion3D;

   /**
    * @brief Get the time matrix block
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @parm eigs     Wave number k
    */
   void timeBlock(const TestTFTBidiffusion3D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::vector<MHDFloat>& eigs);

   /**
    * @brief Get the linear matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @parm eigs     Wave number k
    */
   void linearBlock(const TestTFTBidiffusion3D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs);

   /**
    * @brief Get the boundary condition matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @parm eigs     Wave number k
    */
   void boundaryBlock(const TestTFTBidiffusion3D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs);

}
}

#endif // TESTTFTBIDIFFUSION3D_HPP
