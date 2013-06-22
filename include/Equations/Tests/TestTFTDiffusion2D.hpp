/** \file TestTFTDiffusion2D.hpp
 *  \brief Implementation of the TFT test equation for 2D diffusion (within 3D model)
 */

#ifndef TESTTFTDIFFUSION2D_HPP
#define TESTTFTDIFFUSION2D_HPP

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
    * @brief Implementation of the TFT test equation for 2D diffusion (within 3D model)
    */
   class TestTFTDiffusion2D: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         TestTFTDiffusion2D(SharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~TestTFTDiffusion2D();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Generic operator row dispatcher
          */
         virtual DecoupledZSparse operatorRow(const OperatorRowId opId, FieldComponents::Spectral::Id comp, const int matIdx) const;

         /**
          * @brief Initialise the spectral equation matrices
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary spBcIds);

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

   /// Typedef for a shared TestTFTDiffusion2D
   typedef SharedPtrMacro<TestTFTDiffusion2D> SharedTestTFTDiffusion2D;

   /**
    * @brief Get the time matrix block
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param k       Wave number k
    */
   void timeBlock(const TestTFTDiffusion2D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const MHDFloat k);

   /**
    * @brief Get the linear matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param k       Wave number k
    */
   void linearBlock(const TestTFTDiffusion2D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k);

   /**
    * @brief Get the boundary condition matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param k       Wave number k
    */
   void boundaryBlock(const TestTFTDiffusion2D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k);

}
}

#endif // TESTTFTDIFFUSION2D_HPP
