/** \file TestTFTCoupledBidiffusion2DOne.hpp
 *  \brief Implementation of the TFT test equation for 2D bi-diffusion as coupled system: part one (within 3D model)
 */

#ifndef TESTTFTCOUPLEDBIDIFFUSION2DONE_HPP
#define TESTTFTCOUPLEDBIDIFFUSION2DONE_HPP

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
    * \brief Implementation of the TFT test equation for 2D bi-diffusion as coupled system: part one (within 3D model)
    */
   class TestTFTCoupledBidiffusion2DOne: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         TestTFTCoupledBidiffusion2DOne(SharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~TestTFTCoupledBidiffusion2DOne();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id first, const PhysicalNames::Id second);

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

         /**
          * @brief Storage for the name of the second field
          */
         PhysicalNames::Id mSecondName;

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

   /// Typedef for a shared TestTFTCoupledBidiffusion2DOne
   typedef SharedPtrMacro<TestTFTCoupledBidiffusion2DOne> SharedTestTFTCoupledBidiffusion2DOne;

   /**
    * @brief Get the time matrix block
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param k       Wave number k
    */
   void timeBlock(const TestTFTCoupledBidiffusion2DOne& eq, DecoupledZSparse& mat, const MHDFloat k);

   /**
    * @brief Get the linear matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param k       Wave number k
    */
   void linearBlock(const TestTFTCoupledBidiffusion2DOne& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k);

   /**
    * @brief Get the boundary condition matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param k       Wave number k
    */
   void boundaryBlock(const TestTFTCoupledBidiffusion2DOne& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k);

}
}

#endif // TESTTFTCOUPLEDBIDIFFUSION2DONE_HPP
