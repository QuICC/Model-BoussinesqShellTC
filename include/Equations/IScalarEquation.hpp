/** \file IScalarEquation.hpp
 *  \brief Base for the implementation of a scalar equation
 */

#ifndef ISCALAREQUATION_HPP
#define ISCALAREQUATION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/IEquation.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Base for the implementation of a scalar equation
    */
   class IScalarEquation: public IEquation
   {
      public:
         /**
          * \brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IScalarEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IScalarEquation();

         /**
          * @brief Set the shared pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         virtual void setUnknown(Datatypes::SharedScalarVariableType spUnknown);

         /**
          * @brief Get the unknown variable
          */
         const Datatypes::ScalarVariableType& unknown() const;

         /**
          * @brief Update unknown from dealised data
          *
          * @param rhs     Dealised input
          * @param compId  ID of the component (allows for a more general implementation)
          */
         void updateDealiasedUnknown(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id compId);

         /**
          * @brief Transfer solver solution to equation unknown
          *
          * @param compId  Component ID
          * @param storage Solver solution
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         void storeSolution(FieldComponents::Spectral::Id compId, const DecoupledZMatrix& storage, const int matIdx, const int start);
         void storeSolution(FieldComponents::Spectral::Id compId, const MatrixZ& storage, const int matIdx, const int start);

         /**
          * @brief Initialise the spectral equation matrices
          *
          * It has only a semi-dummy implementation
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary spBcIds);
         
      protected:
         /**
          * @brief Initialise the spectral equation matrices
          *
          * @param spBcIds   List of boundary condition IDs
          */
         void initSpectralMatrices1DPeriodic(const SharedSimulationBoundary spBcIds);

         /**
          * @brief Set the unknown variable
          */
         Datatypes::ScalarVariableType& rUnknown();

      private:
         /**
          * @brief Set the quasi inverse matrix operator
          *
          * It has only a dummy implementation and should never get called!
          */
         virtual void setQuasiInverse(SparseMatrix &mat) const; //= 0;

         /**
          * @brief Set the explicit linear matrix operator
          */
         virtual void setExplicitLinearBlock(DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const; //= 0;

         /**
          * @brief Storage for the shared scalar variable
          */
         Datatypes::SharedScalarVariableType mspUnknown;
   };

   /// Typedef for a shared IScalarEquation
   typedef SharedPtrMacro<IScalarEquation> SharedIScalarEquation;

   /**
    * @brief Transfer equation unknown to solver
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   void copyUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start);
   void copyUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start);

   /**
    * @brief Add source term
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   void addSource(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start);
   void addSource(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start);

   /**
    * @brief General implementation of linear row for equations with a single periodic dimension
    */
   template <typename TEquation> DecoupledZSparse linearRow1DPeriodic(const TEquation& eq, FieldComponents::Spectral::Id comp, const int matIdx);

   /**
    * @brief General implementation of time row for equations with a single periodic dimension
    */
   template <typename TEquation> DecoupledZSparse timeRow1DPeriodic(const TEquation& eq, FieldComponents::Spectral::Id comp, const int matIdx);

   /**
    * @brief General implementation of boundary row for equations with a single periodic dimension
    */
   template <typename TEquation> DecoupledZSparse boundaryRow1DPeriodic(const TEquation& eq, FieldComponents::Spectral::Id comp, const int matIdx);

   /**
    * @brief General implementation of the boundary block for equations with a single periodic dimensions
    */
   void boundaryBlock1DPeriodic(const IScalarEquation& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int pX, const int pZ, const MHDFloat cX, const MHDFloat cZ);

   template <typename TEquation> DecoupledZSparse linearRow1DPeriodic(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get wave number rescale to box size
      MHDFloat k_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx));

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx)), SparseMatrix(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(eq.couplingInfo(compId).blockN(matIdx), eq.couplingInfo(compId).blockN(matIdx)),SparseMatrix(eq.couplingInfo(compId).blockN(matIdx), eq.couplingInfo(compId).blockN(matIdx)));
      SparseMatrix  tmp(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx));

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = eq.couplingInfo(FieldComponents::Spectral::SCALAR).implicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         SparseMatrix   blockMatrix(eq.couplingInfo(compId).nBlocks(),eq.couplingInfo(compId).nBlocks());
         blockMatrix.insert(eq.couplingInfo(compId).fieldIndex(), colIdx) = 1;

         linearBlock(eq, block, *fIt, k_);
         Eigen::kroneckerProduct(blockMatrix, block.first, tmp);
         matrixRow.first += tmp;
         Eigen::kroneckerProduct(blockMatrix, block.second, tmp);
         matrixRow.second += tmp;

         colIdx++;
      }

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   template <typename TEquation> DecoupledZSparse timeRow1DPeriodic(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get wave number rescale to box size
      MHDFloat k_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx));

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx)), SparseMatrix(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(eq.couplingInfo(compId).blockN(matIdx), eq.couplingInfo(compId).blockN(matIdx)),SparseMatrix(eq.couplingInfo(compId).blockN(matIdx), eq.couplingInfo(compId).blockN(matIdx)));
      SparseMatrix  tmp(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx));

      // Create time row
      SparseMatrix   blockMatrix(eq.couplingInfo(compId).nBlocks(),eq.couplingInfo(compId).nBlocks());
      blockMatrix.insert(eq.couplingInfo(compId).fieldIndex(), eq.couplingInfo(compId).fieldIndex()) = 1;
      timeBlock(eq, block, k_);
      Eigen::kroneckerProduct(blockMatrix, block.first, tmp);
      matrixRow.first += tmp;
      Eigen::kroneckerProduct(blockMatrix, block.second, tmp);
      matrixRow.second += tmp;

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   template <typename TEquation> DecoupledZSparse boundaryRow1DPeriodic(const TEquation& eq, FieldComponents::Spectral::Id comp, const int matIdx)
   {
      // Get wave number rescale to box size
      MHDFloat k_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx));

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(eq.couplingInfo(comp).systemN(matIdx), eq.couplingInfo(comp).systemN(matIdx)), SparseMatrix(eq.couplingInfo(comp).systemN(matIdx), eq.couplingInfo(comp).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(eq.couplingInfo(comp).blockN(matIdx), eq.couplingInfo(comp).blockN(matIdx)),SparseMatrix(eq.couplingInfo(comp).blockN(matIdx), eq.couplingInfo(comp).blockN(matIdx)));
      SparseMatrix  tmp(eq.couplingInfo(comp).systemN(matIdx), eq.couplingInfo(comp).systemN(matIdx));

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = eq.couplingInfo(FieldComponents::Spectral::SCALAR).implicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         SparseMatrix   blockMatrix(eq.couplingInfo(comp).nBlocks(),eq.couplingInfo(comp).nBlocks());
         blockMatrix.insert(eq.couplingInfo(comp).fieldIndex(), colIdx) = 1;

         boundaryBlock(eq, block, *fIt, k_);
         Eigen::kroneckerProduct(blockMatrix, block.first, tmp);
         matrixRow.first += tmp;
         Eigen::kroneckerProduct(blockMatrix, block.second, tmp);
         matrixRow.second += tmp;

         colIdx++;
      }

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

}
}

#endif // ISCALAREQUATION_HPP
